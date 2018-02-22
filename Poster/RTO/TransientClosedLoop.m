% Transient Closed-Loop
addpath('../ConvexModel/')
clear
%close all

% plant optimum
optionu = optimoptions('fmincon','Display','off');
xGuess = [0.09, 0.36, 0.1, 0.25, 0.1, 0.1];

[rp_opt] = [0.1061    6.0527];
phip_opt =   -211.5920;
% fmincon(@(x)phiFun(plantController(x)',CSTRplant(plantController(x)',xGuess)),...
%     [0.09, 12],[],[],[],[],[0,0],[1,50],...
%     @(x)conFun(plantController(x)',CSTRplant(plantController(x)',xGuess)),optionu);

% Run model 0
[r0_opt] = fmincon(@(x)phiCR(x'),[0.09, 12],[],[],[],[],...
    [0,0],[1,50],@(x)deal([g1CR(x'),g2CR(x')],[]),optionu);

% Get phi and g
[phi0_opt, dphi0_opt] =  phiCR(r0_opt');
[con0_opt(1), dcon0_opt(1:2)] = g1CR(r0_opt');
[con0_opt(2), dcon0_opt(3:4)] = g2CR(r0_opt');

% Run plant @u0_opt
Kp = -1000;
T0 = 120;
u0 = plantController2(r0_opt,xGuess,Kp,T0)';
[~,a] = ode15s(@(t,y)controlCSTRode(t,y,Kp), [0 10000],[u0, xGuess]);
u0_opt = a(end,1:3);
Xp0 = a(end,4:end);

% Run plant for tau
tau = 500;
[t,Xp] = ode15s(@(t,y)controlCSTRode(t,y,Kp), [0 tau],[u0_opt, Xp0]);
base.t = t-t(end);
base.Xp = Xp(:,4:end);

% Get phi and g
base.phip = phiFun(u0_opt,base.Xp);
base.conp = conFun(u0_opt,base.Xp);

dr = diag([0.001, 0.01]);

for i = 1:2
    r = r0_opt + dr(i,:);
    u = plantController2(r,Xp0,Kp,T0)';
    [c,a] = ode15s(@(t,y)controlCSTRode(t,y,Kp), [0 tau],[u, Xp0]);
    Xp2(i,:) = a(end,:);
    con12 = conFun(u, a(end,4:end));
    
    base.dphip(i) = (phiFun(u, a(end,4:end)) - base.phip(end))/dr(i,i);
    base.dconp(i) = (con12(1) - base.conp(end,1))/dr(i,i);
    base.dconp(i+2) = (con12(2) - base.conp(end,2))/dr(i,i);
end

% Get modifiers
K = 1;
m0phi = K*(base.phip(end) - phi0_opt);
m0con = K*(base.conp(end,:) - con0_opt);

m1phi = K*(base.dphip - dphi0_opt);
m1con = K*(base.dconp - dcon0_opt);

% Make modified model
phiMod = @(r)(phiCR(r) + m0phi + m1phi*(r-r0_opt'));
g1Mod = @(r)(g1CR(r) + m0con(1) + m1con(1:2)*(r-r0_opt'));
g2Mod = @(r)(g2CR(r) + m0con(2) + m1con(3:4)*(r-r0_opt'));

% set-up iterations
unsolved = 1;
k = 1;
rGuess = r0_opt;
xGuess = Xp(end,4:end);

while unsolved
    % Run model i
    [ri_opt(k,:)] = fmincon(@(x)phiMod(x'),rGuess,[],[],[],[],...
        [0,0],[1,50],@(x)deal([g1Mod(x'),g2Mod(x')],[]),optionu);
    
    
    % Get phi and g
    [phii_opt(k), dphii_opt(k,:)] =  phiCR(ri_opt(k,:)');
    [coni_opt(k,1), dconi_opt(k,1:2)] = g1CR(ri_opt(k,:)');
    [coni_opt(k,2), dconi_opt(k,3:4)] = g2CR(ri_opt(k,:)');
    
    % Run plant @u0_opt
    newXp = base.Xp(end,:);
    ui_opt(k,:) = plantController2(ri_opt(k,:),newXp, Kp, T0)';
    [t,Xp] = ode15s(@(t,y)controlCSTRode(t,y,Kp), [0:1:tau],[ui_opt(k,:), newXp]);
    n = numel(t);
    base.t(end+1:end+n) = t+base.t(end);
    base.Xp(end+1:end+n,:) = Xp(:,4:end);
    
    % Get phi and g
    base.phip(end+1:end+n) = phiFun(ui_opt(k,:),Xp(:,4:end));
    base.conp(end+1:end+n,:) = conFun(ui_opt(k,:),Xp(:,4:end));
    
    for i = 1:2
        r = ri_opt(k,:) + dr(i,:);
        u = plantController2(r,Xp2(i,4:end),Kp,T0)';
        [c,a] = ode15s(@(t,y)controlCSTRode(t,y,Kp), [0:1:tau],[u, Xp2(i,4:end)]);
        Xp2(i,:) = a(end,:);
        con12 = conFun(u, a(end,4:end));
        
        base.dphip(i) = (phiFun(u, a(end,4:end)) - base.phip(end))/dr(i,i);
        base.dconp(i) = (con12(1) - base.conp(end,1))/dr(i,i);
        base.dconp(i+2) = (con12(2) - base.conp(end,2))/dr(i,i);
    end
    
    % Get modifiers
    m0phi = (1-K)*m0phi + K*(base.phip(end) - phii_opt(k));
    m0con = (1-K)*m0con + K*(base.conp(end,:) - coni_opt(k,:));
    
    m1phi = (1-K)*m1phi + K*(base.dphip(end,:) - dphii_opt(k,:));
    m1con = (1-K)*m1con + K*(base.dconp(end,:) - dconi_opt(k,:));
    
    % Make modified model
    phiMod = @(r)(phiCR(r) + m0phi + m1phi*(r-ri_opt(k,:)'));
    g1Mod = @(r)(g1CR(r) + m0con(1) + m1con(1:2)*(r-ri_opt(k,:)'));
    g2Mod = @(r)(g2CR(r) + m0con(2) + m1con(3:4)*(r-ri_opt(k,:)'));
    
    
    k = k + 1;
    if k > 20
        unsolved = 0;
    end
end

%plot(base.t,-base.phip)





