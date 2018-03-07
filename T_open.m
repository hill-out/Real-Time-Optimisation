% Transient
addpath('../ConvexModel/')
%clear
%close all

% Run model 0
optionu = optimoptions('fmincon','Display','off');
[u0_opt] = fmincon(@(x)phiCU(x'),[4, 12, 90],[],[],[],[],...
    [0,0,70],[20,50,120],@(x)deal([g1CU(x'),g2CU(x')],[]),optionu);

% Get phi and g
[phi0_opt, dphi0_opt] =  phiCU(u0_opt');
[con0_opt(1), dcon0_opt(1:3)] = g1CU(u0_opt');
[con0_opt(2), dcon0_opt(4:6)] = g2CU(u0_opt');

% Run plant @u0_opt
xGuess = [0.09, 0.36, 0.1, 0.25, 0.1, 0.1];
[Xp0] = CSTRplant(u0_opt, xGuess);

% Run plant for tau
tau = 30;
[t,Xp] = ode15s(@(t,y)CSTRode(t,y), [0:0.01:tau],[u0_opt, Xp0]);
base.t = t-t(end);
base.Xp = Xp(:,4:end);

% Get phi and g
base.phip = phiFun(u0_opt,base.Xp);
base.conp = conFun(u0_opt,base.Xp);

du = diag([0.001, 0.001, 0.01]);

for i = 1:3
    u = u0_opt + du(i,:);
    [c,a] = ode15s(@(t,y)CSTRode(t,y), [0:0.01:tau],[u, Xp0]);
    Xp2(i,:) = a(end,:);
    con12 = conFun(u, a(end,4:end));
    
    base.dphip(i) = (phiFun(u, a(end,4:end)) - base.phip(end))/du(i,i);
    base.dconp(i) = (con12(1) - base.conp(end,1))/du(i,i);
    base.dconp(i+3) = (con12(2) - base.conp(end,2))/du(i,i);
end

% Get modifiers
K = 0.1;
m0phi = K*(base.phip(end) - phi0_opt);
m0con = K*(base.conp(end,:) - con0_opt);

m1phi = K*(base.dphip - dphi0_opt);
m1con = K*(base.dconp - dcon0_opt);

% Make modified model
phiMod = @(u)(phiCU(u) + m0phi + m1phi*(u-u0_opt'));
g1Mod = @(u)(g1CU(u) + m0con(1) + m1con(1:3)*(u-u0_opt'));
g2Mod = @(u)(g2CU(u) + m0con(2) + m1con(4:6)*(u-u0_opt'));

% set-up iterations
unsolved = 1;
k = 1;
uGuess = u0_opt;
xGuess = Xp(end,4:end);

while unsolved
    % Run model i
    [ui_opt(k,:)] = fmincon(@(x)phiMod(x'),uGuess,[],[],[],[],...
        [0,0,70],[20,50,120],@(x)deal([g1Mod(x'),g2Mod(x')],[]),optionu);
    
    % Get phi and g
    [phii_opt(k), dphii_opt(k,:)] =  phiCU(ui_opt(k,:)');
    [coni_opt(k,1), dconi_opt(k,1:3)] = g1CU(ui_opt(k,:)');
    [coni_opt(k,2), dconi_opt(k,4:6)] = g2CU(ui_opt(k,:)');
    
    % Run plant @u0_opt
    newXp = base.Xp(end,:);
    [t,Xp] = ode15s(@(t,y)CSTRode(t,y), [0:0.1:tau],[ui_opt(k,:), newXp]);
    n = numel(t);
    base.t(end+1:end+n) = t+base.t(end);
    base.Xp(end+1:end+n,:) = Xp(:,4:end);
    
    % Get phi and g
    base.phip(end+1:end+n) = phiFun(ui_opt(k,:),Xp(:,4:end));
    base.conp(end+1:end+n,:) = conFun(ui_opt(k,:),Xp(:,4:end));
    
    for i = 1:3
        u = ui_opt(k,:) + du(i,:);
        [c,a] = ode15s(@(t,y)CSTRode(t,y), [0:0.1:tau],[u, Xp2(i,4:end)]);
        Xp2(i,:) = a(end,:);
        con12 = conFun(u, a(end,4:end));
        
        base.dphip(i) = (phiFun(u, a(end,4:end)) - base.phip(end))/du(i,i);
        base.dconp(i) = (con12(1) - base.conp(end,1))/du(i,i);
        base.dconp(i+3) = (con12(2) - base.conp(end,2))/du(i,i);
    end
    
    % Get modifiers
    m0phi = (1-K)*m0phi + K*(base.phip(end) - phii_opt(k));
    m0con = (1-K)*m0con + K*(base.conp(end,:) - coni_opt(k,:));
    
    m1phi = (1-K)*m1phi + K*(base.dphip(end,:) - dphii_opt(k,:));
    m1con = (1-K)*m1con + K*(base.dconp(end,:) - dconi_opt(k,:));
    
    % Make modified model
    phiMod = @(u)(phiCU(u) + m0phi + m1phi*(u-ui_opt(k,:)'));
    g1Mod = @(u)(g1CU(u) + m0con(1) + m1con(1:3)*(u-ui_opt(k,:)'));
    g2Mod = @(u)(g2CU(u) + m0con(2) + m1con(4:6)*(u-ui_opt(k,:)'));
    
    
    k = k + 1;
    if k > 200
        unsolved = 0;
    end
end


plot(base.t,-base.phip)
hold on






