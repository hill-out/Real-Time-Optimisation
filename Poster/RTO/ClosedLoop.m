% Closed-loop
addpath('../ConvexModel/')
clear
%close all

% True optimum
optionu = optimoptions('fmincon','Display','off');
xGuess = [0.09, 0.36, 0.1, 0.25, 0.1, 0.1];

[rp_opt] = [0.12, 6.034];
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
u0_opt = plantController(r0_opt)';
[Xp] = CSTRplant(u0_opt, xGuess);

% Get phi and g
phip0 = phiFun(u0_opt,Xp);
conp0 = conFun(u0_opt,Xp);

dr = diag([0.0001, 0.01]);

for i = 1:2
    r = r0_opt + dr(i,:);
    u = plantController(r)';
    dphip(i) = (phiFun(u, CSTRplant(u,Xp)) - phip0)/dr(i,i);
    a = (conFun(u, CSTRplant(u,Xp)) - conp0)/dr(i,i);
    dconp(i) = a(1);
    dconp(i+2) = a(2);
end

% Get modifiers
K = 0.4;
m0phi = K*(phip0 - phi0_opt);
m0con = K*(conp0 - con0_opt);

m1phi = K*(dphip - dphi0_opt);
m1con = K*(dconp - dcon0_opt);

% Make modified model
phiMod = @(r)(phiCR(r) + m0phi + m1phi*(r-r0_opt'));
g1Mod = @(r)(g1CR(r) + m0con(1) + m1con(1:2)*(r-r0_opt'));
g2Mod = @(r)(g2CR(r) + m0con(2) + m1con(3:4)*(r-r0_opt'));

% set-up iterations
unsolved = 1;
k = 1;
rGuess = r0_opt;
xGuess = Xp;

while unsolved
    % Run model i
    [ri_opt(k,:)] = fmincon(@(x)phiMod(x'),rGuess,[],[],[],[],...
        [0,0],[1,50],@(x)deal([g1Mod(x'),g2Mod(x')],[]),optionu);
    ui_opt(k,:) = plantController(ri_opt(k,:))';
    
    % Get phi and g
    [phii_opt(k), dphii_opt(k,:)] =  phiCR(ri_opt(k,:)');
    [coni_opt(k,1), dconi_opt(k,1:2)] = g1CR(ri_opt(k,:)');
    [coni_opt(k,2), dconi_opt(k,3:4)] = g2CR(ri_opt(k,:)');
    
    % Run plant @ui_opt
    [Xp(k,:)] = CSTRplant(ui_opt(k,:), xGuess);
    
    % Get phi and g
    phip(k) = phiFun(ui_opt(k,:),Xp(k,:));
    conp(k,:) = conFun(ui_opt(k,:),Xp(k,:));
    
    for i = 1:2
        r = ri_opt(k,:) + dr(i,:);
        u = plantController(r)';
        
        dphip(k,i) = (phiFun(u, CSTRplant(u,Xp(k,:))) - phip(k))/dr(i,i);
        a = (conFun(u, CSTRplant(u,Xp(k,:))) - conp(k,:))/dr(i,i);
        dconp(k,i) = a(1);
        dconp(k,i+2) = a(2);
    end
    
    % Get modifiers
    m0phi = (1-K)*m0phi + K*(phip(k) - phii_opt(k));
    m0con = (1-K)*m0con + K*(conp(k,:) - coni_opt(k,:));
    
    m1phi = (1-K)*m1phi + K*(dphip(k,:) - dphii_opt(k,:));
    m1con = (1-K)*m1con + K*(dconp(k,:) - dconi_opt(k,:));
    
    % Make modified model
    phiMod = @(r)(phiCR(r) + m0phi + m1phi*(r-ri_opt(k,:)'));
    g1Mod = @(r)(g1CR(r) + m0con(1) + m1con(1:2)*(r-ri_opt(k,:)'));
    g2Mod = @(r)(g2CR(r) + m0con(2) + m1con(3:4)*(r-ri_opt(k,:)'));
    
    k = k + 1;
    if k > 400
        unsolved = 0;
    end
end

plot([r0_opt(1), ri_opt(:,1)'],[r0_opt(2), ri_opt(:,2)'])

