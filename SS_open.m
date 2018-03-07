% Standard
addpath('../ConvexModel/')
clear
%close all

% True optimum
optionu = optimoptions('fmincon','Display','off');
xGuess = [0.09, 0.36, 0.1, 0.25, 0.1, 0.1];

[up_opt] = fmincon(@(x)phiFun(x,CSTRplant(x,xGuess)),[3.7780   10.0301  102.8543],[],[],[],[],...
    [0,0,70],[20,50,120],@(x)conFun(x,CSTRplant(x,xGuess)),optionu);
Xp_opt = CSTRplant(up_opt, xGuess);
phip_opt = phiFun(up_opt,Xp_opt);
conp_opt = conFun(up_opt,Xp_opt);

% Run model 0
[u0_opt] = fmincon(@(x)phiCU(x'),[4, 12, 90],[],[],[],[],...
    [0,0,70],[20,50,120],@(x)deal([g1CU(x'),g2CU(x')],[]),optionu);


% Get phi and g
[phi0_opt, dphi0_opt] =  phiCU(u0_opt');
[con0_opt(1), dcon0_opt(1:3)] = g1CU(u0_opt');
[con0_opt(2), dcon0_opt(4:6)] = g2CU(u0_opt');

% Run plant @u0_opt
[Xp] = CSTRplant(u0_opt, xGuess);

% Get phi and g
phip = phiFun(u0_opt,Xp);
phip0 = phip;
conp0 = conFun(u0_opt,Xp);

du = diag([0.01, 0.01, 0.1]);

for i = 1:3
    u = u0_opt + du(i,:);
    dphip(i) = (phiFun(u, CSTRplant(u,Xp)) - phip)/du(i,i);
    a = (conFun(u, CSTRplant(u,Xp)) - conp0)/du(i,i);
    dconp(i) = a(1);
    dconp(i+3) = a(2);
end

% Get modifiers
K = 0.8;
m0phi = K*(phip - phi0_opt);
m0con = K*(conp0 - con0_opt);

m1phi = K*(dphip - dphi0_opt);
m1con = K*(dconp - dcon0_opt);

% Make modified model
phiMod = @(u)(phiCU(u) + m0phi + m1phi*(u-u0_opt'));
g1Mod = @(u)(g1CU(u) + m0con(1) + m1con(1:3)*(u-u0_opt'));
g2Mod = @(u)(g2CU(u) + m0con(2) + m1con(4:6)*(u-u0_opt'));

% set-up iterations
unsolved = 1;
k = 1;
uGuess = u0_opt;
xGuess = Xp;

while unsolved
    % Run model i
    [ui_opt(k,:)] = fmincon(@(x)phiMod(x'),uGuess,[],[],[],[],...
        [0,0,70],[20,50,120],@(x)deal([g1Mod(x'),g2Mod(x')],[]),optionu);
    
    % Get phi and g
    [phii_opt(k), dphii_opt(k,:)] =  phiCU(ui_opt(k,:)');
    [coni_opt(k,1), dconi_opt(k,1:3)] = g1CU(ui_opt(k,:)');
    [coni_opt(k,2), dconi_opt(k,4:6)] = g2CU(ui_opt(k,:)');
    
    % Run plant @u0_opt
    [Xp(k,:)] = CSTRplant(ui_opt(k,:), xGuess);
    
    % Get phi and g
    phip(k) = phiFun(ui_opt(k,:),Xp(k,:));
    conp(k,:) = conFun(ui_opt(k,:),Xp(k,:));
    
    for i = 1:3
        u = ui_opt(k,:) + du(i,:);
        dphip(k,i) = (phiFun(u, CSTRplant(u,Xp(k,:))) - phip(k))/du(i,i);
        a = (conFun(u, CSTRplant(u,Xp(k,:))) - conp(k,:))/du(i,i);
        dconp(k,i) = a(1);
        dconp(k,i+3) = a(2);
    end
    
    % Get modifiers
    m0phi = (1-K)*m0phi + K*(phip(k) - phii_opt(k));
    m0con = (1-K)*m0con + K*(conp(k,:) - coni_opt(k,:));
    
    m1phi = (1-K)*m1phi + K*(dphip(k,:) - dphii_opt(k,:));
    m1con = (1-K)*m1con + K*(dconp(k,:) - dconi_opt(k,:));
    
    % Make modified model
    phiMod = @(u)(phiCU(u) + m0phi + m1phi*(u-ui_opt(k,:)'));
    g1Mod = @(u)(g1CU(u) + m0con(1) + m1con(1:3)*(u-ui_opt(k,:)'));
    g2Mod = @(u)(g2CU(u) + m0con(2) + m1con(4:6)*(u-ui_opt(k,:)'));
    
    k = k + 1;
    if k > 100
        unsolved = 0;
    end
end

plot([0 numel(phip)],-[phip_opt phip_opt],'b--')
hold on
plot([0:numel(phip)],-[phip0, phip])
