% Closed-loop
addpath ConvexModel\ CSTR\ OtherFunctions\ PlotFunctions\
clear
%close all

% True optimum
optionu = optimoptions('fmincon','Display','off');
xGuess = [0.09, 0.36, 0.1, 0.25, 0.1, 0.1];

[rp_opt] = [0.12, 6.034];
% fmincon(@(x)phiFun(plantController(x)',openPlant(plantController(x)',xGuess)),...
%     [0.09, 12],[],[],[],[],[0,0],[1,50],...
%     @(x)conFun(plantController(x)',openPlant(plantController(x)',xGuess)),optionu);


% Run model 0
[u0_opt] = fmincon(@(u)phiFun(u,openModel(u, xGuess)),[4, 12, 90],[],[],[],[],...
    [0,0,60],[24,60,150],@(u)deal(conFun(u,openModel(u, xGuess)),[]),optionu);
X0_opt = openModel(u0_opt, xGuess);

% Get phi and g
phi0_opt = phiFun(u0_opt, X0_opt);
con0_opt = phiFun(u0_opt, X0_opt);

dphi0_opt = finDiff(@(u)phiFun(u,openModel(u, xGuess)), u0_opt, 0.00001)';
dcon0_opt = finDiff(@(u)conFun(u,openModel(u, xGuess)), u0_opt, 0.00001)';

dy = finDiff(@(u)yFromUX(u,openModel(u, xGuess)),u0_opt, 0.00001)';

% Run plant @u0_opt
rp0 = yFromUX(u0_opt,X0_opt);
up0 = plantController(rp0)';
Xp0 = openPlant(up0, X0_opt);

% Get phi and g
phip0 = phiFun(up0,Xp0);
conp0 = conFun(up0,Xp0);

dr = diag([0.00001, 0.001]);

for i = 1:2
    r = rp0 + dr(i,:);
    u = plantController(r)';
    dphip(i) = (phiFun(u, openPlant(u,Xp0)) - phip0)/dr(i,i);
    a = (conFun(u, openPlant(u,Xp0)) - conp0)/dr(i,i);
    dconp(i) = a(1);
    dconp(i+2) = a(2);
end

% Get modifiers
K = 0.4;
m0phi = K*(phip0 - phi0_opt);
m0con = K*(conp0 - con0_opt);

m1phi = K*(dphip - dphi0_opt*pinv(dy));
m1con = K*(dconp - reshape((dcon0_opt*pinv(dy))',1,[]));

% Make modified model
phiMod = @(u)(phiFun(u,openModel(u, xGuess)) + m0phi + m1phi*(yFromUX(u,openModel(u, xGuess))-rp0)');
g1Mod = @(u)(g1CR(r) + m0con(1) + m1con(1:2)*(r-r0_opt'));
g2Mod = @(u)(g2CR(r) + m0con(2) + m1con(3:4)*(r-r0_opt'));

% set-up iterations unsolved = 1;
k = 1;
rGuess = r0_opt;
xGuess = Xp;

while unsolved
    % Run model i
    [u0_opt] = fmincon(@(u)phiFun(u,openModel(u, xGuess)),[4, 12, 90],[],[],[],[],...
        [0,0,60],[24,60,150],@(u)deal(conFun(u,openModel(u, xGuess)),[]),optionu);
    X0_opt = openModel(u0_opt, xGuess);
    
    % Get phi and g
    phi0_opt = phiFun(u0_opt, X0_opt);
    con0_opt = phiFun(u0_opt, X0_opt);
    
    dphi0_opt = finDiff(@(u)phiFun(u,openModel(u, xGuess)), u0_opt, 0.00001)';
    dcon0_opt = finDiff(@(u)conFun(u,openModel(u, xGuess)), u0_opt, 0.00001)';
    
    dy = finDiff(@(u)yFromUX(u,openModel(u, xGuess)),u0_opt, 0.00001)';

    % Run plant @ui_opt
    [Xp(k,:)] = openPlant(ui_opt(k,:), xGuess);
    
    % Get phi and g
    phip(k) = phiFun(ui_opt(k,:),Xp(k,:));
    conp(k,:) = conFun(ui_opt(k,:),Xp(k,:));
    
    for i = 1:2
        r = ri_opt(k,:) + dr(i,:);
        u = plantController(r)';
        
        dphip(k,i) = (phiFun(u, openPlant(u,Xp(k,:))) - phip(k))/dr(i,i);
        a = (conFun(u, openPlant(u,Xp(k,:))) - conp(k,:))/dr(i,i);
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

