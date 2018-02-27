% Closed-loop UR
addpath ConvexModel CSTR OtherFunctions PlotFunctions
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
con0_opt = conFun(u0_opt, X0_opt);

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

dr = diag([0.001, 0.01]);

for i = 1:2
    r = rp0 + dr(i,:);
    u = plantController(r)';
    dphip(i) = (phiFun(u, openPlant(u,Xp0)) - phip0)/dr(i,i);
    a = (conFun(u, openPlant(u,Xp0)) - conp0)/dr(i,i);
    dconp(i,1) = a(1);
    %dconp(i,2) = a(2);
end

% Get modifiers
K = 0.4;
m0phi = K*(phip0 - phi0_opt);
m0con = K*(conp0 - con0_opt);

m1phi = K*(dphip*dy - dphi0_opt*pinv(dy)*dy);
m1con = K*(dconp'*dy - dcon0_opt*pinv(dy)*dy);

% Make modified model
phiMod = @(u)(phiFun(u,openModel(u, xGuess)) + m0phi + m1phi*(u-u0_opt)');
conMod = @(u)(conFun(u,openModel(u, xGuess)) + m0con + (m1con*(u-u0_opt)')');

% set-up iterations unsolved = 1;
unsolved = 1;
k = 1;
rGuess = rp0;
xGuess = Xp0;
uGuess = u0_opt;

while unsolved
    % Run model i
    [ui_opt(k,:)] = fmincon(@(u)phiMod(u),uGuess,[],[],[],[],...
        [0,0,60],[24,60,150],@(u)deal(conMod(u),[]),optionu);
    Xi_opt(k,:) = openModel(ui_opt(k,:), xGuess);
    
    % Get phi and g
    phii_opt(k) = phiFun(ui_opt(k,:),openModel(ui_opt(k,:), xGuess));
    coni_opt(k,:) = conFun(ui_opt(k,:),openModel(ui_opt(k,:), xGuess));
    
    dphii_opt = finDiff(@(u)phiFun(u,openModel(u, xGuess)), ui_opt(k,:), 0.00001)';
    dconi_opt = finDiff(@(u)conFun(u,openModel(u, xGuess)), ui_opt(k,:), 0.00001)';
    
    dy = finDiff(@(u)yFromUX(u,openModel(u, xGuess)),ui_opt(k,:), 0.00001)';
    
    % Run plant @ui_opt
    rpi(k,:) = yFromUX(ui_opt(k,:),Xi_opt(k,:));
    upi(k,:) = plantController(rpi(k,:))';
    Xpi(k,:) = openPlant(upi(k,:), Xi_opt(k,:));
    
    % Get phi and g
    phip(k) = phiFun(upi(k,:),Xpi(k,:));
    conp(k,:) = conFun(upi(k,:),Xpi(k,:));
    
    for i = 1:2        
        r = rpi(k,:) + dr(i,:);
        u = plantController(r)';
        dphip(i) = (phiFun(u, openPlant(u,Xpi(k,:))) - phip(k))/dr(i,i);
        a = (conFun(u, openPlant(u,Xpi(k,:))) - conp(k,:))/dr(i,i);
        dconp(i,1) = a(1);
%        dconp(i,2) = a(2);
    end
    
    % Get modifiers
    m0phi = (1-K)*m0phi + K*(phip(k) - phii_opt(k));
    m0con = (1-K)*m0con + K*(conp(k,:) - coni_opt(k,:));
    
    m1phi = (1-K)*m1phi + K*(dphip*dy - dphii_opt*pinv(dy)*dy);
    m1con = (1-K)*m1con + K*(dconp'*dy - dconi_opt*pinv(dy)*dy);
    
    % Make modified model
    phiMod = @(u)(phiFun(u,openModel(u, xGuess)) + m0phi + m1phi*(u-ui_opt(k,:))');
    conMod = @(u)(conFun(u,openModel(u, xGuess)) + m0con + (m1con*(u-ui_opt(k,:))')');
    
    k = k + 1;
    if k > 20
        unsolved = 0;
    end
end

plot([rp0(1), rpi(:,1)'],[rp0(2), rpi(:,2)'])

