% Closed-loop UU
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
    [0,0,60],[24,60,150],@(u)deal([g1Fun(u,openModel(u, xGuess)),g2Fun(u,openModel(u, xGuess))],[]),optionu);
X0_opt = openModel(u0_opt, xGuess);

% Get phi and g
phi0_opt = phiFun(u0_opt, X0_opt);
g10_opt = g1Fun(u0_opt, X0_opt);
g20_opt = g2Fun(u0_opt, X0_opt);

dphi0_opt = finDiff(@(u)phiFun(u,openModel(u, xGuess)), u0_opt, 0.00001)';
dg10_opt = finDiff(@(u)g1Fun(u,openModel(u, xGuess)), u0_opt, 0.00001)';
dg20_opt = finDiff(@(u)g2Fun(u,openModel(u, xGuess)), u0_opt, 0.00001)';

dy = finDiff(@(u)yFromUX(u,openModel(u, xGuess)),u0_opt, 0.00001)';

% Run plant @u0_opt
rp0 = yFromUX(u0_opt,X0_opt);
up0 = plantController(rp0)';
Xp0 = openPlant(up0, X0_opt);

% Get phi and g
phip0 = phiFun(up0,Xp0);
g1p0 = g1Fun(up0,Xp0);
g2p0 = g2Fun(up0,Xp0);

dr = diag([0.001, 0.01]);

for i = 1:2
    r = rp0 + dr(i,:);
    u = plantController(r)';
    dphip(i) = (phiFun(u, openPlant(u,Xp0)) - phip0)/dr(i,i);
    dg1p(i) = (g1Fun(u, openPlant(u,Xp0)) - g1p0)/dr(i,i);
    dg2p(i) = (g2Fun(u, openPlant(u,Xp0)) - g2p0)/dr(i,i);
end

% Get modifiers
K = 0.2;
m0phi = K*(phip0 - phi0_opt);
m0g1 = K*(g1p0 - g10_opt);
m0g2 = K*(g2p0 - g20_opt);

m1phi = K*(dphip - dphi0_opt*pinv(dy))*dy;
m1g1 = K*(dg1p - dg10_opt*pinv(dy))*dy;
m1g2 = K*(dg2p - dg20_opt*pinv(dy))*dy;

% Make modified model
% % phiMod = @(u)(phiFun(u,openModel(u, xGuess)) + m0phi + m1phi*(u-u0_opt)');
% % conMod = @(u)(conFun(u,openModel(u, xGuess)) + m0con + (m1con*(u-u0_opt)')');

phiMod = @(u)(phiFun(u,openModel(u, xGuess)) + m0phi + m1phi*(u-u0_opt)');
g1Mod = @(u)(g1Fun(u,openModel(u, xGuess)) + m0g1 + m1g1*(u-u0_opt)');
g2Mod = @(u)(g2Fun(u,openModel(u, xGuess)) + m0g2 + m1g2*(u-u0_opt)');

% set-up iterations unsolved = 1;
unsolved = 1;
k = 1;
rGuess = rp0;
xGuess = Xp0;
uGuess = u0_opt;

while unsolved
    % Run model i
    [ui_opt(k,:)] = fmincon(@(u)phiMod(u),[4, 12, 90],[],[],[],[],...
        [0,0,60],[10,60,150],@(u)deal([g1Mod(u),g2Mod(u)],[]),optionu);
    Xi_opt(k,:) = openModel(ui_opt(k,:), xGuess);
    
    % Get phi and g
    
    phii_opt(k) = phiFun(ui_opt(k,:), Xi_opt(k,:));
    g1i_opt(k) = g1Fun(ui_opt(k,:), Xi_opt(k,:));
    g2i_opt(k) = g2Fun(ui_opt(k,:), Xi_opt(k,:));
    
    dphii_opt = finDiff(@(u)phiFun(u,openModel(u, xGuess)), ui_opt(k,:), 0.00001)';
    dg1i_opt = finDiff(@(u)g1Fun(u,openModel(u, xGuess)), ui_opt(k,:), 0.00001)';
    dg2i_opt = finDiff(@(u)g2Fun(u,openModel(u, xGuess)), ui_opt(k,:), 0.00001)';

    dy = finDiff(@(u)yFromUX(u,openModel(u, xGuess)),ui_opt(k,:), 0.00001)';

    % Run plant @ui_opt
    rpi(k,:) = yFromUX(ui_opt(k,:),Xi_opt(k,:));
    upi(k,:) = plantController(rpi(k,:))';
    Xpi(k,:) = openPlant(upi(k,:), Xi_opt(k,:));

    % Get phi and g
    phip(k) = phiFun(upi(k,:),Xpi(k,:));
    g1p(k) = g1Fun(upi(k,:),Xpi(k,:));
    g2p(k) = g2Fun(upi(k,:),Xpi(k,:));
    
    for i = 1:2        
        r = rpi(k,:) + dr(i,:);
        u = plantController(r)';
        dphip(i) = (phiFun(u, openPlant(u,Xpi(k,:))) - phip(k))/dr(i,i);
        dg1p(i) = (g1Fun(u, openPlant(u,Xpi(k,:))) - g1p(k))/dr(i,i);
        dg2p(i) = (g2Fun(u, openPlant(u,Xpi(k,:))) - g2p(k))/dr(i,i);
    end
    
    % Get modifiers
    m0phi = (1-K)*m0phi + K*(phip(k) - phii_opt(k));
    m0g1 = (1-K)*m0g1 + K*(g1p(k) - g1i_opt(k));
    m0g2 = (1-K)*m0g2 + K*(g2p(k) - g2i_opt(k));
    
    m1phi = (1-K)*m1phi + K*(dphip - dphii_opt*pinv(dy))*dy;
    m1g1 = (1-K)*m1g1 + K*(dg1p - dg1i_opt*pinv(dy))*dy;
    m1g2 = (1-K)*m1g2 + K*(dg2p - dg2i_opt*pinv(dy))*dy;

    
    % Make modified model
    phiMod = @(u)(phiFun(u,openModel(u, xGuess)) + m0phi + m1phi*(u-ui_opt(k,:))');
    g1Mod = @(u)(g1Fun(u,openModel(u, xGuess)) + m0g1 + m1g1*(u-ui_opt(k,:))');
    g2Mod = @(u)(g2Fun(u,openModel(u, xGuess)) + m0g2 + m1g2*(u-ui_opt(k,:))');

    k = k + 1;
    if k > 20
        unsolved = 0;
    end
end

plot([rp0(1), rpi(:,1)'],[rp0(2), rpi(:,2)'])

