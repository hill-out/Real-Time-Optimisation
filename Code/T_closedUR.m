% Transient Closed-loop UR 
addpath ConvexModel CSTR OtherFunctions PlotFunctions
clear
%close all

% variables
Kp = -1000;
T0 = 120;
tau = 500;
tFinal = 8000;
K = 0.2;

kMax = ceil(tFinal/tau);

% True optimum
optionu = optimoptions('fmincon','Display','off');
xGuess = [0.09, 0.36, 0.1, 0.25, 0.1, 0.1];

[rp_opt] = [0.12, 6.034];
% fmincon(@(x)phiFun(plantController(x)',openPlant(plantController(x)',xGuess)),...
%     [0.09, 12],[],[],[],[],[0,0],[1,50],...
%     @(x)conFun(plantController(x)',openPlant(plantController(x)',xGuess)),optionu);


% Run model 0
[u0_opt] = fmincon(@(u)phiFun(u,openModel(u, xGuess)),[4, 12, 90],...
    [],[],[],[],[0,0,60],[24,60,150],...
    @(u)deal([g1Fun(u,openModel(u, xGuess)),g2Fun(u,openModel(u, xGuess))],[]),optionu);
X0_opt = openModel(u0_opt, xGuess);

% Get phi and g
phi0_opt = phiFun(u0_opt, X0_opt);
g10_opt = g1Fun(u0_opt, X0_opt);
g20_opt = g2Fun(u0_opt, X0_opt);

dphi0_opt = finDiff(@(u)phiFun(u,openModel(u, xGuess)), u0_opt, 0.00001)';
dg10_opt = finDiff(@(u)g1Fun(u,openModel(u, xGuess)), u0_opt, 0.00001)';
dg20_opt = finDiff(@(u)g2Fun(u,openModel(u, xGuess)), u0_opt, 0.00001)';

dy = finDiff(@(u)yFromUX(u,openModel(u, xGuess)),u0_opt, 0.00001)';

% Run plant to steady
rp0 = yFromUX(u0_opt,X0_opt);
up = plantController2(rp0,X0_opt,Kp,T0)';
[~,a] = ode15s(@(t,y)closedPlantODE(t,y,Kp), [0 10000],[up, X0_opt]);
up0 = a(end,1:3);
Xp0 = a(end,4:end);

% Run plant for tau
[t,Xp] = ode15s(@(t,y)closedPlantODE(t,y,Kp), [0 tau],[up0, Xp0]);
base.t = t-t(end);
base.u = Xp(:,1:3);
base.Xp = Xp(:,4:end);

base.phip = phiFun(up0,base.Xp);
base.g1p = g1Fun(up0,base.Xp);
base.g2p = g2Fun(up0,base.Xp);

dr = diag([0.001, 0.01]);

for i = 1:2
    r = rp0 + dr(i,:);
    u = plantController2(r,Xp0,Kp,T0)';
    [c,a] = ode15s(@(t,y)closedPlantODE(t,y,Kp), [0 tau],[u, Xp0]);
    Xp2(i,:) = a(end,:);
    
    dphip(i) = (phiFun(u, a(end,4:end)) - base.phip(end))/dr(i,i);
    dg1p(i) = (g1Fun(u, a(end,4:end)) - base.g1p(end))/dr(i,i);
    dg2p(i) = (g2Fun(u, a(end,4:end)) - base.g2p(end))/dr(i,i);
end

% Get modifiers
K = 0.4;
m0phi = K*(base.phip(end) - phi0_opt);
m0g1 = K*(base.g1p(end) - g10_opt);
m0g2 = K*(base.g2p(end) - g20_opt);

m1phi = K*(dphip - dphi0_opt*pinv(dy));
m1g1 = K*(dg1p - dg10_opt*pinv(dy));
m1g2 = K*(dg2p - dg20_opt*pinv(dy));


% Make modified model
phiMod = @(u)(phiFun(u,openModel(u, xGuess)) + m0phi + m1phi*(yFromUX(u,openModel(u, xGuess))-rp0)');
g1Mod = @(u)(g1Fun(u,openModel(u, xGuess)) + m0g1 + m1g1*(yFromUX(u,openModel(u, xGuess))-rp0)');
g2Mod = @(u)(g2Fun(u,openModel(u, xGuess)) + m0g2 + m1g2*(yFromUX(u,openModel(u, xGuess))-rp0)');

% set-up iterations unsolved = 1;
unsolved = 1;
k = 1;
rGuess = rp0;
xGuess = Xp0;
uGuess = u0_opt;

while unsolved
    % Run model i
    [ui_opt(k,:)] = fmincon(@(u)phiMod(u),uGuess,[],[],[],[],...
        [0,0,60],[24,60,150],@(u)deal([g1Mod(u),g2Mod(u)],[]),optionu);
    Xi_opt(k,:) = openModel(ui_opt(k,:), xGuess);
    
    % Get phi and g
    phii_opt(k) = phiFun(ui_opt(k,:),openModel(ui_opt(k,:), xGuess));
    g1i_opt(k,:) = g1Fun(ui_opt(k,:),openModel(ui_opt(k,:), xGuess));
    g2i_opt(k,:) = g2Fun(ui_opt(k,:),openModel(ui_opt(k,:), xGuess));

    dphii_opt = finDiff(@(u)phiFun(u,openModel(u, xGuess)), ui_opt(k,:), 0.00001)';
    dg1i_opt = finDiff(@(u)g1Fun(u,openModel(u, xGuess)), ui_opt(k,:), 0.00001)';
    dg2i_opt = finDiff(@(u)g2Fun(u,openModel(u, xGuess)), ui_opt(k,:), 0.00001)';
    
    dy = finDiff(@(u)yFromUX(u,openModel(u, xGuess)),ui_opt(k,:), 0.00001)';
    
    % Run plant for tau
    rpi(k,:) = yFromUX(ui_opt(k,:),Xi_opt(k,:));
    newXp = base.Xp(end,:);
    up = plantController2(rpi(k,:),newXp, Kp, T0)';
    [t,Xp] = ode15s(@(t,y)closedPlantODE(t,y,Kp), [0 tau],[up, newXp]);
    n = numel(t);
    base.t(end+1:end+n) = t+base.t(end);
    base.u(end+1:end+n,:) = Xp(:,1:3);
    base.Xp(end+1:end+n,:) = Xp(:,4:end);
    
    % Get phi and g
    base.phip(end+1:end+n) = phiFun(up,Xp(:,4:end));
    base.g1p(end+1:end+n) = g1Fun(up,Xp(:,4:end));
    base.g2p(end+1:end+n) = g2Fun(up,Xp(:,4:end));
    
    for i = 1:2
        r = rpi(k,:) + dr(i,:);
        u = plantController2(r,Xp2(i,4:end),Kp,T0)';
        [c,a] = ode15s(@(t,y)closedPlantODE(t,y,Kp), [0 tau],[u, Xp2(i,4:end)]);
        Xp2(i,:) = a(end,:);
        
        dphip(i) = (phiFun(u, a(end,4:end)) - base.phip(end))/dr(i,i);
        dg1p(i) = (g1Fun(u, a(end,4:end)) - base.g1p(end))/dr(i,i);
        dg2p(i) = (g2Fun(u, a(end,4:end)) - base.g2p(end))/dr(i,i);
    end
    
    % Get modifiers
    m0phi = (1-K)*m0phi + K*(base.phip(end) - phii_opt(k));
    m0g1 = (1-K)*m0g1 + K*(base.g1p(end) - g1i_opt(k));
    m0g2 = (1-K)*m0g2 + K*(base.g2p(end) - g2i_opt(k));
    
    m1phi = (1-K)*m1phi + K*(dphip - dphii_opt*pinv(dy));
    m1g1 = (1-K)*m1g1 + K*(dg1p - dg1i_opt*pinv(dy));
    m1g2 = (1-K)*m1g2 + K*(dg2p - dg2i_opt*pinv(dy));
    
    % Make modified model
    phiMod = @(u)(phiFun(u,openModel(u, xGuess)) + m0phi + m1phi*(yFromUX(u,openModel(u, xGuess))-rpi(k,:))');
    g1Mod = @(u)(g1Fun(u,openModel(u, xGuess)) + m0g1 + m1g1*(yFromUX(u,openModel(u, xGuess))-rpi(k,:))');
    g2Mod = @(u)(g2Fun(u,openModel(u, xGuess)) + m0g2 + m1g2*(yFromUX(u,openModel(u, xGuess))-rpi(k,:))');
    
    k = k + 1;
    if k > kMax
        unsolved = 0;
    end
end

plot([rp0(1), rpi(:,1)'],[rp0(2), rpi(:,2)'])
