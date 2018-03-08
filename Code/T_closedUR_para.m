% Transient Closed-loop UR Better Model 
%addpath ConvexModel CSTR OtherFunctions PlotFunctions
clearvars -except fig
%close all

% variables
Kp = -1000;
T0 = 110;

tau = 500;
tFinal = 15000;
kMax = ceil(tFinal/tau);

K = 0.05;
NE = 1;

dPara = [1,1,1] + 0*[0.2,0.2,-0.2]/100;

% True optimum
optionu = optimoptions('fmincon','Display','off');
xGuess = [0.09, 0.36, 0.1, 0.25, 0.1, 0.1];

% @[T0 = 120, Kp = -1000]
rp_opt   = [0.10605,6.05325];
up_opt   = plantController2(rp_opt,xGuess,Kp,T0)';
[~,a]    = ode15s(@(t,y)closedPlantODE(t,y,Kp), [0 100000],[up_opt, xGuess]);
up_opt   = a(end,1:3);
Xp_opt   = a(end,4:end);
phip_opt = phiFun(up_opt,Xp_opt);
g1p_opt  = g1Fun(up_opt,Xp_opt);
g2p_opt  = g2Fun(up_opt,Xp_opt);

% Run model 0
[u0_opt] = fmincon(@(u)phiFun(u,openModelPara(u, dPara, xGuess)),[4, 12, 90],...
    [],[],[],[],[0,0,60],[24,60,150],...
    @(u)deal([g1Fun(u,openModelPara(u, dPara, xGuess)),g2Fun(u,openModelPara(u, dPara, xGuess))],[]),optionu);
X0_opt = openModelPara(u0_opt, dPara, xGuess);

% Get phi and g
phi0_opt = phiFun(u0_opt, X0_opt);
g10_opt = g1Fun(u0_opt, X0_opt);
g20_opt = g2Fun(u0_opt, X0_opt);

dphi0_opt = finDiff(@(u)phiFun(u,openModelPara(u, dPara, xGuess)), u0_opt, 0.00001)';
dg10_opt = finDiff(@(u)g1Fun(u,openModelPara(u, dPara, xGuess)), u0_opt, 0.00001)';
dg20_opt = finDiff(@(u)g2Fun(u,openModelPara(u, dPara, xGuess)), u0_opt, 0.00001)';

dy = finDiff(@(u)yFromUX(u,openModelPara(u, dPara, xGuess)),u0_opt, 0.00001)';

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

if NE ~= 1
    for i = 1:2
        r = rp0 + dr(i,:);
        u = plantController2(r,Xp0,Kp,T0)';
        [c,a] = ode15s(@(t,y)closedPlantODE(t,y,Kp), [0 tau],[u, Xp0]);
        Xp2(i,:) = a(end,:);
        
        dphip(i) = (phiFun(u, a(end,4:end)) - base.phip(end))/dr(i,i);
        dg1p(i) = (g1Fun(u, a(end,4:end)) - base.g1p(end))/dr(i,i);
        dg2p(i) = (g2Fun(u, a(end,4:end)) - base.g2p(end))/dr(i,i);
    end
else %run NE
    dOpt.du = base.u(end,:) - u0_opt;
    dOpt.dC = base.Xp(end,:) - X0_opt;
    dfun = NEgradPara(base.u(end,:),dPara,dOpt);
    
    dphip = (dfun.dphidu' + dphi0_opt)*pinv(dy);
    dg1p = (dfun.dg1du' + dg10_opt)*pinv(dy);
    dg2p = (dfun.dg2du' + dg20_opt)*pinv(dy);
end

% Get modifiers
m0phi = K*(base.phip(end) - phi0_opt);
m0g1 = K*(base.g1p(end) - g10_opt);
m0g2 = K*(base.g2p(end) - g20_opt);

m1phi = K*(dphip - dphi0_opt*pinv(dy));
m1g1 = K*(dg1p - dg10_opt*pinv(dy));
m1g2 = K*(dg2p - dg20_opt*pinv(dy));


% Make modified model
phiMod = @(u)(phiFun(u,openModelPara(u, dPara, xGuess)) + m0phi + m1phi*(yFromUX(u,openModelPara(u, dPara, xGuess))-rp0)');
g1Mod = @(u)(g1Fun(u,openModelPara(u, dPara, xGuess)) + m0g1 + m1g1*(yFromUX(u,openModelPara(u, dPara, xGuess))-rp0)');
g2Mod = @(u)(g2Fun(u,openModelPara(u, dPara, xGuess)) + m0g2 + m1g2*(yFromUX(u,openModelPara(u, dPara, xGuess))-rp0)');

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
    Xi_opt(k,:) = openModelPara(ui_opt(k,:), dPara, xGuess);
    
    % Get phi and g
    phii_opt(k) = phiFun(ui_opt(k,:),openModelPara(ui_opt(k,:), dPara, xGuess));
    g1i_opt(k,:) = g1Fun(ui_opt(k,:),openModelPara(ui_opt(k,:), dPara, xGuess));
    g2i_opt(k,:) = g2Fun(ui_opt(k,:),openModelPara(ui_opt(k,:), dPara, xGuess));

    dphii_opt = finDiff(@(u)phiFun(u,openModelPara(u, dPara, xGuess)), ui_opt(k,:), 0.00001)';
    dg1i_opt = finDiff(@(u)g1Fun(u,openModelPara(u, dPara, xGuess)), ui_opt(k,:), 0.00001)';
    dg2i_opt = finDiff(@(u)g2Fun(u,openModelPara(u, dPara, xGuess)), ui_opt(k,:), 0.00001)';
    
    dy = finDiff(@(u)yFromUX(u,openModelPara(u, dPara, xGuess)),ui_opt(k,:), 0.00001)';
    
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
    
    % Estimate plant gradient
    if NE == 0 %run MU
        for i = 1:2
            r = rpi(k,:) + dr(i,:);
            u = plantController2(r,Xp2(i,4:end),Kp,T0)';
            [c,a] = ode15s(@(t,y)closedPlantODE(t,y,Kp), [0 tau],[u, Xp2(i,4:end)]);
            Xp2(i,:) = a(end,:);
            
            dphip(i) = (phiFun(u, a(end,4:end)) - base.phip(end))/dr(i,i);
            dg1p(i) = (g1Fun(u, a(end,4:end)) - base.g1p(end))/dr(i,i);
            dg2p(i) = (g2Fun(u, a(end,4:end)) - base.g2p(end))/dr(i,i);
        end
    else %run NE
        dOpt.du = base.u(end,:) - u0_opt;
        dOpt.dC = base.Xp(end,:) - X0_opt;
        dfun = NEgradPara(base.u(end,:),dPara,dOpt);
        
        dphip = (dfun.dphidu' + dphi0_opt)*pinv(dy);
        dg1p = (dfun.dg1du' + dg10_opt)*pinv(dy);
        dg2p = (dfun.dg2du' + dg20_opt)*pinv(dy);
    end
    
    % Get modifiers
    m0phi = (1-K)*m0phi + K*(base.phip(end) - phii_opt(k));
    m0g1 = (1-K)*m0g1 + K*(base.g1p(end) - g1i_opt(k));
    m0g2 = (1-K)*m0g2 + K*(base.g2p(end) - g2i_opt(k));
    
    m1phi = (1-K)*m1phi + K*(dphip - dphii_opt*pinv(dy));
    m1g1 = (1-K)*m1g1 + K*(dg1p - dg1i_opt*pinv(dy));
    m1g2 = (1-K)*m1g2 + K*(dg2p - dg2i_opt*pinv(dy));
    
    % Make modified model
    phiMod = @(u)(phiFun(u,openModelPara(u, dPara, xGuess)) + m0phi + m1phi*(yFromUX(u,openModelPara(u, dPara, xGuess))-rpi(k,:))');
    g1Mod = @(u)(g1Fun(u,openModelPara(u, dPara, xGuess)) + m0g1 + m1g1*(yFromUX(u,openModelPara(u, dPara, xGuess))-rpi(k,:))');
    g2Mod = @(u)(g2Fun(u,openModelPara(u, dPara, xGuess)) + m0g2 + m1g2*(yFromUX(u,openModelPara(u, dPara, xGuess))-rpi(k,:))');
    
    k = k + 1;
    if k > kMax
        unsolved = 0;
    end
end

% solution analysis
solysis(rpi,base,rp_opt,phip_opt,'UR',tau,K)

%plots
try fig = allPlots(rpi, rp0, base, 'UR',tau, K, fig);
catch
    fig = allPlots(rpi, rp0, base, 'UR',tau, K);
end


