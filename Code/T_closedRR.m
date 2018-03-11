% Transient Closed-Loop
%addpath ConvexModel CSTR OtherFunctions PlotFunctions
clearvars -except fig
close all

% variables
Kp = -1000;
T0 = 120;

tau = 500;
tFinal = 10000;
kMax = ceil(tFinal/tau);

K = 0.5;
NE = 1;

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
%[r0_opt] = fmincon(@(x)phiCR(x'),[0.09, 12],[],[],[],[],...
%    [0,0],[1,50],@(x)deal([g1CR(x'),g2CR(x')],[]),optionu);
[u0_opt] = fmincon(@(u)phiFun(u,openModel(u, xGuess)),[4, 12, 90],...
    [],[],[],[],[0,0,60],[24,60,150],...
    @(u)deal([g1Fun(u,openModel(u, xGuess)),g2Fun(u,openModel(u, xGuess))],[]),optionu);
X0_opt = openModel(u0_opt, xGuess);
r0_opt = yFromUX(u0_opt,X0_opt);

% Get phi and g
[phi0_opt, dphi0_opt] =  phiCR(r0_opt');
[g10_opt, dg10_opt] = g1CR(r0_opt');
[g20_opt, dg20_opt] = g2CR(r0_opt');

% Run plant @u0_opt
u0 = plantController2(r0_opt,xGuess,Kp,T0)';
[~,a] = ode15s(@(t,y)closedPlantODE(t,y,Kp), [0 10000],[u0, xGuess]);
u0_opt = a(end,1:3);
Xp0 = a(end,4:end);

% Run plant for tau
[t,Xp] = ode15s(@(t,y)closedPlantODE(t,y,Kp), [0 tau],[u0_opt, Xp0]);
base.t = t-t(end);
base.u = Xp(:,1:3);
base.Xp = Xp(:,4:end);

% Get phi and g
base.phip = phiFun(u0_opt,base.Xp);
base.g1p = g1Fun(u0_opt,base.Xp);
base.g2p = g2Fun(u0_opt,base.Xp);

dr = diag([0.001, 0.01]);

for i = 1:2
    r = r0_opt + dr(i,:);
    u = plantController2(r,Xp0,Kp,T0)';
    [c,a] = ode15s(@(t,y)closedPlantODE(t,y,Kp), [0 tau],[u, Xp0]);
    Xp2(i,:) = a(end,:);
    
    dphip(i) = (phiFun(u, a(end,4:end)) - base.phip(end))/dr(i,i);
    dg1p(i) = (g1Fun(u,a(end,4:end)) - base.g1p(end))/dr(i,i);
    dg2p(i) = (g1Fun(u,a(end,4:end)) - base.g2p(end))/dr(i,i);
end

% Get modifiers
m0phi = K*(base.phip(end) - phi0_opt);
m0g1 = K*(base.g1p(end) - g10_opt);
m0g2 = K*(base.g2p(end) - g20_opt);

m1phi = K*(dphip - dphi0_opt);
m1g1 = K*(dg1p - dg10_opt);
m1g2 = K*(dg2p - dg20_opt);

% Make modified model
phiMod = @(r)(phiCR(r) + m0phi + m1phi*(r-r0_opt'));
g1Mod = @(r)(g1CR(r) + m0g1 + m1g1*(r-r0_opt'));
g2Mod = @(r)(g2CR(r) + m0g2 + m1g2*(r-r0_opt'));

% set-up iterations
unsolved = 1;
k = 1;
rGuess = r0_opt;
xGuess = Xp(end,4:end);

while unsolved
    % Run model i
    [rpi(k,:)] = fmincon(@(x)phiMod(x'),rGuess,[],[],[],[],...
        [0,0],[1,50],@(x)deal([g1Mod(x'),g2Mod(x')],[]),optionu);
    
    
    % Get phi and g
    [phii_opt(k), dphii_opt] =  phiCR(rpi(k,:)');
    [g1i_opt(k), dg1i_opt] = g1CR(rpi(k,:)');
    [g2i_opt(k), dg2i_opt] = g2CR(rpi(k,:)');
    
    % Run plant @u0_opt
    newXp = base.Xp(end,:);
    ui_opt(k,:) = plantController2(rpi(k,:),newXp, Kp, T0)';
    [t,Xp] = ode15s(@(t,y)closedPlantODE(t,y,Kp), [0:1:tau],[ui_opt(k,:), newXp]);
    n = numel(t);
    base.t(end+1:end+n) = t+base.t(end);
    base.u(end+1:end+n,:) = Xp(:,1:3);
    base.Xp(end+1:end+n,:) = Xp(:,4:end);
    
    % Get phi and g
    base.phip(end+1:end+n) = phiFun(ui_opt(k,:),Xp(:,4:end));
    base.g1p(end+1:end+n,:) = g1Fun(ui_opt(k,:),Xp(:,4:end));
    base.g2p(end+1:end+n,:) = g2Fun(ui_opt(k,:),Xp(:,4:end));
    
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
        dfun = NEgrad(base.u(end,:),dOpt);
        
        dphip = (dfun.dphidu' + dphi0_opt)*pinv(dy);
        dg1p = (dfun.dg1du' + dg10_opt)*pinv(dy);
        dg2p = (dfun.dg2du' + dg20_opt)*pinv(dy);
    end
    
    % Get modifiers
    m0phi = (1-K)*m0phi + K*(base.phip(end) - phii_opt(k));
    m0g1 = (1-K)*m0g1 + K*(base.g1p(end) - g1i_opt(k));
    m0g2 = (1-K)*m0g2 + K*(base.g2p(end) - g2i_opt(k));
    
    m1phi = (1-K)*m1phi + K*(dphip - dphii_opt);
    m1g1 = (1-K)*m1g1 + K*(dg1p - dg1i_opt);
    m1g2 = (1-K)*m1g2 + K*(dg2p - dg2i_opt);
    
    % Make modified model
    phiMod = @(r)(phiCR(r) + m0phi + m1phi*(r-rpi(k,:)'));
    g1Mod = @(r)(g1CR(r) + m0g1 + m1g1*(r-rpi(k,:)'));
    g2Mod = @(r)(g2CR(r) + m0g2 + m1g2*(r-rpi(k,:)'));
    
    
    k = k + 1;
    if k > kMax
        unsolved = 0;
    end
end

% solution analysis
solysis(rpi,base,rp_opt,phip_opt,'RR',tau,K)

%plots
try fig = allPlots(rpi, r0_opt, base, 'RR',tau, K, fig);
catch
    fig = allPlots(rpi, r0_opt, base, 'RR',tau, K);
end





