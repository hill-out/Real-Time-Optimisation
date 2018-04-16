% Transient Closed-loop UU conv
%addpath ConvexModel CSTR OtherFunctions PlotFunctions
clearvars -except fig
%close all

% variables
Kp = -1000;
T0 = 120;

tFinal = 5000;
steady = 0;
if steady
    tau = 2000;
else
    tau = 30;
end

K = 0.8;

meth = 0;

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
[u0_opt] = fmincon(@(u)phiFun(u,openModel(u, xGuess)),[4, 12, 90],...
    [],[],[],[],[0,0,60],[24,60,150],...
    @(u)deal([g1Fun(u,openModel(u, xGuess)),g2Fun(u,openModel(u, xGuess))],[]),optionu);
X0_opt = openModel(u0_opt, xGuess);

% Get phi and g
[phi0_opt, dphi0_opt] = phiCU(u0_opt');
[g10_opt, dg10_opt] = g1CU(u0_opt');
[g20_opt, dg20_opt] = g2CU(u0_opt');

% Run model 0
dphi0_NE = finDiff(@(u)phiFun(u,openModel(u, xGuess)), u0_opt, 0.0001)';
dg10_NE = finDiff(@(u)g1Fun(u,openModel(u, xGuess)), u0_opt, 0.00001)';
dg20_NE = finDiff(@(u)g2Fun(u,openModel(u, xGuess)), u0_opt, 0.00001)';

dy = finDiff(@(u)yFromUX(u,openModel(u, xGuess)),u0_opt, 0.00001)';

% Run plant to steady
rp0 = yFromUX(u0_opt,X0_opt);
up = plantController2(rp0,X0_opt,Kp,T0)';
[~,a] = ode15s(@(t,y)closedPlantODE(t,y,Kp), [0 10000],[up, X0_opt]);
up0 = a(end,1:3);
Xp0 = a(end,4:end);

% Do FE at start to estimate dudu
du = diag([0.001, 0.001, 0.01]);

for i = 1:3
    u = u0_opt + du(i,:);
    X = openModel(u, xGuess);
    r = yFromUX(u,X);
    up = plantController2(r,Xp0,Kp,T0)';
    [~,a] = ode15s(@(t,y)closedPlantODE(t,y,Kp), [0 5000],[up, Xp0]);
    dudu(i,:) = (a(end,1:3) - up0)/du(i,i);
end

% Run plant for tau
[t,Xp] = ode15s(@(t,y)closedPlantODE(t,y,Kp), [0 tau],[up0, Xp0]);
base.t = t-t(end);
base.u = Xp(:,1:3);
base.Xp = Xp(:,4:end);

base.ti = base.t(end);

base.phip = phiFun(up0,base.Xp);
base.g1p = g1Fun(up0,base.Xp);
base.g2p = g2Fun(up0,base.Xp);

dr = diag([0.001, 0.01]);

if meth == -1
    % Perfect grad. est.
    xSteady = closedPlant(rp0,X0_opt, @(r, x)plantController2(r, x, Kp, T0));
    uSteady = plantController2(rp0, xSteady, Kp, T0)';
    
    for i = 1:2
        r = rp0 + dr(i,:);
        
        xSteady2 = closedPlant(r,xSteady, @(r, x)plantController2(r, x, Kp, T0));
        uSteady2 = plantController2(r, xSteady2, Kp, T0)';
        
        dphip(i) = (phiFun(uSteady2,xSteady2) - phiFun(uSteady,xSteady))/dr(i,i);
        dg1p(i) = (g1Fun(xSteady2,xSteady2) - g1Fun(uSteady,xSteady))/dr(i,i);
        dg2p(i) = (g2Fun(xSteady2,xSteady2) - g2Fun(uSteady,xSteady))/dr(i,i);
    end
    
    dphip = dphip*dy;
    dg1p = dg1p*dy;
    dg2p = dg2p*dy;
elseif meth == 0 %run MU
    for i = 1:2
        r = rp0 + dr(i,:);
        u = plantController2(r,Xp0,Kp,T0)';
        [c,a] = ode15s(@(t,y)closedPlantODE(t,y,Kp), [0 tau],[u, Xp0]);
        Xp2(i,:) = a(end,:);
        
        dphip(i) = (phiFun(u, a(end,4:end)) - base.phip(end))/dr(i,i);
        dg1p(i) = (g1Fun(u, a(end,4:end)) - base.g1p(end))/dr(i,i);
        dg2p(i) = (g2Fun(u, a(end,4:end)) - base.g2p(end))/dr(i,i);
    end
    
    dphip = dphip*dy;
    dg1p = dg1p*dy;
    dg2p = dg2p*dy;
else
    dOpt.du = base.u(end,:) - u0_opt;
    dOpt.dC = base.Xp(end,:) - X0_opt;
    dfun = NEgrad(u0_opt+dOpt.du/2,dOpt); %+dOpt.du/2
    
    if meth == 1.0 %run NE with perfect dudr
        dudr = truedudr(rp0,base.Xp(end,:),@(r, x)plantController2(r, x, Kp, T0),dr)';
        
    elseif meth == 1.1
        dudr = pinv(dy);
        
    elseif meth == 1.2
        dudr = bsxfun(@times,pinv(bsxfun(@times,dy,u0_opt))',u0_opt)';
        
    elseif meth == 1.3
        dudr = bsxfun(@times,pinv(bsxfun(@times,dy,u0_opt))',u0_opt)';
        
    elseif meth == 1.4
        dudr = [0, 0, 1/dy(1,3); u0_opt(1)/u0_opt(2), 1, 0]';
        
    elseif meth == 1.5 %run NE with FE dudr
        kMax = ceil(tFinal/(tau));
        
        dudr = zeros(3,2);
        ord = [2,1]; %fastest to slowest
        u0 = base.u(end,:);
        
        for i = ord
            r = rp0 + dr(i,:);
            u = plantController2(r,base.Xp(end,:),Kp,T0)';
            
            [t,Xp] = ode15s(@(t,y)closedPlantODE(t,y,Kp), [0 tau/2],[u, base.Xp(end,:)]);
            n = numel(t);
            base.t(end+1:end+n) = t+base.t(end);
            base.u(end+1:end+n,:) = Xp(:,1:3);
            base.Xp(end+1:end+n,:) = Xp(:,4:end);
            
            base.phip = phiFun(up0,base.Xp);
            base.g1p = g1Fun(up0,base.Xp);
            base.g2p = g2Fun(up0,base.Xp);
            
            
            dudr(:,i) = (base.u(end,:) - u0)/dr(i,i);
        end
        
    else
        error('meth needs cannot be %d', meth)
    end
    
    dphip = (dfun.dphidu' + dphi0_NE)*dudr*dy;
    dg1p = (dfun.dg1du' + dg10_NE)*dudr*dy;
    dg2p = (dfun.dg2du' + dg20_NE)*dudr*dy;
end

% Get modifiers
m0phi = K*(base.phip(end) - phi0_opt);
m0g1 = K*(base.g1p(end) - g10_opt);
m0g2 = K*(base.g2p(end) - g20_opt);

m1phi = K*(dphip - dphi0_opt*pinv(dy)*dy);
m1g1 = K*(dg1p - dg10_opt*pinv(dy)*dy);
m1g2 = K*(dg2p - dg20_opt*pinv(dy)*dy);

% Make modified model
phiMod = @(u)(phiCU(u') + m0phi + m1phi*(u-u0_opt)');
g1Mod = @(u)(g1CU(u') + m0g1 + m1g1*(u-u0_opt)');
g2Mod = @(u)(g2CU(u') + m0g2 + m1g2*(u-u0_opt)');

% set-up iterations unsolved = 1;
unsolved = 1;
k = 1;
rGuess = rp0;
xGuess = Xp0;
uGuess = u0_opt;

while unsolved
    % Run model i
    [ui_opt(k,:)] = fmincon(@(u)phiMod(u),[4, 12, 90],...
        [],[],[],[],[0,0,60],[24,60,150],...
        @(u)deal([g1Mod(u),g2Mod(u)],[]),optionu);
    Xi_opt(k,:) = openModel(ui_opt(k,:), xGuess);
    
    % Get phi and g
    [phii_opt(k), dphii_opt] = phiCU(ui_opt(k,:)');
    [g1i_opt(k), dg1i_opt] = g1CU(ui_opt(k,:)');
    [g2i_opt(k), dg2i_opt] = g2CU(ui_opt(k,:)');
    
    dy = finDiff(@(u)yFromUX(u,openModel(u, xGuess)),ui_opt(k,:), 0.00001)';
    
    % Run plant for tau
    
    
    rpi(k,:) = yFromUX(ui_opt(k,:),Xi_opt(k,:));
    newXp = base.Xp(end,:);
    up = plantController2(rpi(k,:),newXp, Kp, T0)';
    if steady
        [t,Xp] = ode15s(@(t,y)closedPlantODE(t,y,Kp), [0:0.1:0.9,logspace(0,log10(tau),100)],[up, newXp]);
        
        % take until steady
        XpS = Xp(end,4:end);
        XpSn = all(bsxfun(@times,abs(Xp(:,4:end) - XpS),1./XpS)<0.0001,2);
        n1 = find(XpSn, 1, 'first');
        n2 = find(t>300, 1, 'first');
        
        n = max(n1,n2);
    else
        [t,Xp] = ode15s(@(t,y)closedPlantODE(t,y,Kp), [0 tau],[up, newXp]);
        
        % take the whole run
        n = numel(t);
    end
    
    base.t(end+1:end+n) = t(1:n)+base.t(end);
    base.u(end+1:end+n,:) = Xp(1:n,1:3);
    base.Xp(end+1:end+n,:) = Xp(1:n,4:end);
    
    base.ti(end+1) = base.t(end);
    % Get phi and g
    base.phip(end+1:end+n) = phiFun(Xp(1:n,1:3),Xp(1:n,4:end));
    base.g1p(end+1:end+n) = g1Fun(Xp(1:n,1:3),Xp(1:n,4:end));
    base.g2p(end+1:end+n) = g2Fun(Xp(1:n,1:3),Xp(1:n,4:end));
    
    % Estimate plant gradient
    if meth == -1
        dphip = zeros(1,2);
        dg1p = zeros(1,2);
        dg2p = zeros(1,2);
        
        % Perfect grad. est.
        xSteady = closedPlant(rpi(k,:),base.Xp(end,:), @(r, x)plantController2(r, x, Kp, T0));
        uSteady = plantController2(rpi(k,:), xSteady, Kp, T0)';
        
        for i = 1:2
            r = rpi(k,:) + dr(i,:);
            
            xSteady2 = closedPlant(r,xSteady, @(r, x)plantController2(r, x, Kp, T0));
            uSteady2 = plantController2(r, xSteady2, Kp, T0)';
            
            dphip(i) = (phiFun(uSteady2,xSteady2) - phiFun(uSteady,xSteady))/dr(i,i);
            dg1p(i) = (g1Fun(uSteady2,xSteady2) - g1Fun(uSteady,xSteady))/dr(i,i);
            dg2p(i) = (g2Fun(uSteady2,xSteady2) - g2Fun(uSteady,xSteady))/dr(i,i);
        end
        
        dphip = dphip*dy;
        dg1p = dg1p*dy;
        dg2p = dg2p*dy;
    elseif meth == 0 %run MU
        dphip = zeros(1,2);
        dg1p = zeros(1,2);
        dg2p = zeros(1,2);
        
        for i = 1:2
            r = rpi(k,:) + dr(i,:);
            u = plantController2(r,Xp2(i,4:end),Kp,T0)';
            [c,a] = ode15s(@(t,y)closedPlantODE(t,y,Kp), [0 tau],[u, Xp2(i,4:end)]);
            Xp2(i,:) = a(end,:);
            
            dphip(i) = (phiFun(u, a(end,4:end)) - base.phip(end))/dr(i,i);
            dg1p(i) = (g1Fun(u, a(end,4:end)) - base.g1p(end))/dr(i,i);
            dg2p(i) = (g2Fun(u, a(end,4:end)) - base.g2p(end))/dr(i,i);
        end
        
        dphip = dphip*dy;
        dg1p = dg1p*dy;
        dg2p = dg2p*dy;
    else
        dOpt.du = base.u(end,:) - u0_opt;%base.u(end,:) - u0_opt;
        dOpt.dC = base.Xp(end,:) - X0_opt;
        dfun0 = NEgrad(u0_opt+dOpt.du/2,dOpt);%+dOpt.du/2+dOpt.du/2
        pn = 3;
        dfun2.dphidu = dfun0.dphidu*0;
        dfun2.dg1du = dfun0.dg1du*0;
        dfun2.dg2du = dfun0.dg2du*0;
        for p = 0:pn-1
            dfun = NEgrad(u0_opt+(p/(pn-1))*dOpt.du,dOpt);%+dOpt.du/2
            dfun2.dphidu = dfun2.dphidu + dfun.dphidu*(1/pn);
            dfun2.dg1du = dfun2.dg1du + dfun.dg1du*(1/pn);
            dfun2.dg2du = dfun2.dg2du + dfun.dg2du*(1/pn);
        end
        dfun = dfun0;
        
        if meth == 1.0 %run NE with perfect dudr
            dudr = truedudr(rpi(k,:),base.Xp(end,:),@(r, x)plantController2(r, x, Kp, T0),dr)';
            
        elseif meth == 1.1
            dudr = pinv(dy);
            
        elseif meth == 1.2
            dudr = bsxfun(@times,pinv(bsxfun(@times,dy,u0_opt))',u0_opt)';
            
        elseif meth == 1.3
            dudr = bsxfun(@times,pinv(bsxfun(@times,dy,ui_opt(k,:)))',ui_opt(k,:))';
            
        elseif meth == 1.4
            dudr = [0, 0, 1/dy(1,3); u0_opt(1)/u0_opt(2), 1, 0]';
            
        elseif meth == 1.5 %run NE with FE dudr
            dudr = zeros(3,2);
            ord = [2,1]; %fastest to slowest
            u0 = base.u(end,:);
            
            for i = ord
                r = rpi(k,:) + dr(i,:);
                u = plantController2(r,base.Xp(end,:),Kp,T0)';
                
                [t,Xp] = ode15s(@(t,y)closedPlantODE(t,y,Kp), [0 tau/2],[u, base.Xp(end,:)]);
                n = numel(t);
                base.t(end+1:end+n) = t+base.t(end);
                base.u(end+1:end+n,:) = Xp(:,1:3);
                base.Xp(end+1:end+n,:) = Xp(:,4:end);
                
                % Get phi and g
                base.phip(end+1:end+n) = phiFun(up,Xp(:,4:end));
                base.g1p(end+1:end+n) = g1Fun(up,Xp(:,4:end));
                base.g2p(end+1:end+n) = g2Fun(up,Xp(:,4:end));
                
                dudr(:,i) = (base.u(end,:) - u0)/dr(i,i);
            end
            
        else
            error('meth needs cannot be %d', meth)
        end
        
        dphip = (dfun.dphidu' + dphi0_NE)*dudr*dy;
        dg1p = (dfun.dg1du' + dg10_NE)*dudr*dy;
        dg2p = (dfun.dg2du' + dg20_NE)*dudr*dy;
    end
    
    % Get modifiers
    m0phi = (1-K)*m0phi + K*(base.phip(end) - phii_opt(k));
    m0g1 = (1-K)*m0g1 + K*(base.g1p(end) - g1i_opt(k));
    m0g2 = (1-K)*m0g2 + K*(base.g2p(end) - g2i_opt(k));
    
    m1phi = (1-K)*m1phi + K*(dphip - dphii_opt*pinv(dy)*dy);
    m1g1 = (1-K)*m1g1 + K*(dg1p - dg1i_opt*pinv(dy)*dy);
    m1g2 = (1-K)*m1g2 + K*(dg2p - dg2i_opt*pinv(dy)*dy);
    
    % Make modified model
    phiMod = @(u)(phiCU(u') + m0phi + m1phi*(u-ui_opt(k,:))');
    g1Mod = @(u)(g1CU(u') + m0g1 + m1g1*(u-ui_opt(k,:))');
    g2Mod = @(u)(g2CU(u') + m0g2 + m1g2*(u-ui_opt(k,:))');
    
    k = k + 1;
    if base.t(end) > tFinal
        unsolved = 0;
    end
end

% solution analysis
solysis(rpi,base,rp_opt,phip_opt, 'UU', tau, K)

%plots
try fig = allPlots(rpi, rp0, base, 'UU',tau, K, fig);
catch
    fig = allPlots(rpi, rp0, base, 'UU',tau, K);
end

