% Transient Closed-Loop
%addpath ConvexModel CSTR OtherFunctions PlotFunctions
clearvars -except fig
%close all

% variables
Kp = -1000;
T0 = 120;

tFinal = 10000;
steady = 0;
if steady
    tau = 2000;
else
    tau = 30;
end

K = 0.4;

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
%[r0_opt] = fmincon(@(x)phiCR(x'),[0.09, 12],[],[],[],[],...
%    [0,0],[1,50],@(x)deal([g1CR(x'),g2CR(x')],[]),optionu);
[u0_opt] = fmincon(@(u)phiFun(u,openModel(u, xGuess)),[4, 12, 90],...
    [],[],[],[],[0,0,60],[24,60,150],...
    @(u)deal([g1Fun(u,openModel(u, xGuess)),g2Fun(u,openModel(u, xGuess))],[]),optionu);
X0_opt = openModel(u0_opt, xGuess);
rp0 = yFromUX(u0_opt,X0_opt);

dphi0_NE = finDiff(@(u)phiFun(u,openModel(u, xGuess)), u0_opt, 0.00001)';
dg10_NE = finDiff(@(u)g1Fun(u,openModel(u, xGuess)), u0_opt, 0.00001)';
dg20_NE = finDiff(@(u)g2Fun(u,openModel(u, xGuess)), u0_opt, 0.00001)';

dy = finDiff(@(u)yFromUX(u,openModel(u, xGuess)),u0_opt, 0.00001)';

% Get phi and g
[phi0_opt, dphi0_opt] =  phiCR(rp0');
[g10_opt, dg10_opt] = g1CR(rp0');
[g20_opt, dg20_opt] = g2CR(rp0');

% Run plant @u0_opt
u0 = plantController2(rp0,xGuess,Kp,T0)';
[~,a] = ode15s(@(t,y)closedPlantODE(t,y,Kp), [0 10000],[u0, xGuess]);
u0_opt = a(end,1:3);
Xp0 = a(end,4:end);

% Run plant for tau
[t,Xp] = ode15s(@(t,y)closedPlantODE(t,y,Kp), [0 tau],[u0_opt, Xp0]);
base.t = t-t(end);
base.u = Xp(:,1:3);
base.Xp = Xp(:,4:end);

base.ti = base.t(end);

% Get phi and g
base.phip = phiFun(u0_opt,base.Xp);
base.g1p = g1Fun(u0_opt,base.Xp);
base.g2p = g2Fun(u0_opt,base.Xp);

dr = diag([0.001, 0.01]);

%Gradient Estimation
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
else
    dOpt.du = base.u(end,:) - u0_opt;
    dOpt.dC = base.Xp(end,:) - X0_opt;
    dfun = NEgrad(u0_opt+dOpt.du/2,dOpt); %
    
    if meth == 1.0 %run NE with perfect dudr
        dudr = truedudr(rp0,base.Xp(end,:),Kp,T0,dr)';
        
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
            
            base.phip = phiFun(u0_opt,base.Xp);
            base.g1p = g1Fun(u0_opt,base.Xp);
            base.g2p = g2Fun(u0_opt,base.Xp);
            
            
            dudr(:,i) = (base.u(end,:) - u0)/dr(i,i);
        end
        
    else
        error('meth needs cannot be %d', meth)
    end
    
    dphip = (dfun.dphidu' + dphi0_NE)*dudr;
    dg1p = (dfun.dg1du' + dg10_NE)*dudr;
    dg2p = (dfun.dg2du' + dg20_NE)*dudr;
end

% Get modifiers
m0phi = K*(base.phip(end) - phi0_opt);
m0g1 = K*(base.g1p(end) - g10_opt);
m0g2 = K*(base.g2p(end) - g20_opt);

m1phi = K*(dphip - dphi0_opt);
m1g1 = K*(dg1p - dg10_opt);
m1g2 = K*(dg2p - dg20_opt);

% Make modified model
phiMod = @(r)(phiCR(r) + m0phi + m1phi*(r-rp0'));
g1Mod = @(r)(g1CR(r) + m0g1 + m1g1*(r-rp0'));
g2Mod = @(r)(g2CR(r) + m0g2 + m1g2*(r-rp0'));

% set-up iterations
unsolved = 1;
k = 1;
rGuess = rp0;
xGuess = Xp(end,4:end);

while unsolved
    % Run model i
    [rpi(k,:)] = fmincon(@(x)phiMod(x'),rGuess,[],[],[],[],...
        [0,0],[1,50],@(x)deal([g1Mod(x'),g2Mod(x')],[]),optionu);
    
    
    % Get phi and g
    [phii_opt(k), dphii_opt] =  phiCR(rpi(k,:)');
    [g1i_opt(k), dg1i_opt] = g1CR(rpi(k,:)');
    [g2i_opt(k), dg2i_opt] = g2CR(rpi(k,:)');
    
        
    newXp = base.Xp(end,:);
    ui_opt(k,:) = plantController2(rpi(k,:),newXp, Kp, T0)';
    
    if steady
        [t,Xp] = ode15s(@(t,y)closedPlantODE(t,y,Kp), [0:0.1:0.9,logspace(0,log10(tau),100)],[ui_opt(k,:), newXp]);
        
        % take until steady
        XpS = Xp(end,4:end);
        XpSn = all(bsxfun(@times,abs(Xp(:,4:end) - XpS),1./XpS)<0.0001,2);
        n1 = find(XpSn, 1, 'first');
        n2 = find(t>300, 1, 'first');
        
        n = max(n1,n2);
    else
        [t,Xp] = ode15s(@(t,y)closedPlantODE(t,y,Kp), [0 tau],[ui_opt(k,:), newXp]);
        
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
        
    elseif meth == 0 %run MU
        for i = 1:2
            r = rpi(k,:) + dr(i,:);
            u = plantController2(r,Xp2(i,4:end),Kp,T0)';
            [c,a] = ode15s(@(t,y)closedPlantODE(t,y,Kp), [0 tau],[u, Xp2(i,4:end)]);
            Xp2(i,:) = a(end,:);
            
            dphip(i) = (phiFun(u, a(end,4:end)) - base.phip(end))/dr(i,i);
            dg1p(i) = (g1Fun(u, a(end,4:end)) - base.g1p(end))/dr(i,i);
            dg2p(i) = (g2Fun(u, a(end,4:end)) - base.g2p(end))/dr(i,i);
        end
    else
        dOpt.du = base.u(end,:) - u0_opt;
        dOpt.dC = base.Xp(end,:) - X0_opt;
        dfun = NEgrad(base.u(end,:) - dOpt.du/2,dOpt);
        
        if meth == 1.0 %run NE with perfect dudr
            dudr = truedudr(rpi(k,:),base.Xp(end,:),Kp,T0,dr)';
            
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
                base.phip(end+1:end+n) = phiFun(base.u(end,:),Xp(:,4:end));
                base.g1p(end+1:end+n) = g1Fun(ui_opt(k,:),Xp(:,4:end));
                base.g2p(end+1:end+n) = g2Fun(ui_opt(k,:),Xp(:,4:end));
                
                dudr(:,i) = (base.u(end,:) - u0)/dr(i,i);
            end
            
        else
            error('meth needs cannot be %d', meth)
        end
        
        dphip = (dfun.dphidu' + dphi0_NE)*dudr;
        dg1p = (dfun.dg1du' + dg10_NE)*dudr;
        dg2p = (dfun.dg2du' + dg20_NE)*dudr;
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
    if base.t(end) > tFinal
        unsolved = 0;
    end
end

% solution analysis
solysis(rpi,base,rp_opt,phip_opt,'RR',tau,K)

%plots
try fig = allPlots(rpi, rp0, base, 'RR',tau, K, fig);
catch
    fig = allPlots(rpi, rp0, base, 'RR',tau, K);
end





