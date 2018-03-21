% Transient Closed-loop UR PI
%addpath ConvexModel CSTR OtherFunctions PlotFunctions
clearvars -except fig
%close all

% variables
Kp = -2000;
Ki = -2.0;
T0 = 120;

tau = 300;
K = 0.8;

tFinal = 6000;
kMax = ceil(tFinal/tau);

meth = 1;

noise = 0;
if noise
    pow = 0.01;
    freq = 100;
    
    [s{1:6}] = RandStream.create('mrg32k3a','NumStreams',6,'seed','shuffle');
end

% globals
global error
error = 0;

% True optimum
optionu = optimoptions('fmincon','Display','off');
xGuess = [0.09, 0.36, 0.1, 0.25, 0.1, 0.1];

% @[T0 = 120, Kp = -1000]
% rp_opt   = [0.10605,6.05325];
% up_opt   = plantController2(rp_opt,xGuess,Kp,T0)';
% [b,a]    = ode15s(@(t,y)closedPlantODE(t,y,Kp), [0:1:10000],[up_opt, xGuess]);
% up_opt   = a(end,1:3);
% Xp_opt   = a(end,4:end);
% phip_opt = phiFun(up_opt,Xp_opt);
% g1p_opt  = g1Fun(up_opt,Xp_opt);
% g2p_opt  = g2Fun(up_opt,Xp_opt);

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
up = plantControllerPI(rp0,X0_opt,Kp,Ki,T0)';
options = @(r)odeset('OutputFcn',@(t,y,flag) myOutputFcn(t,y,flag,  Kp, Ki, T0, r));

[a0,a] = ode45(@(t,y)closedPlantODE_PI(t,y,Kp,Ki,T0,rp0), [0 10000],[X0_opt],options(rp0));
clear closedPlantODE_PI
up0 = uODE(end,:);
Xp0 = a(end,:);

% Run plant for tau
if ~noise
    
    [t,Xp] = ode15s(@(t,y)closedPlantODE_PI(t,y,Kp,Ki,T0,rp0), [0 tau],[Xp0],options(rp0));
    clear closedPlantODE_PI

    base.t = t-t(end);
    base.u = uODE;
    base.Xp = Xp;
else
    
    for no = 1:6
        state(:,no) = s{no}.State;
        Xpn(:,no) = [0;randn(s{no},2*freq,1)]*pow;
    end
    
    a = spline(linspace(0,2*tau,size(Xpn,1)),Xpn',linspace(0,2*tau,size(Xpn,1)*100));
    b = diff(a');
    
    [t,Xp] = ode15s(@(t,y)closedPlantNoiseODE(t,y,Kp,...
        @(t)spline(linspace(0,2*tau,size(b,1)),b',t)),[0 tau],[up0, Xp0, Xpn(1,:)*0]);
    
    base.t = t-t(end);
    base.u = Xp(:,1:3);
    base.Xpn = Xp(:,4:9);
    base.Xp = Xp(:,4:9)+Xp(:,10:end);
end

base.phip = phiFun(up0,base.Xp);
base.g1p = g1Fun(up0,base.Xp);
base.g2p = g2Fun(up0,base.Xp);

dr = diag([0.001, 0.01]);

if meth == 0 %run MU
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
    dfun = NEgrad(u0_opt+dOpt.du/2,dOpt);
    
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
        kMax = ceil(tFinal/(2*tau));
        
        dudr = zeros(3,2);
        ord = [2,1]; %fastest to slowest
        u0 = base.u(end,:);
        if ~noise
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
            for i = ord
                r = rp0 + dr(i,:);
                u = plantController2(r,base.Xp(end,:),Kp,T0)';
                
                Xpn0 = Xpn;
                clear Xpn
                for no = 1:6
                    state(:,no) = s{no}.State;
                    Xpn(:,no) = [Xpn0(end-2*freq:end,no);randn(s{no},freq,1)*pow];
                end
                
                a = spline(linspace(-tau,2*tau,size(Xpn,1)),Xpn',linspace(-tau,2*tau,size(Xpn,1)*100));
                b = diff(a');
                
                [t,Xp] = ode15s(@(t,y)closedPlantNoiseODE(t,y,Kp,...
                    @(t)spline(linspace(-tau,2*tau,size(b,1)),b',t)),[0 tau],[u, base.Xpn(end,:), base.Xp(end,:)-base.Xpn(end,:)]);
                
                n = numel(t);
                base.t(end+1:end+n) = t+base.t(end);
                base.u(end+1:end+n,:) = Xp(:,1:3);
                base.Xpn(end+1:end+n,:) = Xp(:,4:9);
                base.Xp(end+1:end+n,:) = Xp(:,4:9)+Xp(:,10:end);
                
                base.phip(end+1:end+n) = phiFun(u,Xp(:,4:9)+Xp(:,10:end));
                base.g1p(end+1:end+n) = g1Fun(u,Xp(:,4:9)+Xp(:,10:end));
                base.g2p(end+1:end+n) = g2Fun(u,Xp(:,4:9)+Xp(:,10:end));
                
                dudr(:,i) = (base.u(end,:) - u0)/dr(i,i);
            end
        end
        
        
    else
        error('meth needs cannot be %d', meth)
    end
    
    dphip = (dfun.dphidu'*dudr + dphi0_opt*dudr);
    dg1p = (dfun.dg1du'*dudr + dg10_opt*dudr);
    dg2p = (dfun.dg2du'*dudr + dg20_opt*dudr);
end

% Get modifiers
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
    
    up = plantControllerPI(rpi(k,:),X0_opt,Kp,Ki,T0)';
    
    if ~noise
        
        [t,Xp] = ode45(@(t,y)closedPlantODE_PI(t,y,Kp,Ki,T0,rpi(k,:)), [0:1:tau],[newXp],options(rpi(k,:)));
        clear closedPlantODE_PI
        
        base.t = [base.t; base.t(end) + t];
        base.u = [base.u; uODE];
        base.Xp = [base.Xp; Xp];
                
        % Get phi and g
        base.phip = [base.phip; phiFun(up,Xp)];
        base.g1p = [base.g1p; g1Fun(up,Xp)];
        base.g2p = [base.g2p; g2Fun(up,Xp)];
    else
        newXpn = base.Xpn(end,:);
        %[s{1:6}] = RandStream.create('mrg32k3a','NumStreams',6,'seed',time/2);
        Xpn0 = Xpn;
        clear Xpn
        for no = 1:6
            state(:,no) = s{no}.State;
            Xpn(:,no) = [Xpn0(end-2*freq:end,no);randn(s{no},freq,1)*pow];
        end
        
        a = spline(linspace(-tau,2*tau,size(Xpn,1)),Xpn',linspace(-tau,2*tau,size(Xpn,1)*100));
        b = diff(a');
        
        [t,Xp] = ode15s(@(t,y)closedPlantNoiseODE(t,y,Kp,...
            @(t)spline(linspace(-tau,2*tau,size(b,1)),b',t)),[0 tau],[up, newXpn, base.Xp(end,:)-base.Xpn(end,:)]);
        
        n = numel(t);
        base.t(end+1:end+n) = t+base.t(end);
        base.u(end+1:end+n,:) = Xp(:,1:3);
        base.Xpn(end+1:end+n,:) = Xp(:,4:9);
        base.Xp(end+1:end+n,:) = Xp(:,4:9)+Xp(:,10:end);
        
        % Get phi and g
        base.phip(end+1:end+n) = phiFun(up,Xp(:,4:9)+Xp(:,10:end));
        base.g1p(end+1:end+n) = g1Fun(up,Xp(:,4:9)+Xp(:,10:end));
        base.g2p(end+1:end+n) = g2Fun(up,Xp(:,4:9)+Xp(:,10:end));
    end
    

    
    % Estimate plant gradient
    if meth == 0 %run MU
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
        dOpt.du = base.u(end,:) - ui_opt(k,:);
        dOpt.dC = base.Xp(end,:) - Xi_opt(k,:);
        dfun = NEgrad(base.u(end,:) - dOpt.du/2, dOpt);
        
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
            if ~noise
                for i = ord
                    r = rpi(k,:) + dr(i,:);
                    u = plantController2(r,base.Xp(end,:),Kp,T0)';
                    
                    [t,Xp] = ode15s(@(t,y)closedPlantODE(t,y,Kp), [0 tau/2],[u, base.Xp(end,:)]);
                    n = numel(t);
                    base.t(end+1:end+n) = t+base.t(end);
                    base.u(end+1:end+n,:) = Xp(:,1:3);
                    base.Xp(end+1:end+n,:) = Xp(:,4:end);
                    
                    % Get phi and g
                    base.phip(end+1:end+n) = phiFun(u,Xp(:,4:end));
                    base.g1p(end+1:end+n) = g1Fun(u,Xp(:,4:end));
                    base.g2p(end+1:end+n) = g2Fun(u,Xp(:,4:end));
                    
                    dudr(:,i) = (base.u(end,:) - u0)/dr(i,i);
                end
            else
                for i = ord
                    r = rpi(k,:) + dr(i,:);
                    u = plantController2(r,base.Xp(end,:),Kp,T0)';
                    
                    Xpn0 = Xpn;
                    clear Xpn
                    for no = 1:6
                        state(:,no) = s{no}.State;
                        Xpn(:,no) = [Xpn0(end-2*freq:end,no);randn(s{no},freq,1)*pow];
                    end
                    
                    a = spline(linspace(-tau,2*tau,size(Xpn,1)),Xpn',linspace(-tau,2*tau,size(Xpn,1)*100));
                    b = diff(a');
                    
                    [t,Xp] = ode15s(@(t,y)closedPlantNoiseODE(t,y,Kp,...
                        @(t)spline(linspace(-tau,2*tau,size(b,1)),b',t)),[0 tau],[u, base.Xpn(end,:), base.Xp(end,:)-base.Xpn(end,:)]);
                    
                    n = numel(t);
                    base.t(end+1:end+n) = t+base.t(end);
                    base.u(end+1:end+n,:) = Xp(:,1:3);
                    base.Xpn(end+1:end+n,:) = Xp(:,4:9);
                    base.Xp(end+1:end+n,:) = Xp(:,4:9)+Xp(:,10:end);
                    
                    base.phip(end+1:end+n) = phiFun(u,Xp(:,4:9)+Xp(:,10:end));
                    base.g1p(end+1:end+n) = g1Fun(u,Xp(:,4:9)+Xp(:,10:end));
                    base.g2p(end+1:end+n) = g2Fun(u,Xp(:,4:9)+Xp(:,10:end));
                    
                    dudr(:,i) = (base.u(end,:) - u0)/dr(i,i);
                end
            end
        else
            error('meth needs cannot be %d', meth)
        end
        
        dphip = (dfun.dphidu'*dudr + dphii_opt*dudr);
        dg1p = (dfun.dg1du'*dudr + dg1i_opt*dudr);
        dg2p = (dfun.dg2du'*dudr + dg2i_opt*dudr);
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

% solution analysis
%solysis(rpi,base,rp_opt,phip_opt,'UR',tau,K)

%plots
try fig = allPlots(rpi, rp0, base, 'UR',tau, K, fig);
catch
    fig = allPlots(rpi, rp0, base, 'UR',tau, K);
end


