% Transient Closed-loop UR
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
    tau = 300;
end

K = 0.7;

meth = 1.4;

noise = 0;

if noise
    pow = 0.02;
    freq = 60;
    
    [s{1:6}] = RandStream.create('mrg32k3a','NumStreams',6,'seed','shuffle');
end

% True optimum
optionu = optimoptions('fmincon','Display','off');
xGuess = [0.08, 0.37, 0.1, 0.25, 0.1, 0.1];

% @[T0 = 120, Kp = -1000]
rp_opt   = [0.10605,6.05325];
up_opt   = plantController2(rp_opt,xGuess,Kp,T0)';
[b,a]    = ode15s(@(t,y)closedPlantODE(t,y,Kp), [0:1:10000],[up_opt, xGuess]);
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
if ~noise
    
    [t,Xp] = ode15s(@(t,y)closedPlantODE(t,y,Kp), [0 tau],[up0, Xp0]);
    
    % take the whole run
    base.t = t-t(end);
    base.u = Xp(:,1:3);
    base.Xp = Xp(:,4:end);
    
    base.ti = base.t(end);
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
    
    base.phipn = phiFun(up0,base.Xpn);
    base.g1pn = g1Fun(up0,base.Xpn);
    base.g2pn = g2Fun(up0,base.Xpn);
end

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
    dOpt.dC = base.Xp(end,:) - X0_opt;%[base.Xp(end,1), base.u(end,2)];
    dfun = NEgrad(u0_opt+dOpt.du/2,dOpt);
    
    if meth == 1.0 %run NE with perfect dudr
        dudr = truedudr(rp0,base.Xp(end,:),@(r, x)plantController2(r, x, Kp, T0),dr)';
        
    elseif meth == 1.1
        dudr = pinv(dy);
        
    elseif meth == 1.2
        dy_norm = bsxfun(@times,bsxfun(@rdivide,dy',rp0)',u0_opt);
        dudr = bsxfun(@times,bsxfun(@rdivide,pinv(dy_norm),rp0)',u0_opt)';
        
    elseif meth == 1.3
        dy_norm = bsxfun(@times,bsxfun(@rdivide,dy',rp0)',u0_opt);
        dudr = bsxfun(@times,bsxfun(@rdivide,pinv(dy_norm),rp0)',u0_opt)';
        
    elseif meth == 1.4
        dudr = [0, 0, 1/dy(1,3); u0_opt(1)/u0_opt(2), 1, 0]';%-(1/dy(1,3)*dy(1,2)+1/dy(1,3)*dy(1,1)*(u0_opt(1)/u0_opt(2)))]';%u0_opt(1)/u0_opt(2)
        
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
    
    up = plantController2(rpi(k,:),newXp, Kp, T0)';
    
    if ~noise
        
        
        
        
        
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
        
                % Get phi and g
        base.phipn(end+1:end+n) = phiFun(up,Xp(:,4:9));
        base.g1pn(end+1:end+n) = g1Fun(up,Xp(:,4:9));
        base.g2pn(end+1:end+n) = g2Fun(up,Xp(:,4:9));
    end
    

    
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
            [c,a] = ode15s(@(t,y)closedPlantODE(t,y,Kp), [0 base.ti(end)-base.ti(end-1)],[u, Xp2(i,4:end)]);
            Xp2(i,:) = a(end,:);
            
            dphip(i) = (phiFun(u, a(end,4:end)) - base.phip(end))/dr(i,i);
            dg1p(i) = (g1Fun(u, a(end,4:end)) - base.g1p(end))/dr(i,i);
            dg2p(i) = (g2Fun(u, a(end,4:end)) - base.g2p(end))/dr(i,i);
        end
    else
        dOpt.du = base.u(end,:) - u0_opt(end,:);
        dOpt.dC = base.Xp(end,:) - X0_opt;
        dfun0 = NEgrad(u0_opt+dOpt.du/2,dOpt);%+dOpt.du/2
        pn = 21;
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
            dy_norm = bsxfun(@times,bsxfun(@rdivide,dy',rp0)',u0_opt);
            dudr = bsxfun(@times,bsxfun(@rdivide,pinv(dy_norm),rp0)',u0_opt)';
            
        elseif meth == 1.3
            dy_norm = bsxfun(@times,bsxfun(@rdivide,dy',rpi(k,:))',ui_opt(k,:));
            dudr = bsxfun(@times,bsxfun(@rdivide,pinv(dy_norm),rpi(k,:))',ui_opt(k,:))';
            
        elseif meth == 1.4
            dudr = [0, 0, 1/dy(1,3); ui_opt(end,1)/ui_opt(end,2), 1, 0]';%-(1/dy(1,3)*dy(1,2)+1/dy(1,3)*dy(1,1)*(ui_opt(end,1)/ui_opt(end,2)))]';
            
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
                    ui_opt(end,:)
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
    if base.t(end) > tFinal
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

%contours
X_Ac = linspace(min([rpi(:,1);rp0(1)])*0.9,max([rpi(:,1);rp0(1)])*1.1,21);
F_Bc = linspace(min([rpi(:,2);rp0(2)])*0.9,max([rpi(:,2);rp0(2)])*1.1,21);
fig = contPlot(X_Ac, F_Bc, fig);

