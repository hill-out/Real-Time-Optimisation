% ------------------------------| T_closed |------------------------------
% A combined method for all closed algorithms using transient measurements
clc


%% changeable variables
if ~exist('holdVar')
    holdVar = 0;
end

if holdVar == 0 % Change variables
    clearvars -except fig
    
    % Set up the plant controller
    Kp = -1000;
    T0 = 120;
    controller = @(r, x)plantController2(r, x, Kp, T0);
    
    plantODE = @(t,y)closedPlantODE(t,y,Kp);
    
    % Set up the MA iteration variables
    tau = 500; % Iteration time-step
    tFinal = 10000; % Min end time
    kMax = ceil(tFinal/tau); % Number of iterations
    
    K = 0.2; % Gain of MA iterations
    
    % Set up method for gradient estimation
    % -------| Method Key |----------------
    % 0     Direct MU gradient estimation
    % 1     NE gradient estimation over plant
    % 1.0   Perfect dudr
    % 1.1   dudr = pinv(dydu) [fail]
    % 1.11  dudr = *scaled by nominal [u]* pinv(dydu) [fail]
    % 1.12  dudr = *scaled by nominal [u + y]* pinv(dydu) [fail]
    % 1.13  dudr = *scaled by iteration [u]* pinv(dydu) [fail]
    % 1.14  dudr = *scaled by iteration [u + y]* pinv(dydu) [fail]
    % 1.2   Limited controller knowledge estimation of dudr
    % 1.3   Transient FE estimation of dudr
    method = 1;
    
    % Set up noise
    noise = 0; %noise or not
    
    if noise
        pow = 0.01; % Relative amplitude
        freq = 60;  % Relative frequency
        
        % Set up noise
        [s{1:6}] = RandStream.create('mrg32k3a','NumStreams',6,'seed','shuffle');
        for no = 1:6
            Xpn(:,no) = [0;randn(s{no},2*freq,1)]*pow;
        end
    end
    
    % Set up closed loop algorithm
    % -------| Algo Key |---------
    % 0         UR method
    % 1         UU method
    % 2         RR method
    algo = 'UR';
    
    if ~any(strcmpi(algo,{'UR','UU','RR'})) %check if algo is valid
        warning('Invalid value for algo: ''%s''\nUsing ''UR'' instead\n',algo);
        algo = 'UR';
    end
    
else % Hold changeable variables
    clearvar -except fig controller tau tFinal kMax K method noise pow freq s algo
        
    % Reset holdVar
    holdVar = 0;
    
end

% %% Set up log file
% % get current files
% currLog = struct2cell(dir('Logs/*.txt'));
% currName = currLog(1,:);
% 
% % get last files
% last = currName{end};
% lastN = str2double(last(regexp(last,'\d')));
% 
% % set up new file name
% newN = lastN + 1;
% new = sprintf('Logs/log%03d.txt',newN);
% 
% % create file
% fid = fopen(new,'w');
% fprintf(fid,'-----| Log file %03d |-----| T_closed |-----| Created: %s |----',newN,datestr(now,0));
% fid = fclose(fid);

%% True optimum
optionu = optimoptions('fmincon','Display','off');
xGuess = [0.08, 0.37, 0.1, 0.25, 0.1, 0.1];

% Get setpoint
rp_opt = fmincon(@(r)phiFun(controller(r,closedPlant(r,xGuess,controller))',closedPlant(r,xGuess,controller)),[0.09, 6],[],[],[],[],[0,0],[1,60],...
    @(r)deal([g1Fun(controller(r,closedPlant(r,xGuess,controller)),closedPlant(r,xGuess,controller)),...
    g2Fun(controller(r,closedPlant(r,xGuess,controller)),closedPlant(r,xGuess,controller))],[]),optionu);

% Get X and u
Xp_opt   = closedPlant(rp_opt,xGuess,controller);
up_opt   = controller(rp_opt,Xp_opt)';

% Get phi and g
phip_opt = phiFun(up_opt,Xp_opt);
g1p_opt  = g1Fun(up_opt,Xp_opt);
g2p_opt  = g2Fun(up_opt,Xp_opt);

%% Set up for first iteration (get nominal values)
% Find nominal operating point
[u0_opt] = fmincon(@(u)phiFun(u,openModel(u, xGuess)),[4, 12, 90],...
    [],[],[],[],[0,0,60],[24,60,150],...
    @(u)deal([g1Fun(u,openModel(u, xGuess)),g2Fun(u,openModel(u, xGuess))],[]),optionu);
X0_opt = openModel(u0_opt, xGuess);

% Run plant to steady
rp0 = yFromUX(u0_opt,X0_opt);

Xp0   = closedPlant(rp0,xGuess,controller);
up0   = controller(rp0,Xp0)';

% Set up plant variable
plant.t   = -tau;
plant.u   = up0;
plant.Xp  = Xp0;
plant.Xpn = Xp0;

% Set up phi and g
switch algo
    case 'UR'
        phiBase = @(u)(phiFun(u,openModel(u, xGuess)));
        g1Base = @(u)(g1Fun(u,openModel(u, xGuess)));
        g2Base = @(u)(g2Fun(u,openModel(u, xGuess)));
        
        guess = up_opt;
        guessMin = [0,0,60];
        guessMax = [24,60,150];
        
        m1phi = [0,0];
        m1g1 = [0,0];
        m1g2 = [0,0];
        
    case 'UU'
        phiBase = @(u)(phiCU(u'));
        g1Base = @(u)(g1CU(u'));
        g2Base = @(u)(g2CU(u'));
        
        guess = up_opt;
        guessMin = [0,0,60];
        guessMax = [24,60,150];
        
        m1phi = [0,0,0];
        m1g1 = [0,0,0];
        m1g2 = [0,0,0];
        
    case 'RR'
        phiBase = @(r)(phiCR(r'));
        g1Base = @(r)(g1CR(r'));
        g2Base = @(r)(g2CR(r'));
        
        guess = rp_opt;
        guessMin = [0,0];
        guessMax = [1,50];
        
        m1phi = [0,0];
        m1g1 = [0,0];
        m1g2 = [0,0];
        
end

% gradient of mass fractions (y) wrt inputs (u)
switch algo
    case {'UR', 'UU'}
        % opt is u
        dphidu0 = finDiff(@(u)phiBase(u), u0_opt, 0.0001)';
        dg1du0 = finDiff(@(u)g1Base(u), u0_opt, 0.0001)';
        dg2du0 = finDiff(@(u)g2Base(u), u0_opt, 0.0001)';
        
        dydu = finDiff(@(u)yFromUX(u,openModel(u, xGuess)),u0_opt, 0.0001)';
        
        dphidr0 = dphidu0*pinv(dydu);
        dg1dr0 = dg1du0*pinv(dydu);
        dg2dr0 = dg2du0*pinv(dydu);
        
    case {'RR'}
        % opt is r
        dphidr0 = finDiff(@(u)phiBase(u), rp0, 0.0001);
        dg1dr0 = finDiff(@(u)g1Base(u), rp0, 0.0001);
        dg2dr0 = finDiff(@(u)g2Base(u), rp0, 0.0001);
        
        dphidu0 = dphidr0*dydu;
        dg1du0 = dg1dr0*dydu;
        dg2du0 = dg2dr0*dydu;
        
end

% Make modified model
phiMod = @(u)(phiBase(u));
g1Mod = @(u)(g1Base(u));
g2Mod = @(u)(g2Base(u));

% Initialise modifiers
m0phi = 0;
m0g1 = 0;
m0g2 = 0;

unsolved = 1;
k = 0;
dr = diag([0.001, 0.01]);


%% Run iterations
while unsolved
    %% Iteration Update
    k = k + 1;
    
    %% Run MA
    opt = fmincon(@(u)phiMod(u),guess,[],[],[],[],guessMin,guessMax,...
        @(u)deal([g1Mod(u),g2Mod(u)],[]));
    
    %% Run Model
    phi(k) = phiBase(opt);
    g1(k) = g1Base(opt);
    g2(k) = g2Base(opt);
    
    % get r (+ u and X if applicable)
    switch algo
        case {'UR', 'UU'}
            % opt is u
            uOpt(k,:) = opt;
            yOpt(k,:) = openModel(opt, xGuess);
            
            rp(k,:) = yFromUX(opt,yOpt(k,:));
            
        case {'RR'}
            % opt is r
            rp(k,:) = opt;
            
    end
    
    %% Get Model Gradients
    % gradient of mass fractions (y) wrt inputs (u)
    
    
    switch algo
        case {'UR', 'UU'}
            % opt is u
            dphidu = finDiff(@(u)phiBase(u), opt, 0.0001)';
            dg1du = finDiff(@(u)g1Base(u), opt, 0.0001)';
            dg2du = finDiff(@(u)g2Base(u), opt, 0.0001)';
            
            dydu = finDiff(@(u)yFromUX(u,openModel(u, xGuess)),opt, 0.0001)';
            
            dphidr = dphidu*pinv(dydu);
            dg1dr = dg1du*pinv(dydu);
            dg2dr = dg2du*pinv(dydu);
            
        case {'RR'}
            % opt is r
            dphidr = finDiff(@(u)phiBase(u), opt, 0.0001);
            dg1dr = finDiff(@(u)g1Base(u), opt, 0.0001);
            dg2dr = finDiff(@(u)g2Base(u), opt, 0.0001);
            
            dphidu = dphidr*dydu;
            dg1du = dg1dr*dydu;
            dg2du = dg2dr*dydu;
            
    end
    
    %% Run process
    up = controller(rp(k,:), plant.Xp(end,:))';
    
    if ~noise % noiseless
        [t,out] = ode15s(@(t,y)plantODE(t,y), [0 tau],[up, plant.Xp(end,:)]);
        
        % Get time variables
        plant.t  = [plant.t; plant.t(end) + t];
        plant.u  = [plant.u; out(:,1:3)];
        plant.Xp = [plant.Xp; out(:,4:end)];
        plant.Xpn = [plant.Xpn; out(:,4:end)];
        
    else
        % Set up new noise string
        newXpn = plant.Xpn(end,:);
        Xpn0 = Xpn;
        clear Xpn
        for no = 1:6
            Xpn(:,no) = [Xpn0(end-2*freq:end,no);randn(s{no},freq,1)*pow];
            
        end
        
        % Fit line to noise
        a = spline(linspace(-tau,2*tau,size(Xpn,1)),Xpn',linspace(-tau,2*tau,size(Xpn,1)*100));
        
        % Get differential
        b = diff(a');
        
        % Run process
        [t,out] = ode15s(@(t,y)closedPlantNoiseODE(t,y,Kp,...
            @(t)spline(linspace(-tau,2*tau,size(b,1)),b',t)),[0 tau],[up, newXpn, plant.Xp(end,:)-plant.Xpn(end,:)]);
        
        % Get time variables
        plant.t   = [plant.t; plant.t(end) + t];
        plant.u   = [plant.u; out(:,1:3)];
        plant.Xp  = [plant.Xp; out(:,4:9) + out(:,10:end)];
        plant.Xpn = [plant.Xpn; out(:,4:9)]; % measurement noiseless
        
    end
    
    % Get phi and g
    plant.phip = phiFun(plant.u,plant.Xp);
    plant.g1p  = g1Fun(plant.u,plant.Xp);
    plant.g2p  = g2Fun(plant.u,plant.Xp);
    
    % Get measurement noiseless phi and g    
    plant.phipn = phiFun(plant.u,plant.Xpn);
    plant.g1pn  = g1Fun(plant.u,plant.Xpn);
    plant.g2pn  = g2Fun(plant.u,plant.Xpn);
    
    %% Estimate plant gradient
    if method == 0 %run MU
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
        dOpt.du = plant.u(end,:) - uOpt(1,:);
        dOpt.dC = plant.Xp(end,:) - yOpt(1,:);
        dfun = NEgrad(uOpt(1,:)+dOpt.du/2,dOpt);
%         pn = 21;
%         dfun2.dphidu = dfun0.dphidu*0;
%         dfun2.dg1du = dfun0.dg1du*0;
%         dfun2.dg2du = dfun0.dg2du*0;
%         for p = 0:pn-1
%             dfun = NEgrad(u0_opt+(p/(pn-1))*dOpt.du,dOpt);%+dOpt.du/2
%             dfun2.dphidu = dfun2.dphidu + dfun.dphidu*(1/pn);
%             dfun2.dg1du = dfun2.dg1du + dfun.dg1du*(1/pn);
%             dfun2.dg2du = dfun2.dg2du + dfun.dg2du*(1/pn);
%         end
%         dfun = dfun2;
        
        if method == 1.0 %run NE with perfect dudr
            dudr = truedudr(rp(k,:),plant.Xp(end,:),controller,dr)';
            
        elseif method == 1.1
            dudr = pinv(dy);
            
        elseif method == 1.2
            dy_norm = bsxfun(@times,bsxfun(@rdivide,dy',rp0)',u0_opt);
            dudr = bsxfun(@times,bsxfun(@rdivide,pinv(dy_norm),rp0)',u0_opt)';
            
        elseif method == 1.3
            dy_norm = bsxfun(@times,bsxfun(@rdivide,dy',rpi(k,:))',ui_opt(k,:));
            dudr = bsxfun(@times,bsxfun(@rdivide,pinv(dy_norm),rpi(k,:))',ui_opt(k,:))';
            
        elseif method == 1.4
            dudr = [0, 0, 1/dy(1,3); u0_opt(1)/u0_opt(2), 1, 0]';
            
        elseif method == 1.5 %run NE with FE dudr
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
            error('method needs cannot be %d', meth)
        end
        
        dphipdr = (dfun.dphidu'*dudr + dphidu0*dudr);
        dg1pdr = (dfun.dg1du'*dudr + dg1du0*dudr);
        dg2pdr = (dfun.dg2du'*dudr + dg2du0*dudr);
    end
    
    %% Calculate the Modifiers
    % zeroth order
    m0phi = (1-K)*m0phi + K*(plant.phip(end) - phi(k));
    m0g1 = (1-K)*m0g1 + K*(plant.g1p(end) - g1(k));
    m0g2 = (1-K)*m0g2 + K*(plant.g2p(end) - g2(k));
    
    %first order
    switch algo
        case {'UR', 'RR'}
            % modification is in r
            m1phi = (1-K)*m1phi + K*(dphipdr - dphidr);
            m1g1 = (1-K)*m1g1 + K*(dg1pdr - dg1dr);
            m1g2 = (1-K)*m1g2 + K*(dg2pdr - dg2dr);
            
        case {'UU'}
            % modification is in u            
            m1phi = (1-K)*m1phi + K*(dphipdr - dphidr)*dydu;
            m1g1 = (1-K)*m1g1 + K*(dg1pdr - dg1dr)*dydu;
            m1g2 = (1-K)*m1g2 + K*(dg2pdr - dg2dr)*dydu;
            
    end
    
    %% Calculate new modified model
    switch algo
        case {'UR', 'RR'}
            % modification is in r
            phiMod = @(u)(phiBase(u) + m0phi + m1phi*(yFromUX(u,openModel(u, xGuess))-rp(k,:))');
            g1Mod = @(u)(g1Base(u) + m0g1 + m1g1*(yFromUX(u,openModel(u, xGuess))-rp(k,:))');
            g2Mod = @(u)(g2Base(u) + m0g2 + m1g2*(yFromUX(u,openModel(u, xGuess))-rp(k,:))');
            
        case {'UU'}
            % modification is in u
            phiMod = @(u)(phiBase(u) + m0phi + m1phi*(u-opt)');
            g1Mod = @(u)(g1Base(u) + m0g1 + m1g1*(u-opt)');
            g2Mod = @(u)(g2Base(u) + m0g2 + m1g2*(u-opt)');
            
    end
    
    if k > kMax
        unsolved = 0;
    end
    
end










