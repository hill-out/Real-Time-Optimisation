function [model, plant] = runNE(model, plant, gain, tau)
% -------------------------------------------------------------------------
% This function does MU-RTO on the plant using a model
%
% model         - struct        - structure of all the model parameters
% plant         - struct        - structure of all the plant parameters
% gain          - double        - the gain bewteen each iteration
% tau           - double        - time step between each iteration
%
% model         - struct        - updated input
% plant         - struct        - updated input
%
% -------------------------------------------------------------------------

%% inital iteration (plant at steady state, modelOpt conditions)
u0 = model.uOpt;
[mC0, mCost, mCons] = CSTRste(u0, model);
[pC, pCost, pCons] = CSTRste(u0, plant);
delu = 0.01;

% initialise base unit and offset units
plant.base = struct('u', u0,...
    't', 0,...
    'c', pC,...
    'cost', pCost,...
    'cons', pCons);

% calculate zeroth order midifier (cost)
pCost = plant.base.cost(end);

epiCost = pCost - mCost;

% calculate zeroth order midifier (cons)
pCons = plant.base.cons(:,end);

epiCons = pCons - mCons;

% calculate first order midifier
mCostGrad0 = gradFD(@(x)(model.cost(x, CSTRste(x, model))), 1,  u0, delu);
mConsGrad0 = gradFD(@(x)(model.cons(x, CSTRste(x, model))), 2,  u0, delu);

[pCostGrad, pConsGrad] = gradNE;

lamCost = pCostGrad - mCostGrad0';
lamCons = pConsGrad - mConsGrad0;

%apply gain
% epiCost = (gain)*epiCost;
% epiCons = (gain)*epiCons;
% lamCost = (gain)*lamCost;
% lamCons = (gain)*lamCons;

% modify model
modified.cost = @(u)(model.cost(u, CSTRste(u, model)) + epiCost + (u-u0)*lamCost);
modified.cons = @(u)(model.cons(u, CSTRste(u, model)) + epiCons + ((u-u0)*lamCons)');

RTO.i = 0;
RTO.u = reshape(u0,1,[]);
RTO.epiCost = reshape(epiCost,1,[]);
RTO.epiCons = reshape(epiCons,1,[]);
RTO.lamCost = reshape(lamCost,1,[]);
RTO.lamCons = reshape(lamCons,1,[]);

%% run for all other iterations
unsolved = 1;

while unsolved
    
    u0 = CSTRopt(@(x,y)(modified.cost(x)), @(x,y)(modified.cons(x)), model, RTO.u(end,:));
    [~, mCost, mCons] = CSTRste(u0, model);
    
    % run plant
    plant.base.u (end+1,:) = u0;
    
    [pCostGrad, pConsGrad] = gradNE;
    pCost = plant.base.cost(end);
    pCons = plant.base.cons(:,end);
    
    % calculate zeroth order midifier
    epiCost = pCost - mCost;
    epiCons = pCons - mCons;
    
    % calculate first order midifier
    mCostGrad = gradFD(@(x)(model.cost(x, CSTRste(x, model))), 1,  u0, delu);
    mConsGrad = gradFD(@(x)(model.cons(x, CSTRste(x, model))), 2,  u0, delu);
    
    lamCost = pCostGrad - mCostGrad';
    lamCons = pConsGrad - mConsGrad;
    
    % apply gain
    epiCost = RTO.epiCost(end,:)*(1-gain) + (gain)*epiCost;
    epiCons = reshape(RTO.epiCons(end,:),2,[])*(1-gain) + (gain)*epiCons;
    lamCost = reshape(RTO.lamCost(end,:),2,[])*(1-gain) + (gain)*lamCost;
    lamCons = reshape(RTO.lamCons(end,:),[],2)*(1-gain) + (gain)*lamCons;
    
    % modify model
    modified.cost = @(u)(model.cost(u, CSTRste(u, model)) + epiCost + (u-u0)*lamCost);
    modified.cons = @(u)(model.cons(u, CSTRste(u, model)) + epiCons + ((u-u0)*lamCons)');
    
    % save data
    RTO.i(end+1) = RTO.i(end) + 1;
    RTO.u(end+1,:) = reshape(u0,1,[]);
    RTO.epiCost(end+1,:) = reshape(epiCost,1,[]);
    RTO.epiCons(end+1,:) = reshape(epiCons,1,[]);
    RTO.lamCost(end+1,:) = reshape(lamCost,1,[]);
    RTO.lamCons(end+1,:) = reshape(lamCons,1,[]);
    
    if RTO.i(end) == 140
        unsolved = 1;
    end
    %plotConsArea(modified.cons);
    %pause(2)
end

    function [costGrad, consGrad] = gradNE
        % ------------------------------------------------------------------
        % runs the NE method of gradient estimation of a plant
        %
        % costGrad      - double        - gradient of the cost
        % consgrad      - double        - gradient of the cons
        %
        % -----------------------------------------------------------------
        
        % run the plant
        [~, t, c, cost, cons] = ...
            CSTRdyn(plant.base.u(end,:), plant.base.c(:,end), plant, tau);
        
        nIter = numel(t);
        
        plant.base.t(end+1:end+nIter) = plant.base.t(end) + t;
        plant.base.c(:,end+1:end+nIter) = c;
        plant.base.cost(end+1:end+nIter) = cost;
        plant.base.cons(:,end+1:end+nIter) = cons;
        
        % first order gradients of steady state model (du)
        [C0, cost0, cons0] = CSTRste(u0,model);
        
        du = repmat(u0, 2, 1) + 0.001*eye(2); %modified values of u
        Cdu = zeros(2,1,4); %values of c at du
        costdu0 = zeros(2,1); %values of cost at du
        consdu0 = zeros(2,1,2); %values of cons at du
        
        for du_i = 1:2
            [Cdu(du_i,1,:), costdu0(du_i), consdu0(du_i,1,:)] = CSTRste(du(du_i,:),model);
        end
        dcostdu0 = (costdu0 - repmat(cost0,2,1))./0.001;
        dconsdu0 = (consdu0 - repmat(permute(cons0,[3,2,1]),2,1))./0.001;
        dHdu = (Cdu - repmat(reshape(C0,1,1,4),2,1))./0.001;
        
        % first order gradients of steady state model (dp)
        p = [model.k(1), model.k(2), model.c_in(1)];
        dp = repmat(p, 3, 1) + 0.001*eye(3);
        Cdp = zeros(1,3,4);
        
        for dp_i = 1:3
            newModel = model;
            newModel.k = dp(dp_i,1:2);
            newModel.c_in(1) = dp(dp_i,3);
            
            
            
            Cdp(1,dp_i,:) = CSTRste(u0,newModel);
            
        end
        
        dHdp = (Cdp - repmat(reshape(C0,1,1,4),1,3))/0.001;
        
        % second order gradients for cost (dudp)
        dcostdudp = zeros(2,3);
        
        for dp_i = 1:3
            
            costdu = zeros(2,1);
            newModel = model;
            newModel.k = dp(dp_i,1:2);
            newModel.c_in(1) = dp(dp_i,3);
            
            newModel.cost = @(u,c)(costFun(u,c,newModel));
            [~,cost0] = CSTRste(u0,newModel);
            
            for du_j = 1:2
                [~, costdu(du_j,1)] = CSTRste(du(du_j,:),newModel);
            end
            
            dcostdudp(:,dp_i) = (costdu - repmat(cost0,2,1))/0.001;
        end
        
        ddcostdudp = (dcostdudp - repmat(dcostdu0,1,3))./0.001;
        
        % second order gradients for cost (dudu)
        eye3 = repmat(eye(2),1,1,2);
        dudu = repmat(u0, 2, 1, 2) + 0.001*eye3 + 0.001*permute(eye3,[3,2,1]); %modified values of u
        costdudu = zeros(2,1); %values of cost at du
        dcostdudu = zeros(2);
        
        for du_i = 1:2
            
            [~,cost0] = CSTRste(du(du_i,:),model);
            
            for du_j = 1:2
                [~, costdudu(du_j)] = CSTRste(dudu(du_i,:,du_j),model);
            end
            
            dcostdudu(:,du_i) = (costdudu - repmat(cost0,2,1))/0.001;
        end
        
        ddcostdudu = (dcostdudu - repmat(dcostdu0,1,2))./0.001;
        
                
        % second order gradients for cons (dudp)
        dconsdudp = zeros(2,3,2);
        
        for dp_i = 1:3
            
            consdu = zeros(2,1,2);
            newModel = model;
            newModel.k = dp(dp_i,1:2);
            newModel.c_in(1) = dp(dp_i,3);
            
            newModel.cons = @(u,c)(consFun(u,c,newModel));
            
            [~,~,cons0] = CSTRste(u0,newModel);
            
            for du_j = 1:2
                [~, ~, consdu(du_j,1,:)] = CSTRste(du(du_j,:),newModel);
            end
            
            dconsdudp(:,dp_i,:) = (consdu - repmat(permute(cons0,[3,2,1]),2,1))/0.001;
        end
        
        ddconsdudp = (dconsdudp - repmat(dconsdu0,1,3))./0.001;
        
        % second order gradients for cost (dudu)
        consdudu = zeros(2,1,2); %values of cost at du
        dconsdudu = zeros(2,2,2);
        
        for du_i = 1:2
            
            [~,~,cons0] = CSTRste(du(du_i,:),model);
            
            for du_j = 1:2
                [~, ~, consdudu(du_j,1,:)] = CSTRste(dudu(du_i,:,du_j),model);
            end
            
            dconsdudu(:,du_i,:) = (consdudu - repmat(permute(cons0,[3,2,1]),2,1))/0.001;
        end
        
        ddconsdudu = (dconsdudu - repmat(dconsdu0,1,2))./0.001;
        
        invdHdp = pinv(permute(dHdp,[3,2,1]));
        
        dyp = plant.base.c(:,end) - mC0;
        duOpt = u0 - model.uOpt;
        
        costGrad = mCostGrad0' + ddcostdudp*invdHdp*dyp + ...
            (ddcostdudu - ddcostdudp*invdHdp*permute(dHdu,[3,1,2]))*duOpt';
        
        consGrad = zeros(2);
        
        consGrad(1,:) = mConsGrad0(1,:)' + ddconsdudp(:,:,1)*invdHdp*dyp + ...
            (ddconsdudu(:,:,1) - ddconsdudp(:,:,1)*invdHdp*permute(dHdu,[3,1,2]))*duOpt';
        consGrad(2,:) = mConsGrad0(2,:)' + ddconsdudp(:,:,2)*invdHdp*dyp + ...
            (ddconsdudu(:,:,2) - ddconsdudp(:,:,2)*invdHdp*permute(dHdu,[3,1,2]))*duOpt';

    end
end