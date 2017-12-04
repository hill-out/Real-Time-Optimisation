function [model, plant] = runMU(model, plant, gain, tau)
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
[mC, mCost, mCons] = CSTRste(u0, model);
[pC, pCost, pCons] = CSTRste(u0, plant);
du = 0.01;

% initialise base unit and offset units
plant.base = struct('u', u0,...
    't', 0,...
    'c', pC,...
    'cost', pCost,...
    'cons', pCons);

u0_1 = u0;
u0_1(1) = u0_1(1) + du;

u0_2 = u0;
u0_2(2) = u0_2(2) + du;

plant.offset = struct('u', {u0_1, u0_2},...
    't', 0,...
    'c', pC,...
    'cost', pCost,...
    'cons', pCons);

[pCostGrad, pConsGrad] = gradMU;

% calculate zeroth order midifier (cost)
pCost = plant.base.cost(end);

epiCost = pCost - mCost;

% calculate zeroth order midifier (cons)
pCons = plant.base.cons(end);

epiCons = pCons - mCons;

% calculate first order midifier
mCostGrad = gradFD(@(x)(model.cost(x, CSTRste(x, model))), 1,  u0, du);
mConsGrad = gradFD(@(x)(model.cons(x, CSTRste(x, model))), 2,  u0, du);

lamCost = (pCostGrad - mCostGrad)/model.costOpt;
lamCons = pConsGrad - mConsGrad;

% modify model
modified.cost = @(u,c)(model.cost(u,c) + epiCost + lamCost*u');
modified.cons = @(u,c)(model.cons(u,c) + epiCons + lamCons*u');

RTO.i = 0;
RTO.u = reshape(u0,1,[]);
RTO.epiCost = reshape(epiCost,1,[]);
RTO.epiCons = reshape(epiCons,1,[]);
RTO.lamCost = reshape(lamCost,1,[]);
RTO.lamCons = reshape(lamCons,1,[]);

%% run for all other iterations
unsolved = 1;

while unsolved
    
    u0 = CSTRopt(modified.cost, modified.cons, model, RTO.u(end,:));
    [mC, mCost, mCons] = CSTRste(u0, model);
    [pC, pCost, pCons] = CSTRste(u0, plant);
    
    % run plant
    u0_1 = u0;
    u0_1(1) = u0_1(1) + du;

    u0_2 = u0;
    u0_2(2) = u0_2(2) + du;
    
    plant.base.u (end+1,:) = u0;
    plant.offset(1).u(end+1,:) = u0_1;
    plant.offset(2).u(end+1,:) = u0_2;
    
    [pCostGrad, pConsGrad] = gradMU;
    
    % calculate zeroth order midifier
    epiCost = pCost - mCost;
    epiCons = pCons - mCons;
    
    % calculate first order midifier    
    mCostGrad = gradFD(@(x)(model.cost(x, CSTRste(x, model))), 1,  u0, du);
    mConsGrad = gradFD(@(x)(model.cons(x, CSTRste(x, model))), 2,  u0, du);
    
    lamCost = (pCostGrad - mCostGrad)/model.costOpt;
    lamCons = pConsGrad - mConsGrad;
    
    % apply gain
    epiCost = RTO.epiCost(end,:)*gain - (1-gain)*epiCost;
    epiCons = RTO.epiCost(end,:)*gain - (1-gain)*epiCons;
    lamCost = reshape(RTO.lamCost(end,:),[],2)*gain - (1-gain)*lamCost;
    lamCons = reshape(RTO.lamCons(end,:),[],2)*gain - (1-gain)*lamCons;
    
    % modify model
    modified.cost = @(u,c)(model.cost(u,c) + epiCost + lamCost*u');
    modified.cons = @(u,c)(model.cons(u,c) + epiCons + lamCons*u');
    
    RTO.i(end+1) = RTO.i(end) + 1;
    RTO.u(end+1,:) = reshape(u0,1,[]);
    RTO.epiCost(end+1,:) = reshape(epiCost,1,[]);
    RTO.epiCons(end+1,:) = reshape(epiCons,1,[]);
    RTO.lamCost(end+1,:) = reshape(lamCost,1,[]);
    RTO.lamCons(end+1,:) = reshape(lamCons,1,[]);
    
end


    function [costGrad, consGrad] = gradMU
        % -----------------------------------------------------------------
        % runs the MU method of gradient estimation of a plant
        %
        % costGrad      - double        - gradient of the cost
        % consgrad      - double        - gradient of the cons
        %
        % -----------------------------------------------------------------
        
        costGrad = zeros(1:2);
        consGrad = zeros(2);
        
        [~, t, c, cost, cons] = ...
            CSTRdyn(plant.base.u(end,:), plant.base.c(:,end), plant, tau);
        nIter = numel(t);
        
        plant.base.t(end+1:end+nIter) = plant.base.t(end) + t;
        plant.base.c(:,end+1:end+nIter) = c;
        plant.base.cost(end+1:end+nIter) = cost;
        plant.base.cons(:,end+1:end+nIter) = cons;
        
        for nO = 1:numel(plant.offset)
            [~, t, c, cost, cons] = ...
                CSTRdyn(plant.offset(nO).u(end,:), plant.offset(nO).c(:,end), plant, tau);
            
            nIter = numel(t);
            
            plant.offset(nO).t(end+1:end+nIter) = plant.offset(nO).t(end) + t;
            plant.offset(nO).c(:,end+1:end+nIter) = c;
            plant.offset(nO).cost(end+1:end+nIter) = cost;
            plant.offset(nO).cons(:,end+1:end+nIter) = cons;
            
            costGrad(nO) = (plant.offset(nO).cost(end) - plant.base.cost(end))./du;
            consGrad(:,nO) = (plant.offset(nO).cons(:,end) - plant.base.cons(:,end))./du;
        end
        
        
    end

end