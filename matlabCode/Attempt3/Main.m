%% main script
%
% main steps:
% 
% 1. gets parametrs
% 2. sets up functions (cost, cons)
% 3. finds opt of model
% 4. defines conv. approx.
% 5. runs RTO
% 6. plots figures

%% 1. Get parameters

allP = parameterDefiner('mpc');
model = allP{1};
plant = allP{2};
cons = allP{3};

%% 2. Set up functions (cost, cons)

% model functions
model.cost = @(u, c)(costFun(u, c, model));
model.cons = @(~, c)(consFun([], c, model));

% plant functions
plant.cost = @(u, c)(costFun(u, c, plant));
plant.cons = @(~, c)(consFun([], c, plant));

%% 3. Find opt of model

[model.uOpt, model.cOpt, model.costOpt, model.consOpt] = CSTRopt(model,[14,14]);


%% 4. Define Convex Approximation

[model.convCost, model.convCons] = convApprox(model);

%$ 5. Run RTO