%RTO_Dynamics

clear

% #### /Reation Data ####
%
%  A +  B ->  C         (R1)
%      2B ->      D     (R2)
% liquid, constant density reation
% Constant volume CSTR at 500L
% k1 = 0.75, k2 = 1.5                           % L/(mol min)
% cAin = 2, cBin = 1.5, cCin = 0, cDin = 0      % mol/L
% H1 = 3.5, H2 = 1.5                            % kcal/mol
% C is the component of interest
% D is a byproduct
%
% #### Reaction Data/ ####

% #### /Optimisation Parameters ####
%
% model cost function:
% -phi = J = cC^2*(uA + uB)^2/(uA*cAin) - w(uA^2+uB^2)
%
% where: w = 0.004 (mol min)/L^2
%
% model constraints
% G1 = Q/Qmax - 1 = 0
% G2 = D/Dmax - 1 = 0
%
% where: Qmax = 110 kcal(/min)?
%        Dmax = 0.1
%
% #### Optimisation Parameters/ ####

% Model reaction
model_reactOrder = [-1, -1,  1,  0;
                     0, -2,  0,  1];
model_kVal = [0.75; 1.5];       % L/(mol min)
model_cIn = [2, 1.5, 0, 0];     % mol/L
model_V = 500;                  % L
model_H = [3.5; 1.5];           % kcal/mol

% Plant reation
plant_reactOrder = model_reactOrder;
plant_kVal = [1.4; 0.4];        % L/(mol min)
plant_cIn = [2.5, 1.5, 0, 0];   % mol/L
plant_V = model_V;              % L
plant_H = model_H;              % kcal/mol

% Optimum Parameters (plant)
plant_fun = @(x,u,n,t)(dyn_CSTR(plant_reactOrder, plant_kVal, x, plant_cIn, u, plant_V, n, t));

% Cost
w = 0.004;
model_phi = @(u)(CSTR_cost(u, w, model_cIn, model_reactOrder, model_kVal, model_V));
plant_phi = @(u)(CSTR_cost(u, w, plant_cIn, plant_reactOrder, plant_kVal, plant_V));

% Constraint 1
model_Qmax = 110; %kcal(/min)?
plant_Qmax = 110; %kcal(/min)?

model_G{1} = @(u)(CSTR_Qcon(u,model_H,model_Qmax,model_cIn,model_reactOrder,model_kVal,model_V));
plant_G{1} = @(u)(CSTR_Qcon(u,plant_H,plant_Qmax,plant_cIn,plant_reactOrder,plant_kVal,plant_V));

% Constraint 2
model_Dmax = 0.1;
plant_Dmax = 0.1;

model_G{2} = @(u)(CSTR_Dcon(u,model_Dmax,model_cIn,model_reactOrder,model_kVal,model_V));
plant_G{2} = @(u)(CSTR_Dcon(u,plant_Dmax,plant_cIn,plant_reactOrder,plant_kVal,plant_V));

% Initial conditions
u0 = [10, 10, 0, 0]; % L/min

model_opt = CSTR_Opt(u0, model_phi, model_G);
model_cSolOpt = CSTR(model_reactOrder, model_kVal, model_cIn, [model_opt(1), model_opt(2), 0, 0], model_V);

plant_c0 = CSTR(plant_reactOrder, plant_kVal, plant_cIn, [model_opt(1), model_opt(2), 0, 0], plant_V);

% Find opt point
model_funArray = {model_phi,model_G{1},model_G{2}};
model_funOpt = zeros(1,3);
for j = 1:length(model_funArray)
    u = zeros(2,4);
    u(2,:) = model_cSolOpt;
    u(1,[1,2]) = model_opt;
    model_funOpt(j) = model_funArray{j}(u);
end

% Convex approximation parameters
convPara = convexApprox(model_funArray, [1, 1], [50, 50],...
    [2,1,1], model_opt, model_funOpt, model_reactOrder, model_kVal, model_cIn, model_V);
convPara = [-0.8305,-0.9121,0.08,0.0051,0.0126,0,-0.0648,0.0857];
% Set-up MA opt
conv_phi = @(u)(model_funOpt(1)+convPara([1,2])*(u-model_opt)'+...
    0.5*(u-model_opt)*(convPara(3)*eye(2))*(u-model_opt)');
conv_G{1} = @(u)(model_funOpt(2)+convPara([1,2]+3)*(u-model_opt)');
conv_G{2} = @(u)(model_funOpt(3)+convPara([1,2]+6)*(u-model_opt)');
conv_funArray = {conv_phi, conv_G{:}};

conv_dphi = @(u)(convPara([1,2])'+(convPara(3))*(u-model_opt)');
conv_dG{1} = @(u)(convPara([1,2]+3)');
conv_dG{2} = @(u)(convPara([1,2]+6)');
conv_dfunArray = {conv_dphi, conv_dG{:}};

% Run plant optimisation
dynamicRTO(conv_funArray,plant_fun,0.2,model_opt,plant_c0,5);





