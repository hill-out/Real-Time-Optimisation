% script for running steady state real time optimisation

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

% Model reaction
model_reactOrder = [-1, -1,  1,  0;
                     0, -2,  0,  1];
model_kVal = [0.75; 1.5];       % L/(mol min)
model_cIn = [2, 1.5, 0, 0];     % mol/L
model_V = 500;                  % L
model_H = [3.5; 1.5];           % kcal/mol

%Plant reation
plant_reactOrder = model_reactOrder;
plant_kVal = [1.4; 0.4];        % L/(mol min)
plant_cIn = [2.5, 1.5, 0, 0];   % mol/L
plant_V = model_V;              % L
plant_H = model_H;              % kcal/mol

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

% Cost
w = 0.004; %(mol min)/L^2
J = @(c0,c,u)(c(3)^2*(u(1)+u(2))^2/(u(1)*c0(1))-w*(u(1)^2+u(2)^2));
phi = @(c0,c,u)(-J(c0,c,u));

% Constraint 1
model_Qmax = 110; %kcal(/min)?
plant_Qmax = 110; %kcal(/min)?

model_Q = @(c0,c,u)(model_V*(model_kVal(1)*c(1)*c(2)*model_H(1)+model_kVal(2)*c(2)^2*model_H(2)));
plant_Q = @(c0,c,u)(plant_V*(plant_kVal(1)*c(1)*c(2)*plant_H(1)+plant_kVal(2)*c(2)^2*plant_H(2)));

model_G{1} = @(c0,c,u)(model_Q(c0,c,u)/model_Qmax - 1);
plant_G{1} = @(c0,c,u)(plant_Q(c0,c,u)/plant_Qmax - 1);

% Constraint 2
model_Dmax = 0.1;
plant_Dmax = 0.1;

model_D = @(c0,c)(c(4)/sum(c));
plant_D = model_D;

model_G{2} = @(c0,c,u)(model_D(c0,c)/model_Dmax - 1);
plant_G{2} = @(c0,c,u)(plant_D(c0,c)/plant_Dmax - 1);

% Initial conditions
u0 = [10, 10, 0, 0]; % L/min

model_opt = CSTR_Opt(model_reactOrder, model_kVal, model_cIn, u0, model_V, phi, model_G);
model_cSol = CSTR(model_reactOrder, model_kVal, model_cIn, [model_opt(1), model_opt(2), 0, 0], model_V);

convexApprox(phi, [1, 1], [49, 49], 1, model_opt, model_reactOrder, model_kVal, model_cIn, model_V);




