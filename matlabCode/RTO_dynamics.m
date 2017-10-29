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






