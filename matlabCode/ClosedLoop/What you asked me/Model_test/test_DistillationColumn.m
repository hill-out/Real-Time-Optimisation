clc
close all    
clear all

disp(' ')
disp('=================================================================')
disp(' Aris PAPASAVVAS                                                 ')
disp(' Test of the Distillation column                                 ')
disp('=================================================================')

addpath('Plant')
addpath('Functions')

%% Initialization
PBstruct = ProblemStructure();
x0 = PBstruct.x0;
x0_DC = x0(12:end);
fsolve_options = PBstruct.fsolve_options;

Theta_P = UncertainParametersStructure('Plant');

%% Inputs definition 
F_E = 41843.443;  % [kg/h]
T_E = 310.928;    % [K]
A_E = 0.126;      % [kgA/kg] - source [1]
B_E = 0.399;      % [kgB/kg] - source [1]
C_E = 0.026;      % [kgC/kg] - source [1]
E_E = 0.361;      % [kgE/kg] - source [1]
G_E = 0;          % [kgG/kg] - source [1]
P_E = 0.088;      % [kgP/kg] - source [1]

% Link the inputs definition to the inputs:
u_DC(1) = F_E;
u_DC(2) = A_E;
u_DC(3) = B_E;
u_DC(4) = C_E;
u_DC(5) = E_E;
u_DC(6) = G_E;
u_DC(7) = P_E;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    2.-Steady state / Dynamic simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x_DC_sol = fsolve(@(x_DC) DistillationCollumnDyn(0, x_DC, u_DC, Theta_P),x0_DC, fsolve_options);

%%
disp(['A_S   Paper: 0.133   Model: ', num2str(x_DC_sol(6))])
disp(['B_S   Paper: 0.421   Model: ', num2str(x_DC_sol(7))])
disp(['C_S   Paper: 0.027   Model: ', num2str(x_DC_sol(8))])
disp(['E_S   Paper: 0.381   Model: ', num2str(x_DC_sol(9))])
disp(['G_S   Paper: 0       Model: ', num2str(0)])
disp(['P_S   Paper: 0.038   Model: ', num2str(x_DC_sol(10))])

sum(x_DC_sol(6:10))





disp('=================================================================')
disp('                             End                                 ')
disp('=================================================================')