clc
close all    
clear all

disp(' ')
disp('=================================================================')
disp(' Aris PAPASAVVAS                                                 ')
disp('=================================================================')

% c = [F_A, F_B, T_R, alpha]
%
% u = [FF_A; FF_B; FF_C; FF_E; FF_G; FF_P; T_R; alpha];
%
% x = [A_R; B_R; C_R; E_R; G_R; P_R; ...                                   // outputs WOreactorDyn          (1:6)
%      A_E; B_E; C_E; E_E; P_E; ...                                        // outputs DecanterDyn           (7:11)
%      A_C; B_C; C_C; E_C; G_C; P_C; A_S; B_S; C_S; E_S; G_S; P_S;...      // outputs DistillationColumnDyn (12:23)
%

addpath('Plant')
addpath('Functions')


%% Initialization
PBstruct = ProblemStructure();
x0 = PBstruct.x0;
fsolve_options = PBstruct.fsolve_options;

Theta_P = UncertainParametersStructure('Plant');

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Inputs definition
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
F_A = 6577.089366;  % [kg/h]  - source [1]
T_A = 294.26111111; % [°K]    - source [1]
A_A = 1;
%
F_B = 15127.305542; % [kg/h]  - source [1]
T_B = 294.26111111; % [°K]    - source [1]
B_B = 1;
%
F_L = 21822.782516; % [kg/h]   - source [1]
T_L = 310.92777778; % [°K]     - source [1]
A_L = 0.133;        % [kgA/kg] - source [1]
B_L = 0.421;        % [kgB/kg] - source [1]
C_L = 0.027;        % [kgC/kg] - source [1]
E_L = 0.381;        % [kgE/kg] - source [1]
G_L = 0;            % [kgG/kg] - source [1]
P_L = 0.038;        % [kgP/kg] - source [1]

FF_A  = F_A + F_L*A_L;
FF_B  = F_B + F_L*B_L;
FF_C  = F_L*C_L;
FF_E  = F_L*E_L;
FF_G  = F_L*G_L;
FF_P  = F_L*P_L;
T_R   = 355.372;       % [K]     - source [1]
alpha = 0.5499;        % [-]     - source [1]

u = [FF_A; FF_B; FF_C; FF_E; FF_G; FF_P; T_R; alpha];

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x_sol = fsolve(@(x) OpenLoopSystem(0, x, u, Theta_P), x0, fsolve_options);

F_P_demand = 2160;
Theta_M = UncertainParametersStructure('Plant');
[g, phi_u, F_P] = ux2gphi(u, x_sol, F_P_demand, Theta_M);


%%
disp(['A_R   Paper: 0.121   Model: ', num2str(x_sol(1))])
disp(['B_R   Paper: 0.383   Model: ', num2str(x_sol(2))])
disp(['C_R   Paper: 0.025   Model: ', num2str(x_sol(3))])
disp(['E_R   Paper: 0.346   Model: ', num2str(x_sol(4))])
disp(['G_R   Paper: 0.039   Model: ', num2str(x_sol(5))])
disp(['P_R   Paper: 0.085   Model: ', num2str(x_sol(6))])
disp(' ')
sum(x_sol(1:6))
disp(' ')
%%
disp(['A_E   Paper: 0.126   Model: ', num2str(x_sol(7))])
disp(['B_E   Paper: 0.399   Model: ', num2str(x_sol(8))])
disp(['C_E   Paper: 0.026   Model: ', num2str(x_sol(9))])
disp(['E_E   Paper: 0.361   Model: ', num2str(x_sol(10))])
disp(['P_E   Paper: 0.088   Model: ', num2str(x_sol(11))])
disp(' ')
sum(x_sol(7:11))
disp(' ')
%%
disp(['A_S   Paper: 0.133   Model: ', num2str(x_sol(17))])
disp(['B_S   Paper: 0.421   Model: ', num2str(x_sol(18))])
disp(['C_S   Paper: 0.027   Model: ', num2str(x_sol(19))])
disp(['E_S   Paper: 0.381   Model: ', num2str(x_sol(20))])
disp(['G_S   Paper: 0       Model: ', num2str(0)])
disp(['P_S   Paper: 0.038   Model: ', num2str(x_sol(21))])
disp(' ')
sum(x_sol(17:21))
disp(' ')

disp(['F_P     closed: 2160.4512        Model: ', num2str(F_P)])
disp(['phi_u   closed: -13734081.3737   Model: ', num2str(phi_u)])


%%



