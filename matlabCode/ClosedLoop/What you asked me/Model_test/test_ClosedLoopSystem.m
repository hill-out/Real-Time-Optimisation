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
%      A_C; B_C; C_C; E_C; P_C; A_S; B_S; C_S; E_S; P_S;...                // outputs DistillationColumnDyn (12:21)
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
% Manipulable variables
F_A   = 6577.089366;  % [kg/h]  - source [1]
F_B   = 15127.305542; % [kg/h]  - source [1]
T_R   = 355.372;      % [K]     - source [1]
alpha = 0.5499;       % [-]     - source [1]

c = [F_A; F_B; T_R; alpha];

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x_sol = fsolve(@(x) ClosedLoopSystem(0, x, c, Theta_P), x0, fsolve_options);


F_P_demand = 2160;
Theta_P = UncertainParametersStructure('Plant');
[u, g, phi_u, phi_c, phi_p, F_P] = cx2ugphi(c, x_sol, F_P_demand, Theta_P);





%%
disp(['A_R   Paper: 0.121   Model: ', num2str(x_sol(1))])
disp(['B_R   Paper: 0.383   Model: ', num2str(x_sol(2))])
disp(['C_R   Paper: 0.025   Model: ', num2str(x_sol(3))])
disp(['E_R   Paper: 0.346   Model: ', num2str(x_sol(4))])
disp(['G_R   Paper: 0.039   Model: ', num2str(x_sol(5))])
disp(['P_R   Paper: 0.085   Model: ', num2str(x_sol(6))])
% disp(' ')
% sum(x_sol(1:6))
disp(' ')
%%
disp(['A_E   Paper: 0.126   Model: ', num2str(x_sol(7))])
disp(['B_E   Paper: 0.399   Model: ', num2str(x_sol(8))])
disp(['C_E   Paper: 0.026   Model: ', num2str(x_sol(9))])
disp(['E_E   Paper: 0.361   Model: ', num2str(x_sol(10))])
disp(['P_E   Paper: 0.088   Model: ', num2str(x_sol(11))])
% disp(' ')
% sum(x_sol(7:11))
disp(' ')
%%
disp(['A_S   Paper: 0.133   Model: ', num2str(x_sol(17))])
disp(['B_S   Paper: 0.421   Model: ', num2str(x_sol(18))])
disp(['C_S   Paper: 0.027   Model: ', num2str(x_sol(19))])
disp(['E_S   Paper: 0.381   Model: ', num2str(x_sol(20))])
disp(['P_S   Paper: 0.038   Model: ', num2str(x_sol(21))])
% disp(' ')
% sum(x_sol(17:21))
disp(' ')
%%

disp(['F_P     Model: 2155.7513        closed: ', num2str(F_P)])
disp(['phi_u   Model: -13707891.2863   closed: ', num2str(phi_u)])
disp(['phi_p   ', num2str(phi_p)])

disp(' ')
disp('=================================================================')
disp('                             End                                 ')
disp('=================================================================')