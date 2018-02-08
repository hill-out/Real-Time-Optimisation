clc
close all    
clear all

disp(' ')
disp('=================================================================')
disp(' Aris PAPASAVVAS:                                                ')
disp('=================================================================')
tic

addpath('Plant')
addpath('Functions')

%% Initialization
PBstruct = ProblemStructure();
x0 = PBstruct.x0;
x0_D = x0(7:11);
fsolve_options = PBstruct.fsolve_options;


%% Inputs definition 
F_R = 43527;  % [kg/h]
T_R = 360;    % [K]
A_R = 0.121;  % [kgA/kgH]
B_R = 0.383;  % [kgB/kgH]
C_R = 0.025;  % [kgC/kgH]
E_R = 0.346;  % [kgE/kgH]
G_R = 0.039;  % [kgG/kgH]
P_R = 0.085;  % [kgP/kgH]

% Link the inputs definition to the inputs:
u_D(1) = F_R;
u_D(2) = A_R;
u_D(3) = B_R;
u_D(4) = C_R;
u_D(5) = E_R;
u_D(6) = G_R;
u_D(7) = P_R;


x_D_sol = fsolve(@(x_D) DecanterDyn(0, x_D, u_D),x0_D, fsolve_options);

%%
disp(['A_E   Paper: 0.126   Model: ', num2str(x_D_sol(1))])
disp(['B_E   Paper: 0.399   Model: ', num2str(x_D_sol(2))])
disp(['C_E   Paper: 0.026   Model: ', num2str(x_D_sol(3))])
disp(['E_E   Paper: 0.361   Model: ', num2str(x_D_sol(4))])
disp(['P_E   Paper: 0.088   Model: ', num2str(x_D_sol(5))])

sum(x_D_sol)
%% 
disp(' ')
toc
disp(' ')
disp('=================================================================')
disp('                            end                                  ')
disp('=================================================================')





