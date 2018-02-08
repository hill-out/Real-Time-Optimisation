clc
close all    
clear all

disp(' ')
disp('=================================================================')
disp(' Aris PAPASAVVAS                                                 ')
disp('=================================================================')
tic

addpath('Plant')
addpath('Functions')

%% Initialization
PBstruct = ProblemStructure();
x0 = PBstruct.x0;
x0_R = x0(1:6);
fsolve_options = PBstruct.fsolve_options;

Theta_P = UncertainParametersStructure('Plant');

%% Inputs definition 
F_A = 10/2.4*3600; %6577.089366;  % [kg/h]  - source [1]
T_A = 294.26111111; % [°K]    - source [1]
A_A = 1;
%
F_B = 10*3600; %15127.305542; % [kg/h]  - source [1]
T_B = 294.26111111; % [°K]    - source [1]
B_B = 1;
%
F_L = 0; %21822.782516; % [kg/h]   - source [1]
T_L = 310.92777778; % [°K]     - source [1]
A_L = 0.133;        % [kgA/kg] - source [1]
B_L = 0.421;        % [kgB/kg] - source [1]
C_L = 0.027;        % [kgC/kg] - source [1]
E_L = 0.381;        % [kgE/kg] - source [1]
G_L = 0;            % [kgG/kg] - source [1]
P_L = 0.038;        % [kgP/kg] - source [1]



%%

F_R = F_A + F_B + F_L;
X_A = (F_L*A_L + F_A) / F_R;
X_B = (F_L*B_L + F_B) / F_R;
X_C = F_L*C_L / F_R;
X_E = F_L*E_L / F_R;
X_G = F_L*G_L / F_R;
X_P = F_L*P_L / F_R;

T_R =  355.372; % 180°F

u_R = [ F_R;
        X_A;
        X_B;
        X_C;
        X_E;
        X_G;
        X_P;
        T_R]; 

x_R_sol = fsolve(@(x_R) WOreactorDyn(0, x_R, u_R, Theta_P),x0_R, fsolve_options);

%%
disp(['A_R   Paper: 0.121   Model: ', num2str(x_R_sol(1))])
disp(['B_R   Paper: 0.383   Model: ', num2str(x_R_sol(2))])
disp(['C_R   Paper: 0.025   Model: ', num2str(x_R_sol(3))])
disp(['E_R   Paper: 0.346   Model: ', num2str(x_R_sol(4))])
disp(['G_R   Paper: 0.039   Model: ', num2str(x_R_sol(5))])
disp(['P_R   Paper: 0.085   Model: ', num2str(x_R_sol(6))])

%% 
disp(' ')
toc
disp(' ')
disp('=================================================================')
disp('                            end                                  ')
disp('=================================================================')





