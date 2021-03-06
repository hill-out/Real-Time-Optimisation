function dx_D = DecanterDyn(t, x_D, u_D)

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  1.- Parameters
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Decanter design parameters
VE    = 5.23862;
rho_R = 800.923168698;     % [kg/m^3] Density of the reactant mixture
W = rho_R*VE;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  2.- Nomenclature
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
F_R = u_D(1);
A_R = u_D(2);
B_R = u_D(3);
C_R = u_D(4);
E_R = u_D(5);
G_R = u_D(6);
P_R = u_D(7);

F_G = F_R*G_R;
F_E = F_R-F_G;

A_E = x_D(1);
B_E = x_D(2);
C_E = x_D(3);
E_E = x_D(4);
P_E = x_D(5);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  3.- Model
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dA_E_dt = F_R/W*A_R-F_E/W*A_E; 
dB_E_dt = F_R/W*B_R-F_E/W*B_E;
dC_E_dt = F_R/W*C_R-F_E/W*C_E;
dE_E_dt = F_R/W*E_R-F_E/W*E_E;
dP_E_dt = F_R/W*P_R-F_E/W*P_E;

dx_D = [dA_E_dt;
        dB_E_dt;
        dC_E_dt;
        dE_E_dt;
        dP_E_dt];
end

