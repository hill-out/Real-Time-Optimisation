function dx_DC = DistillationCollumnDyn(t, x_DC, u, Parameters)

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  1.- Parameters
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Distillation column design parameters
Vc    = 13.3882;       % [m3]
rho_C = 800.923168698; % [kg/m^3] Density ???
Wc    = Vc*rho_C;

Vs    = 5.7087;        % [m3]
rho_S = 800.923168698; % [kg/m^3] Density ???
Ws    = Vs*rho_S;

eff = Parameters.eff;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  2.- Nomenclature
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
F_E = u(1);
A_E = u(2);
B_E = u(3);
C_E = u(4);
E_E = u(5);
G_E = u(6);
P_E = u(7);

A_C = x_DC(1);
B_C = x_DC(2);
C_C = x_DC(3);
E_C = x_DC(4);
P_C = x_DC(5);

A_S = x_DC(6);
B_S = x_DC(7);
C_S = x_DC(8);
E_S = x_DC(9);
P_S = x_DC(10);

F_P = (P_E-eff*E_E)*F_E;
F_S = F_E - F_P;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  3.- Model
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dA_C_dt = F_E/(Wc)*A_E - (2*F_E-F_P)/(Wc)*A_C + F_E/(Wc)*A_S;
dB_C_dt = F_E/(Wc)*B_E - (2*F_E-F_P)/(Wc)*B_C + F_E/(Wc)*B_S;
dC_C_dt = F_E/(Wc)*C_E - (2*F_E-F_P)/(Wc)*C_C + F_E/(Wc)*C_S;
dE_C_dt = F_E/(Wc)*E_E - (2*F_E-F_P)/(Wc)*E_C + F_E/(Wc)*E_S;
dP_C_dt = eff*F_E/(Wc)*E_E - (2*F_E-F_P)/(Wc)*P_C + F_E/(Wc)*P_S;

dA_S_dt = (2*F_E-F_P)/(Ws)*A_C - F_S/(Ws)*A_S - F_E/(Ws)*A_S;
dB_S_dt = (2*F_E-F_P)/(Ws)*B_C - F_S/(Ws)*B_S - F_E/(Ws)*B_S;
dC_S_dt = (2*F_E-F_P)/(Ws)*C_C - F_S/(Ws)*C_S - F_E/(Ws)*C_S;
dE_S_dt = (2*F_E-F_P)/(Ws)*E_C - F_S/(Ws)*E_S - F_E/(Ws)*E_S;
dP_S_dt = (2*F_E-F_P)/(Ws)*P_C - F_S/(Ws)*P_S - F_E/(Ws)*P_S;

dx_DC = [dA_C_dt;
         dB_C_dt;
         dC_C_dt;
         dE_C_dt;
         dP_C_dt;
         dA_S_dt;
         dB_S_dt;
         dC_S_dt;
         dE_S_dt;
         dP_S_dt];

end