function dx_R = WOreactorDyn(t, x_R, u_R, Parameters)

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  1.- Parameters
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reactor design 
% Vr    = 2.1088;               % [m3]    Reactor volume
% rho_R = 800.923168698;     % [kg/m^3] Density of the reactant mixture
% %
W = 2105; %Vr*rho_R;   % [kg]
%
A1 = Parameters.A1; % [1/h]
B1 = Parameters.B1; % [°K]
A2 = Parameters.A2; % [1/h]
B2 = Parameters.B2; % [°K]
A3 = Parameters.A3; % [1/h]
B3 = Parameters.B3; % [°K]
A4 = Parameters.A4; % [1/h]
B4 = Parameters.B4; % [°K]
A5 = Parameters.A5; % [1/h]
B5 = Parameters.B5; % [°K]

MA = 100;
MB = 100;
MC = 200;
ME = 200;
MG = 300;
MP = 100;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  2.- Nomenclature
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

F   = u_R(1);
X_A = u_R(2);
X_B = u_R(3);
X_C = u_R(4);
X_E = u_R(5);
X_G = u_R(6);
X_P = u_R(7);
T_R = u_R(8);

A_R = x_R(1);
B_R = x_R(2);
C_R = x_R(3);
E_R = x_R(4);
G_R = x_R(5);
P_R = x_R(6);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  3.- Model
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
k1 = A1*exp(-B1/T_R);
k2 = A2*exp(-B2/T_R);
k3 = A3*exp(-B3/T_R);
k4 = A4*exp(-B4/T_R);
k5 = A5*exp(-B5/T_R);
%
% A +  B -> C      (1)
% C +  B -> P + E  (2)
% P +  C -> G      (3)
% A + 2B -> P + E  (4)
% A +  B + P -> G      (5)
%  
r1 = k1*A_R*B_R;
r2 = k2*B_R*C_R;
r3 = k3*C_R*P_R;
r4 = k4*A_R*B_R^2;
r5 = k5*A_R*B_R*P_R;
%
dx_R = zeros(6,1);
dx_R(1) = F*X_A/W - F*A_R/W - r1                              - r4         - r5;
dx_R(2) = F*X_B/W - F*B_R/W - MB/MA*r1 - r2                   - 2*MB/MA*r4 - MB/MA*r5;
dx_R(3) = F*X_C/W - F*C_R/W + MC/MA*r1 - MC/MB*r2 - r3;
dx_R(4) = F*X_E/W - F*E_R/W + ME/MP*r2                        + ME/MA*r4;
dx_R(5) = F*X_G/W - F*G_R/W                       + MG/MC*r3               + MG/MA*r5;
dx_R(6) = F*X_P/W - F*P_R/W            + MP/MB*r2 - MP/MC*r3  + MP/MA*r4   - MP/MA*r5;



