function [u, g, phi_u, phi_c, phi_p, F_P] = cx2ugphi(c, x, F_P_demand, Parameters)

% x = [A_R; B_R; C_R; E_R; G_R; P_R; ...                                   // outputs WOreactorDyn          (1:6)
%      A_E; B_E; C_E; E_E; P_E; ...                                        // outputs DecanterDyn           (7:11)
%      A_C; B_C; C_C; E_C; G_C; P_C;  A_S; B_S; C_S; E_S; G_S; P_S;...     // outputs DistillationColumnDyn (12:23)

% c = [F_A; F_B; T_R; alpha];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

E_R = x(4);
G_R = x(5);
P_R = x(6);


E_E = x(10);
P_E = x(11);

A_S = x(17);
B_S = x(18);
C_S = x(19);
E_S = x(20);
P_S = x(21);

F_A   = c(1);
F_B   = c(2);
T_R   = c(3);
alpha = c(4);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

beta = alpha*(1-G_R)*(1-P_E+0.1*E_E);
F_L = (F_A+F_B)*beta/(1-beta);

F_R = F_A + F_B + F_L; 

F_E = F_R*(1-G_R);

eff = Parameters.eff;

F_G = F_R*G_R;
F_P = (P_E-eff*E_E)*F_E;  
F_S = F_E - F_P;


FF_A = F_A+F_L*A_S;
FF_B = F_B+F_L*B_S;
FF_C = F_L*C_S;
FF_E = F_L*E_S;
FF_P = F_L*P_S;
FF_G = 0;


% alpha = c(4);

u = [FF_A; FF_B; FF_C; FF_E; FF_G; FF_P; T_R; alpha]; 

g = F_P-F_P_demand;
 
phi_u = ( 0.6614*F_P ...
          + 0.0150*(1-alpha)*F_S ...
          - 0.0220*F_G ...
         ) * 8400 ...
          - 4.8943*F_R  - 1.0416e+03*(0.6614*F_P+0.0150*(1-alpha)*F_S); 

phi_c =   (- 0.0441*F_A ...
           - 0.0661*F_B ... 0.0661
         ) * 8400;     
     
phi_u = -phi_u;
phi_c = -phi_c;
     
     
phi_p =  phi_u + phi_c;

end