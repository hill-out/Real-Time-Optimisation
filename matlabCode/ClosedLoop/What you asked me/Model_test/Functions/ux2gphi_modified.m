function  [g, phi_u, phi_c, F_P] = ux2gphi_modified(u, c, x, F_P_demand, Parameters)

% x = [A_R; B_R; C_R; E_R; G_R; P_R; ...                                   // outputs WOreactorDyn          (1:6)
%      A_E; B_E; C_E; E_E; P_E; ...                                        // outputs DecanterDyn           (7:11)
%      A_C; B_C; C_C; E_C; G_C; P_C;     A_S; B_S; C_S; E_S; G_S; P_S];    // outputs DistillationColumnDyn (12:23)
%
% u = [F_A; F_B; F_C; F_E; F_P; F_G; T_R; alpha]; 
%
% c = [F_A; F_B; T_R; alpha];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
F     = sum(u(1:6));
alpha = u(8);

F_A = c(1);
F_B = c(2);

eff = Parameters.eff;

G_R = x(5);
E_E = x(10);
P_E = x(11);

F_R = F;
F_E = F_R*(1-G_R);
F_G = F_R*G_R;
F_P = F_E*(P_E-eff*E_E); 

F_S = F_E-F_P;

% Profit function to maxmize

phi_u = ( 0.6614*F_P ...
          + 0.0150*(1-alpha)*F_S ...
          - 0.0220*F_G ...
         ) * 8400  ...                          
          - 4.8943*F_R - 1.0416e+03*(0.6614*F_P+0.0150*(1-alpha)*F_S); 

phi_c =   (- 0.0441*F_A ...
           - 0.0661*F_B ...
         ) * 8400;     

% Costs function to minimize  
phi_u = -phi_u;
phi_c = -phi_c;

g= F_P-F_P_demand;

end