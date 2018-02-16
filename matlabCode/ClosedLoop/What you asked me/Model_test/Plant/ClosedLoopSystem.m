function dx = ClosedLoopSystem(t, x, c, Parameters)

% c = [F_A, F_B, T_R, alpha]
%
% u = [FF_A; FF_B; FF_C; FF_E; FF_G; FF_P; T_R; alpha];
%
% x = [A_R; B_R; C_R; E_R; G_R; P_R; ...                                   // outputs WOreactorDyn          (1:6)
%      A_E; B_E; C_E; E_E; P_E; ...                                        // outputs DecanterDyn           (7:11)
%      A_C; B_C; C_C; E_C; P_C; A_S; B_S; C_S; E_S; P_S;...      // outputs DistillationColumnDyn (12:23)
%

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Inputs definition
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% c
F_A   = c(1);   
F_B   = c(2);   
T_R   = c(3);   
alpha = c(4);

% x
A_R    = x(1);
B_R    = x(2);
C_R    = x(3);
E_R    = x(4);
G_R    = x(5);
P_R    = x(6);

A_E    = x(7);
B_E    = x(8);
C_E    = x(9);
E_E    = x(10);
P_E    = x(11);

A_S    = x(17);
B_S    = x(18);
C_S    = x(19);
E_S    = x(20);
P_S    = x(21);

eff = Parameters.eff;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%                   *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*

beta = alpha*(1-G_R)*(1-P_E+eff*E_E);
F_L = (F_A+F_B)*beta/(1-beta);


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Mixer
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

F_R = F_A + F_B + F_L;

X_A = (F_L*A_S + F_A) / F_R;
X_B = (F_L*B_S + F_B) / F_R;
X_C = F_L*C_S / F_R;
X_E = F_L*E_S / F_R;
X_G = 0;
X_P = F_L*P_S / F_R;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                  %
                                  %
                                  %
                               %  %  %   
                                % % %
                                  %
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  WO reactor (WOr)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x_R  = x(1:6);

u_R = [ F_R;
        X_A;
        X_B;
        X_C;
        X_E;
        X_G;
        X_P;
        T_R];
    
dx_R = WOreactorDyn(t, x_R, u_R, Parameters);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                  %
                                  %
                                  %
                               %  %  %   
                                % % %
                                  %
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Decanter (D)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x_D = x(7:11);

u_D = zeros(7,1);
u_D(1) = F_R;
u_D(2) = A_R;
u_D(3) = B_R;
u_D(4) = C_R;
u_D(5) = E_R;
u_D(6) = G_R;
u_D(7) = P_R;

dx_D = DecanterDyn(t, x_D, u_D);

F_G = F_R*G_R;
F_E = F_R - F_G;
G_E = 0;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                  %
                                  %
                                  %
                               %  %  %
                                % % %
                                  %
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Distillation Column (DC)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x_DC = x(12:21);

u_DC    = zeros(8,1);
u_DC(1) = F_E;
u_DC(2) = A_E;
u_DC(3) = B_E;
u_DC(4) = C_E;
u_DC(5) = E_E;
u_DC(6) = G_E;
u_DC(7) = P_E;

dx_DC = DistillationCollumnDyn(t, x_DC, u_DC, Parameters);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                  %
                                  %
                                  %
                               %  %  %
                                % % %
                                  %
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Recycling (R)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Already taken into account. See lines 49-50
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%                   *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*


dx = [dx_R; dx_D; dx_DC];


end