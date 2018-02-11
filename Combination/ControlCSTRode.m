function [dy] = ControlCSTRode(t, y, r, r0, rK)
% Takes the time and setpoints of the plant and calculates the outputs
% -------------------------------------------------------------------------
% t         Current time
% y         Plant out [F_Ap, F_Bp, T_R, X_Ap, X_Bp, X_Cp, X_Pp, X_Ep, X_Gp]
% r         Controller setpoints [F_Bs, X_As]
% r0        Controller zeropoints [F_Bs0, X_As0]
% rK        Controller Gain [K_FBs, K_XAs]
%
% dy        time derivative y
% -------------------------------------------------------------------------

dy    = zeros(size(y)); %initialise

F_Ain = y(1);
F_Bin = y(2);
T     = y(3);
X     = y(4:end);

F_Bs  = r(1)+2;
X_As  = r(2)*0.75;
F_Bs0 = r0(1);
X_As0 = r0(2);
K_FBs = rK(1);
K_XAs = rK(2);

% controller
dF_B = K_FBs*(F_Bs - F_Bin);
dF_A = dF_B/2.4;

dT_R = K_XAs*(X_As - X(1));

% CSTR
M = 2105; %kg
F = F_Ain + F_Bin; %kg/s
k1 = 1.66e6*exp(-5.543e4/(8.314*(T+273.15))); %1/s
k2 = 7.21e8*exp(-6.928e4/(8.314*(T+273.15))); %1/s
k3 = 2.67e12*exp(-9.238e4/(8.314*(T+273.15))); %1/s

R1 = k1*X(1)*X(2);
R2 = k2*X(2)*X(3);
R3 = k3*X(3)*X(4);

dF(1) = F_Ain - F*X(1) -   M*R1                    ; %A
dF(2) = F_Bin - F*X(2) -   M*R1 -   M*R2           ; %B
dF(3) = 0     - F*X(3) + 2*M*R1 - 2*M*R2 -     M*R3; %C
dF(4) = 0     - F*X(4)          +   M*R2 - 0.5*M*R3; %P
dF(5) = 0     - F*X(5)          + 2*M*R2           ; %E
dF(6) = 0     - F*X(6)                   + 1.5*M*R3; %G

% output
dy(1:3) = [dF_A, dF_B, dT_R];
dy(4:end) = dF/F;

end

