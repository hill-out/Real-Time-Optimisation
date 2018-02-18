function [Xsol] = CSTRplant(u, xGuess)
% calculates the concentrations of the plant symbolically
% -------------------------------------------------------
% u         Inputs
% xGuess    Guess mass fractions
% 
% Xsol      Mass fractions
% -------------------------------------------------------

if nargin<2
    xGuess = zeros(1,6);
end

F_Ain = u(1);
F_Bin = u(2);
T     = u(3);
E = [8050, 12500];
M = 2105; %kg
F = F_Ain + F_Bin; %kg/s
k1 = 1.66e6*exp(-5.543e4/(8.314*(T+273.15))); %1/s
k2 = 7.21e8*exp(-6.928e4/(8.314*(T+273.15))); %1/s
k3 = 2.67e12*exp(-9.238e4/(8.314*(T+273.15))); %1/s

%% mass balances
optionX = optimoptions('fsolve','Display','off');

Xsol = fsolve(@massBalance, xGuess, optionX);

    function [dF] = massBalance(X)
        
        R1 = k1*X(1)*X(2);
        R2 = k2*X(2)*X(3);
        R3 = k3*X(3)*X(4);
        %       IN    - OUT    +GEN/-CON R1   +GEN/-CON R2
        dF(1) = F_Ain - F*X(1) -   M*R1                    ; %A
        dF(2) = F_Bin - F*X(2) -   M*R1 -   M*R2           ; %B
        dF(3) = 0     - F*X(3) + 2*M*R1 - 2*M*R2 -     M*R3; %C
        dF(4) = 0     - F*X(4)          +   M*R2 - 0.5*M*R3; %P
        dF(5) = 0     - F*X(5)          + 2*M*R2           ; %E
        dF(6) = 0     - F*X(6)                   + 1.5*M*R3; %G
    end
end