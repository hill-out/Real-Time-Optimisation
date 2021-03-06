function [Xsol] = openPlant(u, xGuess, dtheta)
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

% inputs
F_Ain = u(1);
F_Bin = u(2);
T     = u(3);

% uncertain model parameters
k_0 = [1.66e6, 7.21e8, 2.67e12]; %1/s
E = [5.543e4, 6.928e4, 9.238e4]; %1/K
M = 2105; %kg

if nargin == 3
    E   = E + dtheta(1:3);
end

F = F_Ain + F_Bin; %kg/s
k1 = k_0(1)*exp(-E(1)/(8.314*(T+273.15))); %1/s
k2 = k_0(2)*exp(-E(2)/(8.314*(T+273.15))); %1/s
k3 = k_0(3)*exp(-E(3)/(8.314*(T+273.15))); %1/s

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