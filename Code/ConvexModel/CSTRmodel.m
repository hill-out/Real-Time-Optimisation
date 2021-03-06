function [Xsol] = CSTRmodel(u, xGuess)
% Symbolically solves the model 2-reaction CSTR
% ---------------------------------------------
% u         Inputs
% xGuess    Guess of mass fractions (for speed)
% 
% Xsol      Mass Fraction [A,B,C,P,E,G]
% ---------------------------------------------

if nargin<2
    xGuess = zeros(1,6);
end

F_Ain = u(1);
F_Bin = u(2);
T     = u(3);
E = [8050, 12500];
M = 2105; %kg
F = F_Ain + F_Bin; %kg/s
k1 = 2.189*1e8*exp(-E(1)/(T+273.15)); %1/s
k2 = 4.310*1e13*exp(-E(2)/(T+273.15)); %1/s

%% mass balances
optionX = optimoptions('fsolve','Display','off');

Xsol = fsolve(@massBalance, xGuess, optionX);

    function [dF] = massBalance(X)
        
        R1 = k1*X(1)*X(2)^2;
        R2 = k2*X(1)*X(2)*X(4);
        %       IN    - OUT    +GEN/-CON R1   +GEN/-CON R2
        dF(1) = F_Ain - F*X(1) -   M*R1       -   M*R2; %A
        dF(2) = F_Bin - F*X(2) - 2*M*R1       -   M*R2; %B
        dF(3) = 0     - F*X(3)                        ; %C
        dF(4) = 0     - F*X(4) +   M*R1       -   M*R2; %P
        dF(5) = 0     - F*X(5) + 2*M*R1               ; %E
        dF(6) = 0     - F*X(6)                + 3*M*R2; %G
        
    end
end

