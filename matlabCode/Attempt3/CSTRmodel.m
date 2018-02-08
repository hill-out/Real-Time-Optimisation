function [Xsol] = CSTRmodel(u)
% Symbolically solves the model 2-reaction CSTR
% ---------------------------------------------
% u         Inputs
% 
% Xsol      Mass Fraction [A,B,P,E,G]
% ---------------------------------------------

F_Ain = u(1);
F_Bin = u(2);
T     = u(3);

M = 2105; %kg
F = F_Ain + F_Bin; %kg/s
k1 = 1.4091e7*exp(-E(1)/(T)); %1/s
k2 = 2.3376e8*exp(-E(2)/(T)); %1/s

%% mass balances


Xsol = fsolve(@massBalance, zeros(1,5));

    function [dF] = massBalance(X)
        
        R1 = k1*X(1)*X(2)^2;
        R2 = k2*X(1)*X(2)*X(3);
        %         IN    - OUT    +GEN/-CON R1 +GEN/-CON R2
        dF(1) = F_Ain - F*X(1) -   M*R1       -   M*R2;
        dF(2) = F_Bin - F*X(2) - 2*M*R1       -   M*R2;
        dF(3) = 0     - F*X(3) +   M*R1       -   M*R2;
        dF(4) = 0     - F*X(4) + 2*M*R1               ;
        dF(5) = 0     - F*X(5)                + 3*M*R2;
        
    end
        
        
end

