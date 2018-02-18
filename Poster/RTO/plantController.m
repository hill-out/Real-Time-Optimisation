function [u] = plantController(r)
% finds the input values for a given value of r
% ---------------------------------------------
% r         controller setpoint
% 
% u         inputs
% ---------------------------------------------

xGuess = zeros(1,6);
u = zeros(3,1);
u(2) = r(2)+2;
u(1) = u(2)/2.4;

optionU = optimoptions('fsolve','Display','off');

u(3) = fsolve(@(x)(tempFinder(x)),78,optionU);

    function dxa = tempFinder(x)
        a = u;
        a(3) = x;
        
        X = CSTRplant(a, xGuess);
        xGuess = X;
        dxa = X(1) - 0.75*r(1);
    end

end