function [c, cost, cons] = CSTRste(u, stru)
% -------------------------------------------------------------------------
% Calculates the steady state position of the CSTR
%
% u     - [1x2]     - inlet flowrates
% stru  - struct    - structure of parameters
%
% c     - [1x4]     - outlet concentrations
% cost  - double    - cost
% cons  - [1x2]     - constraints

% solves the CSTR at steady state
u = [u(1), u(2), 0, 0]';
v0 = sum(u);
options = optimset('Display','off');
c = fsolve(@massBal,stru.c_in,options);

if nargout > 1
    cost = stru.cost(u, c);
    cons = stru.cons(u, c);
    if nargout > 3
        convCost = stru.convCost(u, c);
        convCost = stru.convCons(u, c);
    end
end

    function [F] = massBal(c)
        %evaluates the species balance
        r = reactantRate(c, stru);
        F = u.*stru.c_in/stru.V+r-v0*c/stru.V;
    end

end