function [xOpt, cOpt, costOpt, consOpt] = CSTRopt(stru, x0)
% -------------------------------------------------------------------------
% finds the optimum point (minimum cost) subject to cons
%
% stru      - struct    - structure of parameters
% x0        - [2x1]     - u values for first guess
%
% xOpt      - [2x1]     - u values for optimum solution
%
% -------------------------------------------------------------------------

% check inputs
msg = 'Expected "%s" to be type ''%s'' but got type ''%s''';

if ~isa(x0, 'double')
    error(msg, 'x0', 'double', class(x0))
end

% run fmincon
options = optimoptions(@fmincon,'Display','Off');
uMin = [0, 0];
uMax = [50, 50];



[xOpt] = fmincon(@(x)(stru.cost(x, CSTRste(x, stru))), x0, [], [], [], [], uMin, uMax, @(x)cons(x), options);

if nargout > 2
    cOpt = CSTRste(xOpt, stru);
    costOpt = stru.cost(xOpt,cOpt);
    consOpt = stru.cons([],cOpt);
end

    function [c,ceq] = cons(x)
        c = stru.cons([],CSTRste(x,stru));
        ceq = [];
    end
end