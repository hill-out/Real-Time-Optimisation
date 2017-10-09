% finds the optimum parameters for a CSTR model for a given set of
% constraints

function [] = CSTR_Opt(reactOrder, kVal, c0, volFlow, V, modelCost, varargin)

% check inputs
if nargin < 7
    error('not enough inputs')
else
    constraints = varargin{:};
end

if ~isa(modelCost,'function_handle')
    error('modelCost is not the expected type, function_handle')
end

x_opt = fmincon(@calcCost, volFlow, [], [], [], [], [0, 0, 0, 0], [50, 50, 0, 0], @conFunc);

cSol = CSTR(reactOrder, kVal, c0, x_opt, V);
cost = modelCost(c0, cSol, x_opt);

    function [cost] = calcCost(x)
        % calculates the cost of a specific set of inputs at steady state
        cSol = CSTR(reactOrder, kVal, c0, x, V);
        cost = modelCost(c0, cSol, x);
    end

    function [c,ceq] = conFunc(x)
        % calculates the values of the constraint functions of a specific 
        % set of inputs at steady state
        % cSol = CSTR(reactOrder, kVal, c0, x, V);
        c = cellFuncHandle2vec(constraints,{c0,cSol,x},[1:3;1:3]);
        ceq = [];
    end

end