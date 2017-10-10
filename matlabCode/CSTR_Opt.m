% finds the optimum parameters for a CSTR model for a given set of
% constraints

function [x_opt] = CSTR_Opt(reactOrder, kVal, c0, volFlow, V, modelCost, varargin)

% check inputs
if nargin < 7
    error('not enough inputs')
else
    constraints = varargin{:};
end

if ~isa(modelCost,'function_handle')
    error('modelCost is not the expected type, function_handle')
end

options =  optimoptions('fmincon','Algorithm','interior-point','SpecifyObjectiveGradient',true,... 
    'SpecifyConstraintGradient',true,'SubproblemAlgorithm','cg','HessianMultiplyFcn',@HessMultFcn);

[x_opt,FVAL,EXITFLAG,OUTPUT,LAMBDA,GRAD,HESSIAN] = fmincon(@calcCost, volFlow(1:2), [], [], [], [], [0, 0], [50, 50], @conFunc);

cSol = CSTR(reactOrder, kVal, c0, [x_opt(1),x_opt(2),0,0], V);
cost = modelCost(c0, cSol, [x_opt(1),x_opt(2),0,0]);
cons = cellFuncHandle2vec(constraints,{c0,cSol,[x_opt(1),x_opt(2),0,0]},[1:3;1:3]);

    function [cost] = calcCost(x)
        % calculates the cost of a specific set of inputs at steady state
        cSol = CSTR(reactOrder, kVal, c0, [x(1),x(2),0,0], V);
        cost = modelCost(c0, cSol, x);
    end

    function [c,ceq] = conFunc(x)
        % calculates the values of the constraint functions of a specific 
        % set of inputs at steady state
        cSol = CSTR(reactOrder, kVal, c0, [x(1),x(2),0,0], V);
        c = cellFuncHandle2vec(constraints,{c0,cSol,[x(1),x(2),0,0]},[1:3;1:3]);
        ceq = [];
    end

    function [W] = 

end