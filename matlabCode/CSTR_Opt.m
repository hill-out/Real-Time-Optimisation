% finds the optimum parameters for a CSTR model for a given set of
% constraints

function [x_opt] = CSTR_Opt(cIn, x0, modelCost, varargin)

% check inputs
if nargin < 4
    error('not enough inputs')
else
    constraints = varargin{:};
end

if ~isa(modelCost,'function_handle')
    error('modelCost is not the expected type, function_handle')
end

[x_opt,FVAL,EXITFLAG,OUTPUT,LAMBDA,GRAD,HESSIAN] = fmincon(@calcCost, x0(1:2), [], [], [], [], [0, 0], [50, 50], @conFunc);

u = x0;
u(size(x_opt)) = x_opt;
cost = modelCost(u);
cons = cellFuncHandle2vec(constraints,{u},ones(size(constraints))');

    function [cost] = calcCost(x)
        % calculates the cost of a specific set of inputs at steady state
        u = x0;
        u(size(x)) = x;
        cost = modelCost(u);
    end

    function [c,ceq] = conFunc(x)
        % calculates the values of the constraint functions of a specific 
        % set of inputs at steady state
        u = x0;
        u(size(x)) = x;
        c = cellFuncHandle2vec(constraints,{u},ones(size(constraints))');
        ceq = [];
    end


end