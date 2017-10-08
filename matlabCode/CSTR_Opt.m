% finds the optimum parameters for a CSTR model for a given set of
% constraints

function [] = CSTR_Opt(reactOrder, kVal, c0, volFlow, V, modelCost, varargin)

% check inputs
if nargin < 7
    error('not enough inputs')
else
    constraints = varargin;
end

if ~isa(modelCost,'function_handle')
    error('modelCost is not the expected type, function_handle')
end


end



