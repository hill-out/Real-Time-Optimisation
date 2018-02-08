function [val] = convGrad(b,c,in,inOpt)
% finds the values of a convex function from the variables
% --------------------------------------------------------
% a         Zeroth order term (1x1)
% b         First order term (nx1)
% c         Second order term (nxn)
% in        Inputs to the convex approximation (nx1)
% 
% out       Outputs to the convex approximation (1x1)
% --------------------------------------------------------

n = numel(in);

b = reshape(b,n,1);         % make a row
c = reshape(c,n,n);         % make a square
in = reshape(in,n,1);       % make a column
inOpt = reshape(inOpt,n,1); % make a column

val = b + 0.5*c*(in-inOpt) + 0.5*c'*(in-inOpt);

end