function [funVal] = convFun(convPara, u)
% calculates the output of the convex approximation from a given set of
% convex paramters at a given set of inputs
% -------------------------------------------------------------------------
% convPara      Paramters of the approx [(0th [1x1]), (1st [1xn]), (2nd 1xn^2])]
% u             Relative inputs [nx1]
%
% funVal        Outputs [1x1]
% -------------------------------------------------------------------------

n = numel(u);
u = reshape(u, n, 1);
m0 = convPara(1);
m1 = convPara(2:n+1);
m2 = reshape(convPara(n+1:end),n,n);

funVal = m0 + m1'*u + u'*m2*u;

end