function [y] = convModel(p, x, x0, y0)
% -------------------------------------------------------------------------
% the convex model outputs
%
% p     - [3x1]     - parameters of model
% x     - [1x2]     - value of x
% x0    - [1x2]     - optimum value of x
% y0    - double    - optimum value of y at x
%
% y     - double    - value of y at x
%
% -------------------------------------------------------------------------

dx = (x-x0)';
y = p(1:2)*dx + 0.5*p(3).*dx'*dx + y0;

end