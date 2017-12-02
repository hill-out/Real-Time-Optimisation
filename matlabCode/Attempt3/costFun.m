function [cost] = costFun(u, c, stru)
% -------------------------------------------------------------------------
% calculates the cost
% 
% u     - [1x2]     - inlet flowrates
% c     - [1x4]     - concentrations
% stru  - struct    - structure of parameters
%
% cost  - double    - cost of the current conditions
%
% -------------------------------------------------------------------------

w = 0.004;
J = c(3,:).^2.*(u(1)+u(2))^2./(u(1)*stru.c_in(1)) - w*(u(1)^2+u(2)^2);
cost = -J;

end