function [u] = plantController3(r, x, Kp, T0)
% finds the input values for a given value of r
% ---------------------------------------------
% r         controller setpoint
% 
% u         inputs
% ---------------------------------------------

u = zeros(3,1);
u(2) = r(2)*2.4/3.4+2;
u(1) = u(2)/2.4;

% Temperature controller
u(3) = T0 + Kp*(r(1) - x(1));
end