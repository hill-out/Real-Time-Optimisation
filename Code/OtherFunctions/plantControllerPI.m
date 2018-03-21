function [u] = plantControllerPI(r, x, Kp, Ki, T0)
% finds the input values for a given value of r
% ---------------------------------------------
% r         controller setpoint
% 
% u         inputs
% ---------------------------------------------
global error

u = zeros(3,1);
u(2) = r(2)+2;
u(1) = u(2)/2.4;

% Temperature controller
u(3) = T0 + Kp*(r(1)- x(1)) + Ki*error;
end