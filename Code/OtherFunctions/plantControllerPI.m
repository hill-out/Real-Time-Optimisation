function [u] = plantControllerPI(r, Kp, Ki, T0, delta)
% finds the input values for a given value of r
% ---------------------------------------------
% r         controller setpoint
% 
% u         inputs
% ---------------------------------------------

u = zeros(3,1);
u(2) = r(2)+2;
u(1) = u(2)/2.4;

% Temperature controller
u(3) = T0 + Kp*(delta(end,2)) + Ki*trapz(delta(:,1),delta(:,2));
end