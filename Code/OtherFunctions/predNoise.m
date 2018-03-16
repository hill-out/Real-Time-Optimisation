function n = predNoise(t,tau,n0)
% creates rand numbers base on s, and gives the number based on t
% ---------------------------------------------------------------
% t         ode time
% tau       time steps
% n0        original noise
%
% n         noise
% ---------------------------------------------------------------

a = spline(linspace(0,tau,size(n0,1)),n0',linspace(0,tau,size(n0,1)*10));
b = diff(a');

n = spline(linspace(0,tau,size(b,1)),b',t);

end