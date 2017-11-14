function dc = CSTRode(c, c_in, q_in, k, V, s)
% the ode of a CSTR
%
% c - Current concentrations [nx1] (mol/L)
% c_in - Inlet concentrations [nx1] (mol/L)
% q_in - Inlet flowrates [nx1] (L/min)
% k - Reaction constants [1xm] L/(mol min)
% V - Volume (L)
% s - Reaction order [nxm] 
%
% dc - gradient of c

%% Gradient of c
% 
% $$\frac{dc}{dt} = \frac{q_{in}}{V}\times(c_{in}-c) + r$$
% 

r = reactantRate(c, k, s);

dc = (q_in/V).*c_in-(sum(q_in)/V).*(c)+r;

end