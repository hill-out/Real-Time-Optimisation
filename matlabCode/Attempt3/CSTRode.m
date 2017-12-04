function dc = CSTRode(u, c, stru)
% -------------------------------------------------------------------------
% the ode of a CSTR
%
% u         - [1x2]     - current flowrates in
% c         - [1x4]     - current concentrations
% stru      - struct    - structure of parameters
%
% dc        - [1x4]     - gradient of c
%
% -------------------------------------------------------------------------

%% Gradient of c

r = reactantRate(c, stru);
q = [u(1), u(2), 0, 0]'; % all flowrates
dc = (q/stru.V).*stru.c_in-(sum(q)/stru.V).*(c)+r;

end