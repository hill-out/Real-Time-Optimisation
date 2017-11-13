function [dt, newC] = CSTRdyn(c_in, q_in, k, V, s, t, c0)

[dt, newC] = ode45(@(t,c)(CSTRode(c, c_in, q_in, k, V, s)),[0 t], c0);

end