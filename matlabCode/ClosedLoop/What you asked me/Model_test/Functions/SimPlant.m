function [up, gp, phi_up, phi_c, phi_p, x, F_P] = SimPlant(c, x0_P, F_P_demand, Parameters, fsolve_options)

x = fsolve(@(x_P) ClosedLoopSystem(0, x_P, c, Parameters), x0_P, fsolve_options);
[up, gp, phi_up, phi_c, phi_p, F_P] = cx2ugphi(c, x, F_P_demand, Parameters);

end