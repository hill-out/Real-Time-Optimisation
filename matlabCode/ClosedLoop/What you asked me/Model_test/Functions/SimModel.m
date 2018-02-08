function [g, phi_u, x_M, F_P] = SimModel(u, x0_M, F_P_demand, Parameters, fsolve_options)

x_M = fsolve(@(x_M) OpenLoopSystem(0, x_M, u, Parameters), x0_M, fsolve_options);
[g, phi_u, F_P] = ux2gphi(u, x_M, F_P_demand, Parameters);


end