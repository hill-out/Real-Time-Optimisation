syms XA XB XP XE XG FA FB F T E1 E2

F = FA + FB;

k_0 = [2.189*1e8, 4.310*1e13]; %1/s
M = 2105;

k1 = k_0(1)*exp(-E1/(T+273.15)); %1/s
k2 = k_0(2)*exp(-E2/(T+273.15)); %1/s

R1 = k1*XA*XB*XB;
R2 = k2*XA*XB*XP;

eq1 = FA/M - XA*F/M + (-R1 - R2);
eq2 = FB/M - XB*F/M + (-2*R1 - R2);
eq3 = - XP*F/M + (R1 - R2);
eq4 = - XE*F/M + 2*(R1);
eq5 = - XG*F/M + 3*(R2);

sol = solve(eq1,eq2,eq3,eq4,eq5, XA, XB, XP, XE, XG);