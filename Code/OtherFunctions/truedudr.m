function [dudr] = truedudr(r0,Xp,Kp,T0,dr)

u = plantController3(r0,Xp,Kp,T0)';
[~,a] = ode15s(@(t,y)closedPlantODE(t,y,Kp), [0 10000],[u, Xp]);
u0 = a(end,1:3);

for i = 1:2
    r = r0 + dr(i,:);
    u = plantController3(r,Xp,Kp,T0)';
    [~,a] = ode15s(@(t,y)closedPlantODE(t,y,Kp), [0 10000],[u, Xp]);
    u1 = a(end,1:3);
    dudr(i,:) = (u1-u0)/dr(i,i);
end

end