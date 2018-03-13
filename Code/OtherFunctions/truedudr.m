function [dudr] = truedudr(r0,Xp,Kp,T0,dr)

r = r0;
u = plantController2(r,Xp,Kp,T0)';
[~,a] = ode15s(@(t,y)closedPlantODE(t,y,Kp), [0 3000],[u, Xp]);
u0 = a(end,1:3);
dudr = zeros(2,3);

for i = 2:-1:1
    Xp = a(end,4:end);
    r = r0 + dr(i,:);
    u = plantController2(r,Xp,Kp,T0)';
    [~,a] = ode15s(@(t,y)closedPlantODE(t,y,Kp), [0 3000],[u, Xp]);
    u1 = a(end,1:3);
    dudr(i,:) = (u1-u0)/dr(i,i);
end

end