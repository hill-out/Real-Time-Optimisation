function [dudr] = truedudr(r0,Xp,cont,dr)


Xp0 = closedPlant(r0,Xp,cont);
u0 = cont(r0, Xp0);
dudr = zeros(2,3);

for i = 1:2
    r = r0 + dr(i,:);
    Xp = closedPlant(r,Xp0,cont);
    u = cont(r, Xp);
    
    dudr(i,:) = (u-u0)/dr(i,i);
    
end

end