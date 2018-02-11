function [J] = costregconv(alp,param)

u1 = param.ureg1;
u2 = param.ureg2;
phi = param.cvJ;
refU = param.refU;
cvJinit = param.cvJinit;
J = 0;
testG = param.testG;
if testG == 0
for i=1:length(u1),
     J = J + ((cvJinit + [alp(1), alp(2)]*[u1(i)-refU(1);u2(i)-refU(2)] + 0.5*[u1(i)-refU(1);u2(i)-refU(2)]'*[alp(3) alp(4);alp(4) alp(5)]*[u1(i)-refU(1);u2(i)-refU(2)])-phi(i))^2;
end
else
 for i=1:length(u1),
     J = J + ((cvJinit + [alp(1), alp(2)]*[u1(i)-refU(1);u2(i)-refU(2)])-phi(i))^2;
end   
end

