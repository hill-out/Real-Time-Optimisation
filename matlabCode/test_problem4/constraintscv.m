function [G, Geq] = constraintscv(u,param)
alp1 = param.alp1;
alp2 = param.alp2;
cvJinit = param.cvJinit;
refU = param.refU;
refUinit = param.refUinit;
G1init = param.G1init;
G2init = param.G2init;
epsG1 = param.epsconst1;
epsG2 = param.epsconst2;
lambdaG1 = param.lambdaconst1;
lambdaG2 = param.lambdaconst2;
G1 =G1init + epsG1 + [alp1(1), alp1(2)]*[u(1)-refUinit(1);u(2)-refUinit(2)] + lambdaG1'*[u(1)-refU(1);u(2)-refU(2)] + 0.5*[u(1)-refUinit(1);u(2)-refUinit(2)]'*[alp1(3) alp1(4);alp1(4) alp1(5)]*[u(1)-refUinit(1);u(2)-refUinit(2)];
G2 =G2init + epsG2 + [alp2(1), alp2(2)]*[u(1)-refUinit(1);u(2)-refUinit(2)] + lambdaG2'*[u(1)-refU(1);u(2)-refU(2)] + 0.5*[u(1)-refUinit(1);u(2)-refUinit(2)]'*[alp2(3) alp2(4);alp2(4) alp2(5)]*[u(1)-refUinit(1);u(2)-refUinit(2)];
G = [G1; G2];
Geq = [];
end
