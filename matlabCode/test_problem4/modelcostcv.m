function [J] = modelcostcv(u,param)
alp = param.alp;
cvJinit = param.cvJinit;
refU = param.refU;
refUinit = param.refUinit;
epsJ = param.epscost;
lambda = param.lambdacost;
J = cvJinit + epsJ + [alp(1), alp(2)]*[u(1)-refUinit(1);u(2)-refUinit(2)] + lambda'*[u(1)-refU(1);u(2)-refU(2)] + 0.5*[u(1)-refUinit(1);u(2)-refUinit(2)]'*[alp(3) alp(4);alp(4) alp(5)]*[u(1)-refUinit(1);u(2)-refUinit(2)];
end
