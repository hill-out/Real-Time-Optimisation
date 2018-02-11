function [C, Ceq] = constrregconv(alp, param)

Q = eig([alp(3) alp(4);alp(4) alp(5)]);


C(1) =  param.eigmin1 - Q(1);
C(2) =  param.eigmin2 - Q(2);
Ceq = [];
    
end
