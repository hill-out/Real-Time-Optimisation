function [f] = convCalc(convPara, u)
% Calcutes the value of the f from u
% ----------------------------------
% convPara  Parameters of the approx
% u         Inputs to approx [nxm]
%
% f         Outputs of approx [nx1]
% ----------------------------------

n = size(u,1);

% get variables
m0 = convPara(1);

m1 = convPara(2:n+1);

q = convPara(n+2:end);
locat = symMat(n);
m2 = q(locat);

% calc approx
f = m0 + m1*u + sum(u'*m2.*u',2)';

end