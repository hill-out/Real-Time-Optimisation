function [df] = convdCalc(convPara, u)
% Calcutes the value of the df from u
% -----------------------------------
% convPara  Parameters of the approx
% u         Inputs to approx
%
% df        Gradients of outputs
% -----------------------------------

n = size(u,1);

% get variables
m0 = convPara(1);

m1 = convPara(2:n+1);

q = convPara(n+2:end);
locat = symMat(n);
m2 = q(locat);

% calc approx
df = m1 + 2*u'*m2;

end
