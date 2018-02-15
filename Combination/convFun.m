function [funVal] = convFun(convPara, u)
% calculates the output of the convex approximation from a given set of
% convex paramters at a given set of inputs
% -------------------------------------------------------------------------
% convPara      Paramters of the approx [(0th [1x1]), (1st [1xnu]), (2nd 1xsum(1:nu)])]
% u             Relative inputs [nu x m]
%
% funVal        Outputs [1x1]
% -------------------------------------------------------------------------

nu = size(u,1);
m = size(u,2);

m0 = convPara(1);
m1 = convPara(2:nu+1);
locat = symMat(nu);
num = convPara(nu+2:end);
m2 = num(locat);
if any(eig(m2)<0)
    %m2 = m2*0;
end

quad = 0.5*u'*m2*u;

funVal = m0 + (m1*u)' + quad(logical(eye(m)));

end