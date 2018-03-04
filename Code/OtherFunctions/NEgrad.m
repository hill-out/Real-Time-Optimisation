function [dPlant] = NEgrad(u,dOpt)
% Approximates the values of the function using NE
% ------------------------------------------------
% u         Value to approximate at
% dOpt      Structure of difference to optimums 
%
% dPlant    Gradient at plant
% ------------------------------------------------

[dfun] = openModelGrad(u);

dPlant.dphidu = dfun.ddphidudu*dOpt.du' + ...
    dfun.ddphidudT*(dOpt.du*dfun.dCdu*pinv(dfun.dCdT) + dOpt.dC*pinv(dfun.dCdT))';

dPlant.dg1du = dfun.ddg1dudu*dOpt.du' + ...
    dfun.ddg1dudT*(dOpt.du*dfun.dCdu*pinv(dfun.dCdT) + dOpt.dC*pinv(dfun.dCdT))';

dPlant.dg2du = dfun.ddg2dudu*dOpt.du' + ...
    dfun.ddg2dudT*(dOpt.du*dfun.dCdu*pinv(dfun.dCdT) + dOpt.dC*pinv(dfun.dCdT))';
end