function [g2, dg2] = g2CR(r)
% approximates phi from r
% -----------------------
% r         Setpoints
%
% phi       Cost
% -----------------------

convPara = [-0.504946517445388,...
    0.00937816996753801,...
    -0.00570648001659968,...
    0.0258551806615696,...
    0.00151620022531883,...
    8.90421709455261e-05];

rShift = bsxfun(@minus, r, [0.0899999997322743,12.7813268918873]');
g2 = convCalc(convPara, rShift);
dg2 = convCalc(convPara, rShift);

end