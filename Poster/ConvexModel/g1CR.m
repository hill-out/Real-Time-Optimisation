function [g1, dg1] = g1CR(r)
% approximates phi from r
% -----------------------
% r         Setpoints
%
% phi       Cost
% -----------------------

convPara = [-2.67725716596168e-10,...
    0.999711980974046,...
    -2.95591996994342e-06,...
    0.000738802863488618,...
    7.68657811142979e-06,...
    3.76431135942724e-07];

rShift = bsxfun(@minus, r, [0.0899999997322743,12.7813268918873]');
g1 = convCalc(convPara, rShift);
dg1 = convdCalc(convPara, rShift);

end