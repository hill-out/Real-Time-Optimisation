function [phi, dphi] = phiCR(r)
% approximates phi from r
% -----------------------
% r         Setpoints
%
% phi       Cost
% -----------------------

convPara = [-328.976366706007,...
    -3364.30954123673,...
    7.42139674704701,...
    3691.19417491436,...
    -129.457031444715,...
    4.54029839567477];

rShift = bsxfun(@minus, r, [0.0899999997322743,12.7813268918873]');
phi = convCalc(convPara, rShift);
dphi = convdCalc(convPara, rShift);


end