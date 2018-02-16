% main
% run the RTO

% load conv para
load('convPara.mat')
y0 = [0.09, 12.7813];

% run plant
[a,b] = fmincon(@(y)(convFun(convPhi, (y-y0)')),[0.09,15],[],[],[],[],[-1 -10], [1, 50], ...
    @(y)(deal([convFun(convG1, (y-y0)'), convFun(convG2, (y-y0)')],[])));

B = [y0(2)/2.4, y0(2), (120-1000*(y0(1)-B(end,4))), 0.09, 0.51, 0.1, 0.1, 0.1, 0.1];
[A,B] = ode45(@(t,y)(ControlCSTRode(t,y,-1000)),[0 10],B);








