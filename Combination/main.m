% main
% run the RTO

% load conv para
load('convPara.mat')
y0 = [0.09, 12.7813];
tau = 900;
gain = 0.3;
K = -100;
T0 = 120;
dr = [0.0001, 0; 0, 0.01];

% run plant
[a,b] = fmincon(@(y)(convFun(convPhi, (y-y0)')),[0.09,15],[],[],[],[],[-1 -10], [1, 50], ...
    @(y)(deal([convFun(convG1, (y-y0)'), convFun(convG2, (y-y0)')],[])));

model.yOpt = a;
model.phi = b;
model.g = [convFun(convG1, (a-y0)'), convFun(convG2, (a-y0)')];

plant.base.r = y0;
B = [(plant.base.r(2)+2)/2.4, plant.base.r(2)+2, (T0+K*(plant.base.r(1)-0.09)), 0.09, 0.51, 0.1, 0.1, 0.1, 0.1];
[Ab,B] = ode45(@(t,y)(ControlCSTRode(t,y,K)),[0 200000],B);

% get plant values
plant.base.u = B(end,1:3);
plant.base.X = B(end,4:end);

plant.base.phi = phiFun(plant.base.u,plant.base.X);
plant.base.g = conFun(plant.base.u,plant.base.X);

% run other units
for ri = 1:2
    r = plant.base.r + dr(ri,:);
    B = [(r(2)+2)/2.4, r(2)+2, (T0+K*(r(1)-plant.base.X(1))),plant.base.X];
    [Ao,C] = ode45(@(t,y)(ControlCSTRode(t,y,K)),[0 tau],B);
    plant.other(ri).u = C(end,1:3);
    plant.other(ri).X = C(end,4:end);
    
    plant.other(ri).phi = phiFun(plant.other(ri).u,plant.other(ri).X);
    plant.other(ri).g = conFun(plant.other(ri).u,plant.other(ri).X);
end

% get plant gradients
plant.dphi = [(plant.other(1).phi - plant.base.phi)/dr(1,1),...
    (plant.other(2).phi - plant.base.phi)/dr(2,2)];
plant.dg = [(plant.other(1).g(1) - plant.base.g(1))/dr(1,1),...
    (plant.other(2).g(1) - plant.base.g(1))/dr(2,2);
    (plant.other(1).g(2) - plant.base.g(2))/dr(1,1),...
    (plant.other(2).g(2) - plant.base.g(2))/dr(2,2)];

% get model gradients
model.dphi = [(convFun(convPhi, (a+dr(1,:)-y0)')-b)/dr(1,1),...
    (convFun(convPhi, (a+dr(2,:)-y0)')-b)/dr(2,2)];
model.dg = [(convFun(convG1, (a+dr(1,:)-y0)')-model.g(1))/dr(1,1),...
    (convFun(convG1, (a+dr(2,:)-y0)')-model.g(1))/dr(2,2);
    (convFun(convG2, (a+dr(1,:)-y0)')-model.g(2))/dr(1,1),...
    (convFun(convG2, (a+dr(2,:)-y0)')-model.g(2))/dr(2,2)];

% modifiers

model.ephi = gain*(plant.base.phi - model.phi);
model.eg = gain*(plant.base.g - model.g);
model.lphi = gain*(plant.dphi - model.dphi);
model.lg = gain*(plant.dg(:) - model.dg(:))';

model.fun = @(y)(convFun(convPhi, (y-y0)')+model.ephi+(y-y0)*model.lphi');
model.con = @(y)deal([(convFun(convG1, (y-y0)')+model.eg(1)+(y-y0)*model.lg([1,3])'),...
    (convFun(convG2, (y-y0)')+model.eg(2)+(y-y0)*model.lg([2,4])')],[]);

unsolved = 1;
i=2;

while unsolved
    [a,b] = fmincon(model.fun,[0.09,15],[],[],[],[],[0 0], [1, 50], ...
        model.con);
    
    model.yOpt(i,:) = a;
    
    plant.base.r(i,:) = a;
    B = [(plant.base.r(i,2)+2)/2.4, plant.base.r(i,2)+2, (T0-K*(plant.base.r(i,1)-plant.base.X(1))), plant.base.X];
    [~,C] = ode45(@(t,y)(ControlCSTRode(t,y,-K)),[0 tau],B);
    
    % get plant values
    plant.base.u = C(end,1:3);
    plant.base.X = C(end,4:end);
    
    plant.base.phi = phiFun(plant.base.u,plant.base.X);
    plant.base.g = conFun(plant.base.u,plant.base.X);
    
    % run other units
    for ri = 1:2
        r = plant.base.r(i,:) + dr(ri,:);
        B = [(r(2)+2)/2.4, r(2)+2, (T0-K*(r(1)-plant.base.X(1))),plant.base.X];
        [~,C] = ode45(@(t,y)(ControlCSTRode(t,y,-K)),[0 tau],B);
        plant.other(ri).u = C(end,1:3);
        plant.other(ri).X = C(end,4:end);
        
        plant.other(ri).phi = phiFun(plant.other(ri).u,plant.other(ri).X);
        plant.other(ri).g = conFun(plant.other(ri).u,plant.other(ri).X);
    end
    
    % get plant gradients
    plant.dphi = [(plant.other(1).phi - plant.base.phi)/dr(1,1),...
        (plant.other(2).phi - plant.base.phi)/dr(2,2)];
    plant.dg = [(plant.other(1).g(1) - plant.base.g(1))/dr(1,1),...
        (plant.other(2).g(1) - plant.base.g(1))/dr(2,2);
        (plant.other(1).g(2) - plant.base.g(2))/dr(1,1),...
        (plant.other(2).g(2) - plant.base.g(2))/dr(2,2)];
    
    % get model gradients/values
    model.phi = convFun(convPhi, (a-y0)');
    model.g = [convFun(convG1, (a-y0)'),convFun(convG2, (a-y0)')];
    model.dphi = [(convFun(convPhi, (a+dr(1,:)-y0)')-convFun(convPhi, (a-y0)'))/dr(1,1),...
        (convFun(convPhi, (a+dr(2,:)-y0)')-convFun(convPhi, (a-y0)'))/dr(2,2)];
    model.dg = [(convFun(convG1, (a+dr(1,:)-y0)')-convFun(convG1, (a-y0)'))/dr(1,1),...
        (convFun(convG1, (a+dr(2,:)-y0)')-convFun(convG1, (a-y0)'))/dr(2,2);...
        (convFun(convG2, (a+dr(1,:)-y0)')-convFun(convG2, (a-y0)'))/dr(1,1),...
        (convFun(convG2, (a+dr(2,:)-y0)')-convFun(convG2, (a-y0)'))/dr(2,2)];
    
    % modifiers
    
    model.ephi = model.ephi*(1-gain) + gain*(plant.base.phi - model.phi);
    model.eg = model.eg*(1-gain) + gain*(plant.base.g - model.g);
    model.lphi = model.lphi*(1-gain) + gain*(plant.dphi - model.dphi);
    model.lg = model.lg*(1-gain) + gain*(plant.dg(:) - model.dg(:))';
    
    model.fun = @(y)(convFun(convPhi, (y-y0)')+model.ephi+(y-a)*model.lphi');
    model.con = @(y)deal([(convFun(convG1, (y-a)')+model.eg(1)+(y-a)*model.lg([1,3])'),...
        (convFun(convG2, (y-a)')+model.eg(2)+(y-a)*model.lg([2,4])')],[]);
    
    i=i+1
end


