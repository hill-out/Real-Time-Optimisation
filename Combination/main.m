% main
% run the RTO

% load conv para
load('convPara.mat')
y0 = [0.09, 12.7813];
tau = 10;
K = 0.5;
dr = [0.0001, 0; 0, 0.01];

% run plant
[a,b] = fmincon(@(y)(convFun(convPhi, (y-y0)')),[0.09,15],[],[],[],[],[-1 -10], [1, 50], ...
    @(y)(deal([convFun(convG1, (y-y0)'), convFun(convG2, (y-y0)')],[])));

model.yOpt = a;
model.phi = b;
model.g = [convFun(convG1, (a-y0)'), convFun(convG2, (a-y0)')];

plant.base.r = y0;
B = [(plant.base.r(2)+2)/2.4, plant.base.r(2)+2, (120-1000*(plant.base.r(1)-0.09)), 0.09, 0.51, 0.1, 0.1, 0.1, 0.1];
[~,B] = ode45(@(t,y)(ControlCSTRode(t,y,-1000)),[0 120],B);

% get plant values
plant.base.u = B(end,1:3);
plant.base.X = B(end,4:end);

plant.base.phi = phiFun(plant.base.u,plant.base.X);
plant.base.g = conFun(plant.base.u,plant.base.X);

% run other units
for ri = 1:2
    r = plant.base.r + dr(ri,:);
    B = [(r(2)+2)/2.4, r(2)+2, (120-1000*(r(1)-B(end,4))),...
        B(end,4:end)];
    [~,C] = ode45(@(t,y)(ControlCSTRode(t,y,-1000)),[0 tau],B);
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

model.ephi = K*(plant.base.phi - model.phi);
model.eg = K*(plant.base.g - model.g);
model.lphi = K*(plant.dphi - model.dphi);
model.lg = K*(plant.dg(:) - model.dg(:))';

model.fun = @(y)(convFun(convPhi, (y-y0)')+model.ephi+(y-y0)*model.lphi');
model.con = @(y)deal([(convFun(convG1, (y-y0)')+model.eg(1)+(y-y0)*model.lg([1,3])'),...
    (convFun(convG2, (y-y0)')+model.eg(2)+(y-y0)*model.lg([2,4])')],[]);

unsolved = 1;

while unsolved
    [a,b] = fmincon(model.fun,[0.09,15],[],[],[],[],[0 0], [1, 50], ...
        model.con);
    
    model.yOpt = a;
    model.phi = b;
    model.g = [convFun(convG1, (a-y0)'), convFun(convG2, (a-y0)')];
    
    plant.base.r = a;
    B = [(plant.base.r(2)+2)/2.4, plant.base.r(2)+2, (120-1000*(plant.base.r(1)-0.09)), 0.09, 0.51, 0.1, 0.1, 0.1, 0.1];
    [~,B] = ode45(@(t,y)(ControlCSTRode(t,y,-1000)),[0 tau],B);
    
    % get plant values
    plant.base.u = B(end,1:3);
    plant.base.X = B(end,4:end);
    
    plant.base.phi = phiFun(plant.base.u,plant.base.X);
    plant.base.g = conFun(plant.base.u,plant.base.X);
    
    % run other units
    for ri = 1:2
        r = plant.base.r + dr(ri,:);
        B = [(r(2)+2)/2.4, r(2)+2, (120-1000*(r(1)-B(end,4))),...
            B(end,4:end)];
        [~,C] = ode45(@(t,y)(ControlCSTRode(t,y,-1000)),[0 tau],B);
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
    model.dphi = [(convFun(convPhi, (a+dr(1,:)-y0)')-(convFun(convPhi, (a-y0)'))/dr(1,1)),...
        (convFun(convPhi, (a+dr(2,:)-y0)')-convFun(convPhi, (a-y0)'))/dr(2,2)];
    model.dg = [(convFun(convG1, (a+dr(1,:)-y0)')-convFun(convG1, (a-y0)'))/dr(1,1),...
        (convFun(convG1, (a+dr(2,:)-y0)')-convFun(convG1, (a-y0)'))/dr(2,2);...
        (convFun(convG2, (a+dr(1,:)-y0)')-convFun(convG2, (a-y0)'))/dr(1,1),...
        (convFun(convG2, (a+dr(2,:)-y0)')-convFun(convG2, (a-y0)'))/dr(2,2)];
    
    % modifiers
    
    model.ephi = K*(plant.base.phi - model.phi);
    model.eg = K*(plant.base.g - model.g);
    model.lphi = K*(plant.dphi - model.dphi);
    model.lg = K*(plant.dg(:) - model.dg(:))';
    
    model.fun = @(y)(convFun(convPhi, (y-y0)')+model.ephi+(y-y0)*model.lphi');
    model.con = @(y)deal([(convFun(convG1, (y-y0)')+model.eg(1)+(y-y0)*model.lg([1,3])'),...
        (convFun(convG2, (y-y0)')+model.eg(2)+(y-y0)*model.lg([2,4])')],[]);
end


