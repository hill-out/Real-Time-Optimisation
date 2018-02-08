% Main
% Runs the RTO of a closed loop plant with an open loop model

% 1. load convex variables
[var_phi, var_g1, var_g2] = convApprox;
%load('convVar.mat')

var_g = [var_g1;var_g2];

% 2. run RTO

yOpt = fmincon(@(x)(convFun(var_phi(1),var_phi(2:3),var_phi(4:7),x,var_phi(8:end))),...
    var_phi(8:end),[],[],[],[],[0,-2],[1,50],...
    @(x)(convCon(var_g(:,1),var_g(:,2:3),var_g(:,4:7),x,var_g(1,8:end))));

r = yOpt;
u = plantController(r);
X_p = CSTRplant(u);
y_p = y_from_Xu(X_p,u);
phi_p = phi_from_Xu(X_p,u);
g_p = g_from_Xu(X_p,u);

[~,~,dphi_p,dg_p] = plantGrad(r);

g = convCon(var_g(:,1),var_g(:,2:3),var_g(:,4:7),r,var_g(1,8:end))';
dphi = convGrad(var_phi(2:3),var_phi(4:7),r,var_phi(8:end))';
dg(1,:) = convGrad(var_g1(2:3),var_g1(4:7),r,var_g1(8:end));
dg(2,:) = convGrad(var_g2(2:3),var_g2(4:7),r,var_g2(8:end));

K = 0.5; %gain

eg = K*(g_p - g);
lphi =  K*(dphi_p - dphi);
lg = K*(dg_p - dg);

modFun = @(x)(convFun(var_phi(1),var_phi(2:3),var_phi(4:7),x,var_phi(8:end)) + ...
    lphi*(x-r)');
modCon = @(x)(convCon(var_g(:,1),var_g(:,2:3),var_g(:,4:7),x,var_g(1,8:end)) + ...
    eg' + lg*(x-r)');


unsolved = 1;
i = 2;

while unsolved
    
    yOpt(i,:) = fmincon(@(x)(modFun(x)),...
        var_phi(8:end),[],[],[],[],[0,-5],[1,50],...
        @(x)deal((modCon(x)),[]));
    
    r = yOpt(i,:);
    u = plantController(r);
    X_p = CSTRplant(u);
    y_p = y_from_Xu(X_p,u);
    phi_p = phi_from_Xu(X_p,u);
    g_p = g_from_Xu(X_p,u);
    
    [~,~,dphi_p,dg_p] = plantGrad(r);
    
    g = convCon(var_g(:,1),var_g(:,2:3),var_g(:,4:7),r,var_g(1,8:end))';
    dphi = convGrad(var_phi(2:3),var_phi(4:7),r,var_phi(8:end))';
    dg(1,:) = convGrad(var_g1(2:3),var_g1(4:7),r,var_g1(8:end));
    dg(2,:) = convGrad(var_g2(2:3),var_g2(4:7),r,var_g2(8:end));
    
    eg = (1-K)*eg + K*(g_p - g);
    lphi =  (1-K)*lphi + K*(dphi_p - dphi);
    lg = (1-K)*lg + K*(dg_p - dg);
    
    modFun = @(x)(convFun(var_phi(1),var_phi(2:3),var_phi(4:7),x,var_phi(8:end)) + ...
        lphi*(x-r)');
    modCon = @(x)(convCon(var_g(:,1),var_g(:,2:3),var_g(:,4:7),x,var_g(1,8:end)) + ...
        eg' + lg*(x-r)');
    
    if all(abs((yOpt(i,:) - yOpt(i-1,:))./(yOpt(i,:))) < 1e-3)
        unsolved = 0;
    else
        i = i + 1;
    end
end