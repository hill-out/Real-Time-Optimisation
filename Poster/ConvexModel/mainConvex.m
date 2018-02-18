% mainConvex

% opt
[u_opt, phi_opt] = fmincon(@(x)phiFun(x,CSTRmodel(x)),[4, 10, 80], ...
    [],[],[],[],[0,0,70],[20,50,120],@(x)(conFun(x,CSTRmodel(x))));
X_opt = CSTRmodel(u_opt);
g_opt = conFun(u_opt,X_opt);

% set up u
uRange = {linspace(u_opt(1)*0.2,u_opt(1)*2,11),...
    linspace(u_opt(2)*0.2,u_opt(2)*2,11),...
    linspace(u_opt(3)-6,u_opt(3)+6, 5)};
[u1,u2,u3] = ndgrid(uRange{:});
u = [u1(:), u2(:), u3(:)]';

% run 
[phi, X] = funRun(@phiFun, u);
g1 = (X(1,:)-0.09);
g2 = (X(6,:)-0.6);

gUpper = (X(1,:)<0.15);
gLower = (X(1,:)>0.03);

% gCut = (gUpper+gLower) == 2;
% 
% phi = phi(gCut);
% X = X(:,gCut);
% u = u(:,gCut);
% g1 = g1(gCut);
% g2 = g2(gCut);

% conv para U
convPhiU = convexParaU(phi, u, phi_opt, u_opt');
convG1U = convexParaU(g1, u, g_opt(1), u_opt');
convG2U = convexParaU(g2, u, g_opt(2), u_opt');

uShift = bsxfun(@minus, u, u_opt');

[uc_optU, phic_optU] = fmincon(@(x)convCalc(convPhiU, (x-u_opt)'),[4, 10, 80],...
    [],[],[],[],[0,0,70],[20,50,120],...
    @(x)deal([convCalc(convG1U, (x-u_opt)'),convCalc(convG2U, (x-u_opt)')],[]));
Xc_optU = CSTRmodel(uc_optU);
gc_optU = conFun(uc_optU,Xc_optU);

% conv para R
r = [X(1,:); u(2,:)];
r_opt = [X_opt(1), u_opt(2)];

convPhiR = convexParaR(phi, r, phi_opt, r_opt');
convG1R = convexParaR(g1, r, g_opt(1), r_opt');
convG2R = convexParaR(g2, r, g_opt(2), r_opt');

rShift = bsxfun(@minus, r, r_opt');

[rc_optR, phic_optR] = fmincon(@(x)convCalc(convPhiR, (x-r_opt)'),r_opt,...
    [],[],[],[],[0,0],[1,50],...
    @(x)deal([convCalc(convG1R, (x-r_opt)'),convCalc(convG2U, (x-r_opt)')],[]));
Xc_optR = CSTRmodel(uc_optU);
gc_optR = [convCalc(convG1R, (rc_optR-r_opt)'),convCalc(convG2U, (rc_optR-r_opt)')];

% plot
a = (g1<0);
b1 = phi(a);
b2 = phi(~a);
c = u(1,:);
c1 = c(a);
c2 = c(~a);
d = u(2,:);
d1 = d(a);
d2 = d(~a);
scatter3(c1,d1,b1)
hold on
scatter3(c2,d2,b2)






