% get convex parameters of model

FArange = linspace(2,10,11);
FBrange = linspace(4,15,11);
Trange = linspace(78,92,3);
[FA, FB, T] = ndgrid(FArange, FBrange, Trange);

xGuess = zeros(1,6);
X = zeros(numel(FA),6);

for i = 1:numel(FA)
    u = [FA(i), FB(i), T(i)];
    X(i,:) = CSTRmodel(u, xGuess);
    xGuess = X(i,:);
end

% y = [XA, FB]
y = zeros(numel(FA),2);
y(:,2) = FB(:);
y(:,1) = X(:,1);

% phi = phiFun(u,x)
u = [FA(:),FB(:),T(:)];
phi = phiFun(u,X);

% constraints
g1 = X(:,1) - 0.09;
g2 = X(:,6) - 0.6;

% min
[u_opt,phi_opt] = fmincon(@(u)(phiFun(u,CSTRmodel(u))),[3,6,84],[],[],[],[],...
    [0 0 70], [20 48 150], @(u)(conFun(u,CSTRmodel(u))));
X_opt = CSTRmodel(u_opt);
y_opt(1) = X_opt(1);
y_opt(2) = u_opt(2);
g_opt = conFun(u_opt,X_opt);

% convPara
convPhi = convCalc(y',phi,y_opt',phi_opt);
convG1 = convCalc(y',g1,y_opt',g_opt(1));
convG2 = convCalc(y',g2,y_opt',g_opt(2));

%reassign
convPhi = [phi_opt,convPhi];
convG1 = [g_opt(1),convG1];
convG2 = [g_opt(2),convG2];