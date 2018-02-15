% get convex parameters of model

FArange = linspace(2,10,11);
FBrange = linspace(4,15,11);
Trange = linspace(78,92,11);
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
phi_opt = fminsearch(@(u)(phiFun(u,CSTRmodel(u))),[3,6,84]);

% convPara
convPhi = convCalc(y',phi);



%plot
plotConv(y,phi,convPhi)