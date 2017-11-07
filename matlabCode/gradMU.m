function [base, other, gradJ, gradG] = gradMU(fun, c0B, c0O, u0, n, tau, delta, base, other, gradJ, gradG)
% calculates the gradient of the plant using the multiple units method
%
% n - number of time steps
% tau - time to run for
%
%

% Base Case
if nargin < 8 || isempty(base)
    supB = 0; %not supplied base 
    base = struct('C', zeros(n+1, size(c0B,2)),...
        'J', zeros(n+1, 1),...
        'G', zeros(n+1, 2));
else 
    supB = 1; %supplied base
    lenB = size(base.C,1); %current length
    c0B = base.C(end,:);
end

if nargin < 9 || isempty(other)
    supO = 0; %not supplied base 
    other = struct('C', zeros(n+1, size(c0O,2), numel(u0)),...
        'J', zeros(n+1, 1, numel(u0)),...
        'G', zeros(n+1, 2, numel(u0)));
    
else 
    supO = 1; %supplied base
    lenO = size(other.C,1); %current length
    c0O = permute(other.C(end,:,:),[3,2,1]);
end


if nargin < 10 || isempty(gradJ)
    supJ = 0;
else
    supJ = 1;
    lenJ = size(gradJ,1);
end

if nargin < 11 || isempty(gradG)
    supG = 0;
else
    supG = 1;
    lenG = size(gradG,1);
end

[bC, bJ, bG] = fun(c0B,u0,n,tau);

if supB;
    base.C(lenB:lenB+n,:) = bC;
    base.J(lenB:lenB+n,:) = bJ;
    base.G(lenB:lenB+n,:) = bG;
else
    base.C = bC;
    base.J = bJ;
    base.G = bG;
end    

% Other Units

for i = 1:numel(u0)
    uNew = u0;
    uNew(i) = uNew(i) + delta;
    
    [oC, oJ, oG] = fun(c0O(i,:),uNew,n,tau);
    
    if supO
        other.C(lenO:lenO+n,:,i) = oC;
        other.J(lenO:lenO+n,:,i) = oJ;
        other.G(lenO:lenO+n,:,i) = oG;
    else
        other.C(:,:,i) = oC;
        other.J(:,:,i) = oJ;
        other.G(:,:,i) = oG;
    end
end

% Grad
if supJ
    gradJ(lenJ+1,:,:) = (other.J(end,:,:) - repmat(base.J(end,:),1,1,numel(u0)))/(delta);
else
    gradJ = (other.J(end,:,:) - repmat(base.J(end,:),1,1,numel(u0)))/(delta);
end

if supG
    gradG(lenJ+1,:,:) = (other.G(end,:,:) - repmat(base.G(end,:),1,1,numel(u0)))/(delta);
else
    gradG = (other.G(end,:,:) - repmat(base.G(end,:),1,1,numel(u0)))/(delta);
end

end