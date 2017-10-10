function [c] = dyn_CSTR(reactOrder, kVal, c0, cIn, volFlow, V, nSteps, tau)


% Check inputs
kSize = size(kVal);
rSize = size(reactOrder);
c0Size = size(c0);
cInSize = size(cIn);
vSize = size(volFlow);

if kSize(1) ~= rSize(1) || rSize(2) ~= c0Size(2) || c0Size(2) ~= vSize(2) || ~all(c0Size == cInSize)
    error('Function input dimensional mismatch')
end

% calc number of steps
tStep = (tau/nSteps);
c = zeros(nSteps, c0Size(2));
c(1,:) = c0;
v0 = sum(volFlow);
nReact = size(reactOrder);

% run for 1:steps
for i = 1:nSteps
    nMol = c(i,:)*V;
    nLeave = c(i,:)*v0*tStep;
    nIn = cIn.*volFlow*tStep;
    
    cMat = repmat(c(i,:),nReact(1),1);
    elemOrder = reactOrder;
    elemOrder(elemOrder>0)=0;
    baseM = cMat.^(-elemOrder);
    rReact = kVal.*prod(baseM,2);
    r = rReact'*reactOrder;
    
    nGen = r*V*tStep;
    
    nNew = nMol + nIn - nLeave + nGen;
    c(i+1,:) = nNew/V;
end


end