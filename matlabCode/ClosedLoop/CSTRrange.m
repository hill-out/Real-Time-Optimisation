function [u,y,phi,g] = CSTRrange(uRange)
% Runs the CSTR model for a range of u's
% ----------------------------------------------
% uRange        Range of u values used as a cell
% 
% y             Model outputs
% phi           Model cost
% g             Model constraints
% ----------------------------------------------

%% make grids, where grids(:,[ni...]) is the nth combination of u
uRange = {uRange{1},[0],uRange{2:end}};
[grids{1:numel(uRange)}] = ndgrid(uRange{:});
grids = [grids{:}];
grids(:,2,:,:) = [];
u = permute(grids,[2,1,3,4]);

%% test to get the size of y, phi and g
testU = grids(:,1,1,1); %first combination
testX = CSTRmodel(testU);

testY = y_from_Xu(testX,testU);
testPhi = phi_from_Xu(testX,testU);
testG = g_from_Xu(testX,testU);

%% set up variables
nU = size(u);
nU = nU(2:end);

y = zeros([numel(testY),nU]);
phi = zeros([numel(testPhi),nU]);
g = zeros([numel(testG),nU]);

%% run for all combinations of u
xi = testX;

for i = 1:prod(nU)
    ui = u(:,i);
    xi = CSTRmodel(ui,xi);
    
    y(:,i) = y_from_Xu(xi,ui);
    phi(:,i) = phi_from_Xu(xi,ui);
    g(:,i) = g_from_Xu(xi,ui);
end    

end