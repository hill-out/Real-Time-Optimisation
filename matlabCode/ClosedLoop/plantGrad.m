function [dX, dY, dPhi, dG] = plantGrad(r)
% Symbolically solves the model 2-reaction CSTR
% ---------------------------------------------
% r         Inputs [nr x 1]
% 
% dX        Gradient of X [nx x nr]
% dY        Gradient of Y [ny x nr]
% dPhi      Gradient of Phi [nphi x nr]
% dG        Gradient of G [ng x nr]
% ---------------------------------------------



n = numel(r); % nr
rMesh = meshgrid(r);
dr = 0.001; % percent change of r

rdrMesh = rMesh.*(ones(n)+eye(n)*dr); % r's with changes
drMesh = rMesh.*eye(n)*dr; % dr's

%base case
baseU = plantController(r);
baseX = CSTRplant(baseU);
baseY = y_from_Xu(baseX,baseU);
basePhi = phi_from_Xu(baseX,baseU);
baseG = g_from_Xu(baseX,baseU);

%initialise output
dX = zeros(numel(baseX),n);
dY = zeros(numel(baseY),n);
dPhi = zeros(numel(basePhi),n);
dG = zeros(numel(baseG),n);

for i = 1:n
    %run new plant
    newU = plantController(rdrMesh(i,:));
    newX = CSTRplant(newU, baseX);
    newY = y_from_Xu(newX,newU);
    newPhi = phi_from_Xu(newX,newU);
    newG = g_from_Xu(newX,newU);
    
    dX(:,i) = (newX - baseX)/drMesh(i,i);
    dY(:,i) = (newY - baseY)/drMesh(i,i);
    dPhi(:,i) = (newPhi - basePhi)/drMesh(i,i);
    dG(:,i) = (newG - baseG)/drMesh(i,i);
end

end