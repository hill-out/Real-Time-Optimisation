function [dfun] = openModelGrad(u)
% Calculates the gradients of C at u
% ----------------------------------
% u         Inputs to the model
%
% dfun      Sturcture of gradients
% ----------------------------------
global dudu
baseC = openModel(u);
basePhi = phiFun(u,baseC);
baseg1 = g1Fun(u,baseC);
baseg2 = g2Fun(u,baseC);

T = [8075, 12400];
dT = diag(T)/1000;
dU = diag(u)/1000;

for i = 1:2
    % dCdT
    newC = openModel(u,baseC,dT(i,:));
    dCdT(i,:) = (newC - baseC)/dT(i,i);
    
    % dPhidT
    newPhi = phiFun(u,newC);
    dPhidT(i,:) = (newPhi - basePhi)/dT(i,i);
    
    % dg1dT
    newg1 = g1Fun(u,newC);
    dg1dT(i,:) = (newg1 - baseg1)/dT(i,i);    
    
    % dg2dT
    newg2 = g2Fun(u,newC);
    dg2dT(i,:) = (newg2 - baseg2)/dT(i,i);  
end

for j = 1:numel(u)
    % dCdu
    newC = openModel(u+dU(j,:),baseC);
    dCdu(j,:) = (newC - baseC)/dU(j,j);
    
    % dPhidu
    newPhi = phiFun(u+dU(j,:),newC);
    dPhidu(j,:) = (newPhi - basePhi)/dU(j,j);
    
    % dg1du
    newg1 = g1Fun(u+dU(j,:),newC);
    dg1du(j,:) = (newg1 - baseg1)/dU(j,j);    
    
    % dg2du
    newg2 = g2Fun(u+dU(j,:),newC);
    dg2du(j,:) = (newg2 - baseg2)/dU(j,j);    
end

for k = 1:numel(u)
    newCk = openModel(u+dU(k,:),baseC);
    newPhik = phiFun(u+dU(k,:),newCk);
    newg1k = g1Fun(u+dU(k,:),newCk);
    newg2k = g2Fun(u+dU(k,:),newCk);
    
    for j = 1:numel(u)
        newCj = openModel(u+dU(k,:)+dU(j,:),baseC);
        newPhij = phiFun(u+dU(k,:)+dU(j,:),newCj);
        newg1j = g1Fun(u+dU(k,:)+dU(j,:),newCj);
        newg2j = g2Fun(u+dU(k,:)+dU(j,:),newCj);
        
        dPhidudu(j,:) = (newPhij - newPhik)/dU(j,j);
        dg1dudu(j,:) = (newg1j - newg1k)/dU(j,j);
        dg2dudu(j,:) = (newg2j - newg2k)/dU(j,j);        
    end
    
    ddPhidudu(k,:) = (dPhidudu - dPhidu)/dU(k,k);
    ddg1dudu(k,:) = (dg1dudu - dg1du)/dU(k,k);
    ddg2dudu(k,:) = (dg2dudu - dg2du)/dU(k,k);
    
    for i = 1:2
        newCi = openModel(u+dU(k,:),baseC,dT(i,:));
        newPhii = phiFun(u+dU(k,:),newCi);
        newg1i = g1Fun(u+dU(k,:),newCi);
        newg2i = g2Fun(u+dU(k,:),newCi);
        
        dPhidudT(i,:) = (newPhii - newPhik)/dT(i,i);
        dg1dudT(i,:) = (newg1i - newg1k)/dT(i,i);
        dg2dudT(i,:) = (newg2i - newg2k)/dT(i,i);        
    end
    
    ddPhidudT(k,:) = (dPhidudT - dPhidT)/dU(k,k);
    ddg1dudT(k,:) = (dg1dudT - dg1dT)/dU(k,k);
    ddg2dudT(k,:) = (dg2dudT - dg2dT)/dU(k,k);
end

% output
dfun.ddphidudu = ddPhidudu*dudu*dudu;
dfun.ddg1dudu = ddg1dudu*dudu*dudu;
dfun.ddg2dudu = ddg2dudu*dudu*dudu;

dfun.ddphidudT = (ddPhidudT'*dudu)';
dfun.ddg1dudT = (ddg1dudT'*dudu)';
dfun.ddg2dudT = (ddg2dudT'*dudu)';

dfun.dCdu = (dCdu'*dudu)';
dfun.dCdT = dCdT;

% dfun.dCdu = [zeros(3,1), ones(3,1)];
% dfun.dCdT = [dCdT(:,1), zeros(2,1)];
end