function solysis(rpi, base, rp_opt, phip_opt, algo, tau, K)
% solution analysis (sol-ysis)
% -----------------------------------
% rpi       Iterative plant setpoints
% base      Struct of plant time data
% rp_opt    Optimum setpoint of plant
% phip_opt  Optimum plant phi @rp_opt
% algo      Algorithm for closed loop
% tau, K    Freq and gain of solution
% -----------------------------------

% rpi analysis
a = any(abs(bsxfun(@rdivide,bsxfun(@minus,rpi,rp_opt),rp_opt))>0.05,2);
b = any(abs(bsxfun(@rdivide,bsxfun(@minus,rpi,rp_opt),rp_opt))>0.01,2);
c = any(abs(bsxfun(@rdivide,bsxfun(@minus,rpi,rp_opt),rp_opt))>0.002,2);

% phi analysis
phi = base.phip(rem(base.t,tau)==0);
phi = phi(2:2:numel(phi));

d = any(abs(bsxfun(@rdivide,bsxfun(@minus,phi,phip_opt),phip_opt))>0.05,2);
e = any(abs(bsxfun(@rdivide,bsxfun(@minus,phi,phip_opt),phip_opt))>0.01,2);
f = any(abs(bsxfun(@rdivide,bsxfun(@minus,phi,phip_opt),phip_opt))>0.002,2);

% print analysis
fprintf('-----------------------------------------------\n')
fprintf('---|  %2s  |---| Tau = %4.0f |---| K = %5.2f |---\n',algo,tau,K);
fprintf('-----------------------------------------------\n')
fprintf('  Iterations to optimum r   : [%3d, %3d, %3d]\n',max(find(a)+1),max(find(b)+1),max(find(c)+1));
fprintf('  Iterations to optimum phi : [%3d, %3d, %3d]\n',max(find(d)+1),max(find(e)+1),max(find(f)+1));
fprintf('-----------------------------------------------\n')

end