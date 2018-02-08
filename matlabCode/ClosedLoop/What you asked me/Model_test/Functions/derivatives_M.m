function [dphi_du, dg_du, F_P0, phi_u0] = derivatives_M(upk_s, upk, delta_M, PBstruct, Theta_M, F_P_demand, A_us, b_u)
% ck = ck{k}
% upk = upk{k}

% To compute outside
%      if norm(dh_dc(:,i)) ~= 0
%         d{k,i} = dh_dc(:,i)/norm(dh_dc(:,i));
%      else
%         d{k,i} = dh_dc(:,i)*0;
%      end

x0              = PBstruct.x0; 
fsolve_options  = PBstruct.fsolve_options;

%% At the point
[~, phi_u0, ~, F_P0] = SimModel(upk, x0, F_P_demand, Theta_M, fsolve_options);


%%
nu = length(upk_s);

for i = 1:nu
     u_perturbed = A_us\(upk_s+[zeros(i-1,1); delta_M(i) ; zeros(nu-i,1)]) + b_u;
     [~, phi_u, ~, F_P] = SimModel(u_perturbed, x0, F_P_demand, Theta_M, fsolve_options);
     phi_x(i) = phi_u;
     g_x(i)   = F_P;
     
     dphi_du(i,1) = ( phi_x(i)-phi_u0)/norm(u_perturbed-upk);
     dg_du(i,1)    = (g_x(i)-F_P0)/norm(u_perturbed-upk);
end
end