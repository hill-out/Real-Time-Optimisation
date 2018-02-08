function [dphi_dc, dgp_dc, dh_dc, upk, F_Pk, phi_pk, phi_upk, upk_NoNoise, F_Pk_NoNoise, phi_upk_NoNoise] = derivatives_P(ck_s, ck, delta_P, PBstruct, Theta_P, F_P_demand, A_cs, b_c, pc, mu, horizon)

x0              = PBstruct.x0; 
fsolve_options  = PBstruct.fsolve_options;

value = 1;



%% At the point
[upk, ~, phi_upk, ~, phi_pk, ~, F_Pk] = SimPlant(ck, x0, F_P_demand, Theta_P, fsolve_options);

upk_NoNoise     = upk;
phi_upk_NoNoise = phi_upk;
F_Pk_NoNoise    = F_Pk;

error_upk1  = 0;
error_upk2  = 0;
error_upk3  = 0;
error_upk4  = 0;
error_upk5  = 0;
error_upk6  = 0;
error_phiup = 0;
error_F_Pk  = 0;
for i = 1:horizon
%     error_upk1  = error_upk1  + min(max(normrnd(mu(1), abs(upk(1)*pc(1)/100)), -abs(upk(1)*pc(1)*value/100)), abs(upk(1)*pc(1)*value/100));
%     error_upk2  = error_upk2  + min(max(normrnd(mu(2), abs(upk(2)*pc(2)/100)), -abs(upk(2)*pc(2)*value/100)), abs(upk(2)*pc(2)*value/100));
%     error_upk3  = error_upk3  + min(max(normrnd(mu(3), abs(upk(3)*pc(3)/100)), -abs(upk(3)*pc(3)*value/100)), abs(upk(3)*pc(3)*value/100));
%     error_upk4  = error_upk4  + min(max(normrnd(mu(4), abs(upk(4)*pc(4)/100)), -abs(upk(4)*pc(4)*value/100)), abs(upk(4)*pc(4)*value/100));
%     error_upk5  = error_upk5  + min(max(normrnd(mu(5), abs(upk(5)*pc(5)/100)), -abs(upk(5)*pc(5)*value/100)), abs(upk(5)*pc(5)*value/100));
%     error_upk6  = error_upk6  + min(max(normrnd(mu(6), abs(upk(6)*pc(6)/100)), -abs(upk(6)*pc(6)*value/100)), abs(upk(6)*pc(6)*value/100));
%     error_phiup = error_phiup + min(max(normrnd(mu(7), abs(phi_pk*pc(7)/100)), -abs(phi_pk*pc(7)*value/100)), abs(phi_pk*pc(7)*value/100));
%     error_F_Pk  = error_F_Pk  + min(max(normrnd(mu(8), abs(F_Pk*pc(8)/100)), -abs(F_Pk*pc(8)*value/100)), abs(F_Pk*pc(8)*value/100));

    error_upk1  = error_upk1  + normrnd(mu(1), abs(upk(1)*pc(1)/100));
    error_upk1  = min([error_upk1,  3*abs(upk(1)*pc(1)/100)]);
    error_upk1  = max([error_upk1, -3*abs(upk(1)*pc(1)/100)]);
    
    error_upk2  = error_upk2  + normrnd(mu(2), abs(upk(2)*pc(2)/100));
    error_upk2  = min([error_upk2,  3*abs(upk(2)*pc(2)/100)]);
    error_upk2  = max([error_upk2, -3*abs(upk(2)*pc(2)/100)]);
    
    error_upk3  = error_upk3  + normrnd(mu(3), abs(upk(3)*pc(3)/100));
    error_upk3  = min([error_upk3,  3*abs(upk(3)*pc(3)/100)]);
    error_upk3  = max([error_upk3, -3*abs(upk(3)*pc(3)/100)]);
    
    
    error_upk4  = error_upk4  + normrnd(mu(4), abs(upk(4)*pc(4)/100));
    error_upk4  = min([error_upk4,  3*abs(upk(4)*pc(4)/100)]);
    error_upk4  = max([error_upk4, -3*abs(upk(4)*pc(4)/100)]);
    
    error_upk5  = 0;
    
    error_upk6  = error_upk6  + normrnd(mu(6), abs(upk(6)*pc(6)/100));
    error_upk6  = min([error_upk6,  3*abs(upk(6)*pc(6)/100)]);
    error_upk6  = max([error_upk6, -3*abs(upk(6)*pc(6)/100)]);
    
    
    error_phiup = error_phiup + normrnd(mu(7), abs(phi_pk*pc(7)/100));
    error_phiup  = min([error_phiup,  3*abs(phi_pk*pc(7)/100)]);
    error_phiup  = max([error_phiup, -3*abs(phi_pk*pc(7)/100)]);
    
    error_F_Pk  = error_F_Pk  + normrnd(mu(8), abs(F_Pk*pc(8)/100));
    error_F_Pk  = min([error_F_Pk,  3*abs(F_Pk*pc(8)/100)]);
    error_F_Pk  = max([error_F_Pk, -3*abs(F_Pk*pc(8)/100)]);

end
error_upk1  = error_upk1/horizon;
error_upk2  = error_upk2/horizon;
error_upk3  = error_upk3/horizon;
error_upk4  = error_upk4/horizon;
error_upk5  = error_upk5/horizon;
error_upk6  = error_upk6/horizon;
error_phiup = error_phiup/horizon;
error_F_Pk  = error_F_Pk/horizon;

upk = upk + [error_upk1; error_upk2; error_upk3; error_upk4; error_upk5; error_upk6; 0; 0];
F_Pk    = F_Pk  + error_F_Pk;
phi_upk = phi_upk + error_phiup;

%%
nc = length(ck_s);

% dphi_dc = zeros(nc,1);
% dgp_dc  = zeros(nc,1); 
% dh_dc   = zeros(8,nc);

for i = 1:nc
     ck_perturbed_p = A_cs\(ck_s+[zeros(i-1,1); delta_P(i) ; zeros(nc-i,1)]) + b_c;
     [up, ~, phi, ~, ~, ~, F_P] = SimPlant(ck_perturbed_p, x0, F_P_demand, Theta_P, fsolve_options);

     error_up1  = 0;
     error_up2  = 0;
     error_up3  = 0;
     error_up4  = 0;
     error_up5  = 0;
     error_up6  = 0;
     error_phi  = 0;
     error_F_P  = 0;
     for j = 1:horizon
         error_up1  = error_up1 + min(max(normrnd(mu(1), abs(upk(1)*pc(1)/100)), -abs(upk(1)*pc(1)*value/100)), abs(upk(1)*pc(1)*value/100));
         error_up2  = error_up2 + min(max(normrnd(mu(2), abs(upk(2)*pc(2)/100)), -abs(upk(2)*pc(2)*value/100)), abs(upk(2)*pc(2)*value/100));
         error_up3  = error_up3 + min(max(normrnd(mu(3), abs(upk(3)*pc(3)/100)), -abs(upk(3)*pc(3)*value/100)), abs(upk(3)*pc(3)*value/100));
         error_up4  = error_up4 + min(max(normrnd(mu(4), abs(upk(4)*pc(4)/100)), -abs(upk(4)*pc(4)*value/100)), abs(upk(4)*pc(4)*value/100));
         error_up5  = error_up5 + min(max(normrnd(mu(5), abs(upk(5)*pc(5)/100)), -abs(upk(5)*pc(5)*value/100)), abs(upk(5)*pc(5)*value/100));
         error_up6  = error_up6 + min(max(normrnd(mu(6), abs(upk(6)*pc(6)/100)), -abs(upk(6)*pc(6)*value/100)), abs(upk(6)*pc(6)*value/100));
         error_phi  = error_phi + min(max(normrnd(mu(7), abs(phi_pk*pc(7)/100)), -abs(phi_pk*pc(7)*value/100)), abs(phi_pk*pc(7)*value/100));
         error_F_P  = error_F_P + min(max(normrnd(mu(8), abs(F_Pk*pc(8)/100)), -abs(F_Pk*pc(8)*value/100)), abs(F_Pk*pc(8)*value/100));
     end
     error_up1  = error_up1/horizon;
     error_up2  = error_up2/horizon;
     error_up3  = error_up3/horizon;
     error_up4  = error_up4/horizon;
     error_up5  = error_up5/horizon;
     error_up6  = error_up6/horizon;
     error_phi  = error_phi/horizon;
     error_F_P  = error_F_P/horizon;
     
     phi_x_p(i) = phi + error_phi;
     g_x_p(i)   = F_P + error_F_P;
     up_x_p(:,i)  = up + [error_up1; error_up2; error_up3; error_up4; error_up5; error_up6; 0; 0];
     
     
     
     ck_perturbed_m = A_cs\(ck_s+[zeros(i-1,1); -delta_P(i) ; zeros(nc-i,1)]) + b_c;
     [up, ~, phi, ~, ~, ~, F_P] = SimPlant(ck_perturbed_m, x0, F_P_demand, Theta_P, fsolve_options);

     error_up1  = 0;
     error_up2  = 0;
     error_up3  = 0;
     error_up4  = 0;
     error_up5  = 0;
     error_up6  = 0;
     error_phi  = 0;
     error_F_P  = 0;
     for j = 1:horizon
         error_up1  = error_up1 + min(max(normrnd(mu(1), abs(upk(1)*pc(1)/100)), -abs(upk(1)*pc(1)*value/100)), abs(upk(1)*pc(1)*value/100));
         error_up2  = error_up2 + min(max(normrnd(mu(2), abs(upk(2)*pc(2)/100)), -abs(upk(2)*pc(2)*value/100)), abs(upk(2)*pc(2)*value/100));
         error_up3  = error_up3 + min(max(normrnd(mu(3), abs(upk(3)*pc(3)/100)), -abs(upk(3)*pc(3)*value/100)), abs(upk(3)*pc(3)*value/100));
         error_up4  = error_up4 + min(max(normrnd(mu(4), abs(upk(4)*pc(4)/100)), -abs(upk(4)*pc(4)*value/100)), abs(upk(4)*pc(4)*value/100));
         error_up5  = error_up5 + min(max(normrnd(mu(5), abs(upk(5)*pc(5)/100)), -abs(upk(5)*pc(5)*value/100)), abs(upk(5)*pc(5)*value/100));
         error_up6  = error_up6 + min(max(normrnd(mu(6), abs(upk(6)*pc(6)/100)), -abs(upk(6)*pc(6)*value/100)), abs(upk(6)*pc(6)*value/100));
         error_phi  = error_phi + min(max(normrnd(mu(7), abs(phi_pk*pc(7)/100)), -abs(phi_pk*pc(7)*value/100)), abs(phi_pk*pc(7)*value/100));
         error_F_P  = error_F_P + min(max(normrnd(mu(8), abs(F_Pk*pc(8)/100)), -abs(F_Pk*pc(8)*value/100)), abs(F_Pk*pc(8)*value/100));
     end
     error_up1  = error_up1/horizon;
     error_up2  = error_up2/horizon;
     error_up3  = error_up3/horizon;
     error_up4  = error_up4/horizon;
     error_up5  = error_up5/horizon;
     error_up6  = error_up6/horizon;
     error_phi  = error_phi/horizon;
     error_F_P  = error_F_P/horizon;
     
     phi_x_m(i) = phi + error_phi;
     g_x_m(i)   = F_P + error_F_P;
     up_x_m(:,i)  = up + [error_up1; error_up2; error_up3; error_up4; error_up5; error_up6; 0; 0];
     
     
     
%      % Central
%      dphi_dc(i,1) = ( phi_x_p(i)-phi_x_m(i))/(2*norm(ck_perturbed_p-ck));
%      dgp_dc(i,1)  = (g_x_p(i)-g_x_m(i))/(2*norm(ck_perturbed_p-ck));
%      dh_dc(:,i)   = (up_x_p(:,i)-up_x_m(:,i))/(2*norm(ck_perturbed_p-ck));

     % Forward
     dphi_dc(i,1) = ( phi_x_p(i)-phi_upk)/(norm(ck_perturbed_p-ck));
     dgp_dc(i,1)  = (g_x_p(i)-F_Pk)/(norm(ck_perturbed_p-ck));
     dh_dc(:,i)   = (up_x_p(:,i)-upk)/(norm(ck_perturbed_p-ck));
end



end