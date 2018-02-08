function [u_c_x_opt, f, eflag, outpt, F_P] = RunOptimization_M_modified_UC(u_x_0, upk_s, ck, F_P_demand, Theta, A_us, b_u, Modified_Opt_parameters, opts, Constrain)

% x = [A_R; B_R; C_R; E_R; G_R; P_R; ...                                   // outputs WOreactorDyn          (1:6)
%      A_E; B_E; C_E; E_E; P_E; ...                                        // outputs DecanterDyn           (7:11)
%      A_C; B_C; C_C; E_C; P_C;     A_S; B_S; C_S; E_S; P_S];              // outputs DistillationColumnDyn (12:21)
%
% u = [F_A; F_B; F_C; F_E; F_P; F_G; T_R; alpha]; 
%
% c = [F_A; F_B; T_R; alpha];

%% %%%%%%%%%%%%%%%%%%%%%%%%%%

c_max  = Modified_Opt_parameters.c_max;
c_min  = Modified_Opt_parameters.c_min;
nu     = Modified_Opt_parameters.nu;
nc     = Modified_Opt_parameters.nc;
nx     = Modified_Opt_parameters.nx;
Nk     = Modified_Opt_parameters.Nk;
dh_dc  = Modified_Opt_parameters.dh_dc;
pinv_dh_dc    = Modified_Opt_parameters.pinv_dh_dc;
epsilon_phi_k = Modified_Opt_parameters.epsilon_phi_k;
lambda_phi_k  = Modified_Opt_parameters.lambda_phi_k;
epsilon_g_k   = Modified_Opt_parameters.epsilon_g_k;
Lambda_g_k    = Modified_Opt_parameters.Lambda_g_k;




inv_A_us = inv(A_us);



%% %%%%%%%%%%%%%%%%%%%%%%%%%%

Inu= eye(nu);


%% NoMismatch
A = [ pinv_dh_dc*inv_A_us, zeros(nc,nx);  % c_max
     -pinv_dh_dc*inv_A_us, zeros(nc,nx)];  % c_min
%       [0,0,1,0]*pinv_dh_dc*inv_A_us, zeros(1,nx);  % DTR_max sup
%      -[0,0,1,0]*pinv_dh_dc*inv_A_us, zeros(1,nx)];  % DTR_max lo
 
b = [ c_max-ck+pinv_dh_dc*inv_A_us*upk_s;
     -(c_min-ck+pinv_dh_dc*inv_A_us*upk_s)];
%        1000+[0,0,1,0]*pinv_dh_dc*inv_A_us*upk_s;
%        1000-[0,0,1,0]*pinv_dh_dc*inv_A_us*upk_s];

Aeq = [Nk'*inv_A_us, zeros(size(Nk',1),nx)];
beq = Nk'*inv_A_us*upk_s;
 
% if Constrain.FA == 1
% Aeq = [Aeq; [1,zeros(1,nc-1)]*pinv_dh_dc*inv_A_us, zeros(1,nx) ];
% beq = [beq; F_A-ck(1)+[1,zeros(1,nc-1)]*pinv_dh_dc*inv_A_us*upk_s];
% end
% if Constrain.FB == 1
% Aeq = [Aeq; [0,1,zeros(1,nc-2)]*pinv_dh_dc*inv_A_us, zeros(1,nx) ];
% beq = [beq; F_B-ck(2)+[0,1,zeros(1,nc-2)]*pinv_dh_dc*inv_A_us*upk_s];
% end
% if Constrain.TR == 1
% Aeq = [Aeq; [0,0,1,zeros(1,nc-3)]*pinv_dh_dc*inv_A_us, zeros(1,nx) ];
% beq = [beq; T_R-ck(3)+[0,0,1,zeros(1,nc-3)]*pinv_dh_dc*inv_A_us*upk_s];
% end
% if Constrain.alpha == 1
% Aeq = [Aeq; [0,0,0,1,zeros(1,nc-4)]*pinv_dh_dc*inv_A_us, zeros(1,nx) ];
% beq = [beq; T_R-ck(3)+[0,0,0,1,zeros(1,nc-4)]*pinv_dh_dc*inv_A_us*upk_s];
% end

%%

ub = ones(1, nu+nx);
lb = zeros(1,nu+nx);

Last_u_x  = [];
Last_f    = [];
Last_c    = []; 
Last_ceq  = [];

%%%%%%%%%%%%%%%  

[u_c_x_opt, f, eflag, outpt] = fmincon(@(u_x)objfun(u_x, F_P_demand),...
                                u_x_0, A, b, Aeq, beq, lb, ub, ...
                                @(u_x)constr(u_x, F_P_demand),opts);

    function y = objfun(u_x, F_P_demand)
        if ~isequal(u_x,Last_u_x) % Check if computation is necessary
            
            u_s = u_x(1:nu);
            x = u_x(nu+1:end);
            
            cc = ck + pinv_dh_dc*inv_A_us*(u_s-upk_s);
            
            u = inv_A_us*u_s+b_u;
                       
            dx = OpenLoopSystem(0, x, u , Theta);
            [g, phi_u, phi_c, F_P] = ux2gphi_modified(u, cc, x, F_P_demand, Theta);
            
            Last_u_x  = u_x;
            F_P_constraint = g + epsilon_g_k + Lambda_g_k'*(cc-ck);
            
            Last_f   = phi_u + phi_c + epsilon_phi_k + lambda_phi_k'*(cc-ck) ;
            Last_c   = [];
            if Constrain.FP == 1
                Last_ceq = [F_P_constraint; dx]; 
            else
                Last_ceq = dx;
            end
        end
        % Now compute objective function
        y = Last_f;
    end

    function [c,ceq] = constr(u_x, F_P_demand)
        if ~isequal(u_x,Last_u_x) % Check if computation is necessary
            
            u_s = u_x(1:nu);
            x = u_x(nu+1:end);
            
            cc = ck + pinv_dh_dc*inv_A_us*(u_s-upk_s);
            
            u = inv_A_us*u_s+b_u;
                       
            dx = OpenLoopSystem(0, x, u , Theta);
            [g, phi_u, phi_c, F_P] = ux2gphi_modified(u, cc, x, F_P_demand, Theta);
            
            Last_u_x  = u_x;
           F_P_constraint = g + epsilon_g_k + Lambda_g_k'*(cc-ck);
            
            Last_f   = phi_u + phi_c + epsilon_phi_k + lambda_phi_k'*(cc-ck) ;
            Last_c   = [];
            
            if Constrain.FP == 1
                
                Last_ceq = [F_P_constraint; dx]; 
            else
                Last_ceq = dx;
            end
        end
        % Now compute constraint functions
        c   = Last_c;
        ceq = Last_ceq;
    end

end