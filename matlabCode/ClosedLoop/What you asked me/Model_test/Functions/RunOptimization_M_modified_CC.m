function [c_x_opt, f, eflag, outpt, F_P] = RunOptimization_M_modified_CC(c_x_0, ck_s, upk, F_P_demand, F_A, F_B, T_R, alpha_k, Theta, A_cs, b_c, Modified_Opt_parameters, opts, Constrain)

% x = [A_R; B_R; C_R; E_R; G_R; P_R; ...                                   // outputs WOreactorDyn          (1:6)
%      A_E; B_E; C_E; E_E; P_E; ...                                        // outputs DecanterDyn           (7:11)
%      A_C; B_C; C_C; E_C; P_C;     A_S; B_S; C_S; E_S; P_S];              // outputs DistillationColumnDyn (12:21)
%
% u = [F_A; F_B; F_C; F_E; F_P; F_G; T_R; alpha]; 
%
% c = [F_A; F_B; T_R; alpha];

%% %%%%%%%%%%%%%%%%%%%%%%%%%%

u_max  = Modified_Opt_parameters.u_max;
u_min  = Modified_Opt_parameters.u_min;
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

inv_A_cs = inv(A_cs);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%

%% NoMismatch
A = [ dh_dc*inv_A_cs, zeros(nu,nx);   % u_max
     -dh_dc*inv_A_cs, zeros(nu,nx)];  % u_min
 
b = [ u_max-upk+dh_dc*inv_A_cs*ck_s;
     -(u_min-upk+dh_dc*inv_A_cs*ck_s)];

Aeq = [];
beq = [];
if Constrain.FA == 1
Aeq = [Aeq; [1,0,0,0]*inv_A_cs, zeros(1,nx)];
beq = [beq; F_A-b_c(1)];
end
if Constrain.FB == 1
Aeq = [Aeq; [0,1,0,0]*inv_A_cs, zeros(1,nx)];
beq = [beq; F_B-b_c(2)];
end
if Constrain.TR == 1
Aeq = [Aeq; [0,0,1,0]*inv_A_cs, zeros(1,nx)];
beq = [beq; T_R-b_c(3)];
end
if Constrain.alpha == 1
Aeq = [Aeq; [0,0,0,1]*inv_A_cs, zeros(1,nx)];
beq = [beq; alpha_k-b_c(4)];
end
%%

ub = ones(1, nc+nx);
lb = zeros(1,nc+nx);

Last_c_x  = [];
Last_f    = [];
Last_c    = []; 
Last_ceq  = [];

%%%%%%%%%%%%%%%  

[c_x_opt, f, eflag, outpt] = fmincon(@(c_x)objfun(c_x, F_P_demand),...
                                c_x_0, A, b, Aeq, beq, lb, ub, ...
                                @(c_x)constr(c_x, F_P_demand),opts);

    function y = objfun(c_x, F_P_demand)
        if ~isequal(c_x,Last_c_x) % Check if computation is necessary
            
            c_s = c_x(1:nc);
            x = c_x(nc+1:end);
            
            cc = inv_A_cs*c_s+b_c;
            u = upk + dh_dc*inv_A_cs*(c_s-ck_s);
                       
            dx = OpenLoopSystem(0, x, u , Theta);
            [g, phi_u, phi_c, F_P] = ux2gphi_modified(u, cc, x, F_P_demand, Theta);
            
            Last_c_x  = c_x;
            F_P_constraint = g + epsilon_g_k + Lambda_g_k'*inv_A_cs*(c_s-ck_s);% + 1e4*1/2*(c_s-ck_s)'*(c_s-ck_s);
            
            Last_f   = phi_u + phi_c + epsilon_phi_k + lambda_phi_k'*inv_A_cs*(c_s-ck_s);% + 1e7*1/2*(c_s-ck_s)'*(c_s-ck_s);
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

    function [c,ceq] = constr(c_x, F_P_demand)
        if ~isequal(c_x,Last_c_x) % Check if computation is necessary
 
            c_s = c_x(1:nc);
            x = c_x(nc+1:end);
            
            cc = inv_A_cs*c_s+b_c;
            u = upk + dh_dc*inv_A_cs*(c_s-ck_s);
                       
            dx = OpenLoopSystem(0, x, u , Theta);
            [g, phi_u, phi_c, F_P] = ux2gphi_modified(u, cc, x, F_P_demand, Theta);
            
            Last_c_x  = c_x;
            F_P_constraint = g + epsilon_g_k + Lambda_g_k'*inv_A_cs*(c_s-ck_s);% + 1e4*1/2*(c_s-ck_s)'*(c_s-ck_s);
            
            Last_f   = phi_u + phi_c + epsilon_phi_k + lambda_phi_k'*inv_A_cs*(c_s-ck_s);% + 1e7*1/2*(c_s-ck_s)'*(c_s-ck_s);
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