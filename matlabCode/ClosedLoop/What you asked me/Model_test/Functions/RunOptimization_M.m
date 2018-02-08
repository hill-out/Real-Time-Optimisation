function [ux, f, eflag, outpt] = RunOptimization_M(ux0, F_P_demand, Theta, opts)

Last_ux  = [];
Last_f   = [];
Last_c   = []; % Attention c is not the input of the plant but the NL constraint !!
Last_ceq = [];

%%%%%%%%%%%%%%%  

% x = [A_R; B_R; C_R; E_R; G_R; P_R; ...                                   // outputs WOreactorDyn          (1:6)
%      A_E; B_E; C_E; E_E; P_E; ...                                        // outputs DecanterDyn           (7:11)
%      A_C; B_C; C_C; E_C; P_C;     A_S; B_S; C_S; E_S; P_S];    // outputs DistillationColumnDyn (12:21)
%
% u = [F_A; F_B; F_C; F_E; F_P; F_G; T_R; alpha]; 


lb = [zeros(6,1)     ; 323; 0; zeros(21,1)];
ub = [50000*ones(6,1); 378; 1;  ones(21,1)];

Aeq = [];  
beq = [];

A = [];
b = [];

%%%%%%%%%%%%%%%  

[ux, f, eflag, outpt] = fmincon(@(ux)objfun(ux, F_P_demand),...
                                ux0, A, b, Aeq, beq, lb, ub, ...
                                @(ux)constr(ux, F_P_demand),opts);

    function y = objfun(ux,F_P_demand)
        if ~isequal(ux,Last_ux) % Check if computation is necessary
            u = ux(1:8);
            x = ux(9:end);
            dx = OpenLoopSystem(0, x, u, Theta);
            [g, phi_u, ~] = ux2gphi(u, x, F_P_demand, Theta);
            Last_ux  = ux;
            Last_f   = phi_u;
            Last_c   = [];
%             F_P_in   = u(5);
            Last_ceq = [g; dx];
        end
        % Now compute objective function
        y = Last_f;
    end

    function [c,ceq] = constr(ux,F_P_demand)
        if ~isequal(ux,Last_ux) % Check if computation is necessary
            u = ux(1:8);
            x = ux(9:end);
            dx = OpenLoopSystem(0, x, u, Theta);
            [g, phi_u, ~] = ux2gphi(u,x, F_P_demand, Theta);
            Last_ux  = ux;
            Last_f   = phi_u;
            Last_c   = [];
%             F_P_in   = u(5);
            Last_ceq = [g; dx];
        end
        % Now compute constraint functions
        c   = Last_c;
        ceq = Last_ceq;
    end

end