function [cx, f, eflag, outpt] = RunOptimization_P(cx0_s, F_P_demand, F_A, F_B, T_R, alpha_k, Theta, opts, Constrain, inv_A_cs, b_c)

Last_cx  = [];
Last_f   = [];
Last_c   = []; % Attention c is not the input of the plant but the NL constraint !!
Last_ceq = [];

nc = 4;
nx = 21;

% Constrain.FA
% Constrain.FB
% Constrain.TR
% Constrain.alpha

%%%%%%%%%%%%%%% 

lb = [  0;   0; 0; 0; zeros(21,1)];
ub = [1; 1; 1; 1; ...
      1; 1; 1; 1; 1; 1;   1; 1; 1; 1; 1;     1; 1; 1; 1; 1;     1; 1; 1; 1; 1];
  

A = [];
b = [];

Aeq = [];
beq = [];

if Constrain.FA == 1
Aeq = [Aeq; [1, zeros(1,nc-1)]*inv_A_cs, zeros(1,nx)];
beq = [beq; F_A-[1, zeros(1,nc-1)]*b_c];
end

if Constrain.FB == 1
Aeq = [Aeq; [0,1, zeros(1,nc-2)]*inv_A_cs, zeros(1,nx)];
beq = [beq; F_B-[0, 1, zeros(1,nc-2)]*b_c];
end

if Constrain.TR == 1
Aeq = [Aeq; [0,0,1, zeros(1,nc-3)]*inv_A_cs, zeros(1,nx)];
beq = [beq; T_R-[0,0, 1, zeros(1,nc-3)]*b_c];
end

if Constrain.alpha == 1
Aeq = [Aeq; [0,0,0,1, zeros(1,nc-4)]*inv_A_cs, zeros(1,nx)];
beq = [beq; alpha_k-[0,0,0, 1, zeros(1,nc-4)]*b_c];
end


%%%%%%%%%%%%%%%  

[cx, f, eflag, outpt] = fmincon(@(cx)objfun(cx, F_P_demand),...
                                  cx0_s, A, b, Aeq, beq, lb, ub, ...
                                @(cx)constr(cx, F_P_demand),opts);

    function y = objfun(cx,F_P_demand)
        if ~isequal(cx,Last_cx) % Check if computation is necessary
            cc_s = cx(1:4);
            xx = cx(5:end);
            dx = ClosedLoopSystem(0, xx, inv_A_cs*cc_s+b_c, Theta);
            [~, g, ~, ~, phi_p] = cx2ugphi(inv_A_cs*cc_s+b_c, xx, F_P_demand, Theta);
            Last_cx  = cx;
            Last_f   = phi_p;
            Last_c   = [];  
            if Constrain.FP == 1
                Last_ceq = [g; dx];
            else
                Last_ceq = dx;
            end
        end
        % Now compute objective function
        y = Last_f;
    end

    function [c,ceq] = constr(cx,F_P_demand)
        if ~isequal(cx,Last_cx) % Check if computation is necessary
            cc_s = cx(1:4);
            xx = cx(5:end);
            dx = ClosedLoopSystem(0, xx, inv_A_cs*cc_s+b_c, Theta);
            [~, g, ~, ~, phi_p] = cx2ugphi(inv_A_cs*cc_s+b_c, xx, F_P_demand, Theta);
            Last_cx  = cx;
            Last_f   = phi_p;
            Last_c   = [];  
            if Constrain.FP == 1
                Last_ceq = [g; dx];
            else
                Last_ceq = dx;
            end
        end
        % Now compute constraint functions
        c   = Last_c;
        ceq = Last_ceq;
    end

end