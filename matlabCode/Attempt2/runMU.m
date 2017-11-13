function [a] = runMU()
% runs the multiUnitTransRTO function
%

s = [-1,  0;
     -1, -2;
      1,  0;
      0,  1]; %stoichiometry
V = 500;                        % volume of reactor [L]

m_k = [0.75, 1.5];              % k value for model [L/(mol min)]
m_c_in = [2; 1.5; 0; 0];        % con. in for model [mol/L]

p_k = [1.4, 0.4];               % k value for plant [L/(mol min)]
p_c_in = [2.5; 1.5; 0; 0];      % con. in for plant [mol/L]
p_model = @(u, c0, t)(CSTRdyn(p_c_in, u, p_k, V, s, t, c0);


    function [phi] = CSTRcost(u, c, c_in)
        % calculates the cost (-J)
        w = 0.004;
        J = c(3)^2*(u(1)+u(2))^2/(u(1)*c_in(1)) - w*(u(1)^2+u(2)^2);
        phi = -J;
    end
end