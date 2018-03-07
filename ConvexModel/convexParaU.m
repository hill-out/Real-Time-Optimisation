function [convPara] = convexParaU(f, u, f_opt, u_opt)
% finds the parameters for the convex approximation in u
% ------------------------------------------------------
% f         Values of the function
% u         Inputs to the function
% f_opt     Optimum value
% u_opt     Optimum inputs
% 
% convPara  Parameters of the convex aproximation
% ------------------------------------------------------

% initialise convPara
convPara = [1, 1, 1, 1, 0, 0, 1, 0, 1];
uShift = bsxfun(@minus, u, u_opt);
n = size(u,1);

convPara = fmincon(@(x)lsr([f_opt, x]), convPara, [], [], [], [], [], [], @(x)conFun(x));
convPara = [f_opt, convPara];

    function [square] = lsr(x)
        fitF = convCalc(x,uShift);
        square = (fitF-f).^2;
        while numel(square)>1
            square = sum(square);
        end
    end

    function [c, ceq] = conFun(a)
        
        q = a(n+1:end);
        locat = symMat(n);
        m2 = q(locat);
        
        c = -eig(m2);
        ceq = 0;
    end
end