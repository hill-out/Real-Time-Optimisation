function [convPara] = convexParaR(f, r, f_opt, r_opt)
% finds the parameters for the convex approximation in r
% ------------------------------------------------------
% f         Values of the function
% r         Setpoints to the function
% f_opt     Optimum value
% r_opt     Optimum setpoints
% 
% convPara  Parameters of the convex aproximation
% ------------------------------------------------------

% initialise convPara
convPara = [1, 1, 1, 0, 1];
rShift = bsxfun(@minus, r, r_opt);
n = size(r,1);

convPara = fmincon(@(x)lsr([f_opt, x]), convPara, [], [], [], [], [], [], @(x)conFun(x));
convPara = [f_opt, convPara];

    function [square] = lsr(x)
        fitF = convCalc(x,rShift);
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