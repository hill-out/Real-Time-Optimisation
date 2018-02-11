function [var_phi, var_g1, var_g2] = convApprox
% Fits a convex model to phi and g in terms of y
% ----------------------------------------------
% uRange        Range of u values used as a cell
% 
% alpha         Linear approx term
% Q             Quadratic approx term
% ----------------------------------------------

%% set up range to make approximation
uRange = {linspace(2,10,11),linspace(8.5,22,11),linspace(78,92,11)};

%% find optimum point of model
optionOpt = optimoptions('fmincon','Display','off');
uGuess = [uRange{1}(6), uRange{2}(6), uRange{3}(6)];

[u_opt, phi_opt] = fmincon(@(u)(phi_from_Xu(CSTRmodel(u),u)),uGuess,[],...
    [],[],[],[0,0,50],[200,600,200],@(u)(nonLinearCon(u)));

X_opt = CSTRmodel(u_opt);
y_opt = y_from_Xu(X_opt, u_opt);
g_opt = g_from_Xu(X_opt, u_opt);

%% run lsr
[u,y,phi,g] = CSTRrange(uRange); % run CSTR at all points [n x nu1 x nu2 x nu3]

nU = size(u);
nU = nU(2:end);

dy = bsxfun(@minus, y, y_opt'); % (y-y*)
size_y = size(y,1);

[var_phi] = fminsearch(@(x)(leastSquare(x,phi_opt,phi)),[1,1,2e4,1]);
[var_g1] = fminsearch(@(x)(leastSquare(x,g_opt(1),reshape(g(1,:),[1,nU]))),[1,1,1,1]);
[var_g2] = fminsearch(@(x)(leastSquare(x,g_opt(2),reshape(g(2,:),[1,nU]))),[1,1,1,1]);

% [var_phi] = fminsearch(@(x)(leastSquare(x,phi_opt,phi)),[var_phi(1:3),0,0,var_phi(4)]);
% [var_g1] = fminsearch(@(x)(leastSquare(x,g_opt(1),reshape(g(1,:),[1,nU]))),[var_g1(1:3),0,0,var_g1(4)]);
% [var_g2] = fminsearch(@(x)(leastSquare(x,g_opt(2),reshape(g(2,:),[1,nU]))),[var_g2(1:3),0,0,var_g2(4)]);

var_phi = [phi_opt,var_phi(1:3),0,0,var_phi(4),y_opt];
var_g1 = [g_opt(1),var_g1(1:3),0,0,var_g1(4),y_opt];
var_g2 = [g_opt(2),var_g2(1:3),0,0,var_g2(4),y_opt];



%% functions

    function [c,ceq] = nonLinearCon(u)
        c = g_from_Xu(CSTRmodel(u),u);
        ceq = [];
    end

    function square = leastSquare(var,opt,fun_val)
        
        fun_c = zeros(size(phi)); %  1 x nu1 x nu2 x nu3
        a = var(1:size_y)';
        q = var(size_y+1:end)';
        
        if numel(q) == numel(a)^2
            q = reshape(q,numel(a),numel(a));        
            if any(eig(q) < 0)
                q = q*0;
            end
        else
            q(q<0) = 0;
        end
        

        
        A = repmat(a,[1,nU]);
        lin_c = sum(A.*dy,1); %linear term
        
        if all(size(q) == [size_y, size_y]) %square provided
            Q = q;
        elseif size(q,1) == size_y || size(q,2) == size_y %row/col provided
            Q = diag(q);
        else %just use first q value
            Q = eye(size_y).*q(1);
        end
        
        dy_long = reshape(dy,size_y,prod(nU)); %reshape into 2D
        quad_c = 0.5*dy_long'*Q*dy_long;
        quad_c = quad_c(logical(eye(prod(nU))));
        quad_c = reshape(quad_c,[1,nU]);
        
        fun_c = repmat(opt,size(fun_c)) + lin_c + quad_c;
        
        square = (fun_c-fun_val).^2;
        
        while numel(square) ~= 1
            square = sum(square);
        end
        
        
    end
end