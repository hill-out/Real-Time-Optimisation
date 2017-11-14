function [] = runMU()
% runs the multiUnitTransRTO function
%


% parameters for both model and plant
s = [-1,  0;
     -1, -2;
      1,  0;
      0,  1]; %stoichiometry
V = 500;                        % volume of reactor [L]
dH = [3.5, 1.5];

% constraint parameters
Qmax = 110;
Dmax = 0.1;

% model parameters
m.k = [0.75, 1.5];              % k value for model [L/(mol min)]
m.c_in = [2; 1.5; 0; 0];        % conc_in for model [mol/L]
m.dynModel = @(u, c0, t)(CSTRdyn(m.c_in, u, m.k, V, s, t, c0));
m.steModel = @(u)(CSTRsteady(u,  m.c_in, m.k));

% plant parameters
p.k = [1.4, 0.4];               % k value for plant [L/(mol min)]
p.c_in = [2.5; 1.5; 0; 0];      % conc_in for plant [mol/L]
p.dynModel = @(u, c0, t)(CSTRdyn(p.c_in, u, p.k, V, s, t, c0));
p.steModel = @(u)(CSTRsteady(u,  p.c_in, p.k));

% model optimum
options_fmincon = optimoptions(@fmincon,'Display','Off');
x0 = [10, 10];
uMin = [0, 0];
uMax = [50, 50];

[m.opt] = fmincon(@(x)(calcCost(x,m.steModel,m.c_in)), x0, [], [], [],...
    [], uMin, uMax, @(x)(calcCon(x,m.steModel,m.k)), options_fmincon);

% plant optimum
[p.opt] = fmincon(@(x)(calcCost(x,p.steModel,p.c_in)), x0, [], [], [],...
    [], uMin, uMax, @(x)(calcCon(x,p.steModel,p.k)), options_fmincon);

% convex approximation of model
m.convCostPara = convexApprox(@(x)(CSTRcost(x, m.steModel(x), m.c_in)),...
    uMin+1, uMax, 31, 2, m.opt);
m.convCon1Para = convexApprox(@(x)(nE(CSTRcon(m.steModel(x), m.k),1)),...
    uMin+1, uMax, 31, 1, m.opt);
m.convCon2Para = convexApprox(@(x)(nE(CSTRcon(m.steModel(x), m.k),2)),...
    uMin+1, uMax, 31, 1, m.opt);

m.optFull = [m.opt(1); m.opt(2); 0; 0];
m.optCost = CSTRcost(m.optFull, m.steModel(m.optFull), m.c_in);
m.optCon = CSTRcon(m.steModel(m.optFull), m.k);

m.convCost = @(x)(convModel(m.convCostPara, x, m.opt', m.optCost));
m.convCon1 = @(x)(convModel([m.convCon1Para,0], x, m.opt', m.optCon(1)));
m.convCon2 = @(x)(convModel([m.convCon2Para,0], x, m.opt', m.optCon(2)));

m.convCostGrad = @(x)(convGrad(m.convCostPara, x, m.opt'));
m.convCon1Grad = @(x)(convGrad([m.convCon1Para,0], x, m.opt'));
m.convCon2Grad = @(x)(convGrad([m.convCon2Para,0], x, m.opt'));

m.conv = {m.convCost, m.convCon1, m.convCon2};
m.cGrad = {m.convCostGrad, m.convCon1Grad, m.convCon2Grad};

p.c0 = p.steModel(m.optFull);

p.dynFunc = {p.dynModel,...
    @(x,c0,t)(CSTRcost(x, p.dynModel(x,c0,t), p.c_in)),...
    @(x,c0,t)(nE(CSTRcon(p.dynModel(x,c0,t), p.k),1)),...
    @(x,c0,t)(nE(CSTRcon(p.dynModel(x,c0,t), p.k),2))};

multiUnitTransRTO(m.conv, m.cGrad, p.dynFunc, 0.5, m.optFull, p.c0, 30);

    function [cSol] = CSTRsteady(u, c_in, k)
        % solves the CSTR at steady state
        v0 = sum(u);
        options_fsolve = optimset('Display','off');
        cSol = fsolve(@massBal,c_in,options_fsolve);
        
        function [F] = massBal(c)
            %evaluates the species balance
            r = reactantRate(c, k, s);
            F = u.*c_in/V+r-v0*c/V;
    	end
    end

    function [phi] = CSTRcost(u, c, c_in)
        % calculates the cost (-J)
        w = 0.004;
        J = c(3)^2*(u(1)+u(2))^2/(u(1)*c_in(1)) - w*(u(1)^2+u(2)^2);
        phi = -J;
    end

    function [G] = CSTRcon(c, k)
        % calculates the constraints (G) [2x1]
        G = zeros(2,1);
        
        % constraint 1: Q/Qmax - 1 < 0
        % Q = sum(r*dH*V)
        r = reactionRate(c, k, s);
        Q = sum(r.*dH.*V);
        G(1) = Q/Qmax-1;
        
        % constraint 1: D/Dmax - 1 < 0
        D = c(4)/sum(c);
        G(2) = D/Dmax - 1;
    end

    function [cost] = calcCost(x, steModel, c_in)
        % used by fmincon
        % x is the flowrates of A and B into the CSTR
        u = [x(1); x(2); 0; 0];
        c = steModel(u);
        cost = CSTRcost(u, c, c_in);
    end

    function [con, con_eq] = calcCon(x, steModel, k)
        % used by fmincon
        % x is the flowrates of A and B into the CSTR
        u = [x(1); x(2); 0; 0];
        c = steModel(u);
        con = CSTRcon(c, k);
        con_eq = [];
    end
    
    function [n] = nE(m,i)
        % outputs the ith element of matrix m
        n = m(i);
    end

    function [p] = convexApprox(fun, umin, umax, nP, deg, opt)
        % makes a convex approximation of "fun" between "umin" and "umax"
        % of degree "deg" around the point "opt"
        
        [u1Mesh, u2Mesh] = meshgrid(linspace(umin(1), umax(1), nP),...
                                    linspace(umin(2), umax(2), nP));
        funMesh = zeros(size(u1Mesh));
        
        for i = 1:numel(u1Mesh)
            funMesh(i) = fun([u1Mesh(i);u2Mesh(i);0;0]);
        end
        
        du1Mesh = u1Mesh - opt(1);
        du2Mesh = u2Mesh - opt(2);
        
        funOpt = fun([opt(1); opt(2); 0; 0]);
        
        p = fminsearch(@lsr,zeros(1,deg+1));
        
        
        function [S] = lsr(x)
            % funtion to find lsr for polyfit
            first = zeros(size(u1Mesh));
            second = zeros(size(u1Mesh));
            
            for j = 1:numel(u1Mesh)
                first(j) = [x(1), x(2)]*[du1Mesh(j);du2Mesh(j)];
                if deg == 2
                    second(j) = 0.5*[du1Mesh(j);du2Mesh(j)]'*(x(3)*eye(2))*...
                        [du1Mesh(j);du2Mesh(j)];
                else
                    second(j) = 0;
                end
            end
            
            funC = funOpt + first + second;
            residual = funMesh - funC;
            
            S = sum(sum(residual.^2));
        end
    end
        
    function [y] = convModel(p, x, x0, y0)
        dx = (x-x0);
        y = p(1:2)*dx + 0.5*p(3).*dx'*dx + y0;
    end

    function [dy] = convGrad(p, x, x0)
        dx = (x-x0);
        dy = p(1:2)' + p(3)*dx;
    end

    function [] = multiUnitTransRTO(convFunc, cGrad, plantFunc, gain, u0, c0, tau)
        % solve the transient RTO problem using the multiple units approach
        %
        % convFunc - convex approximation of the model
        % convGrad - gradient of the convex approximation of the model
        % plantFunc - real plant responce
        % gain - gain of between each step
        % u0 - inital condition
        % c0 - inital concentration
        % tau - time between each step
        %
        %
        
        du = 50e-3;
        %% first modification
        % first order modifier
        conv0 = convCell2vec(convFunc, u0(1:2));
        [c,plant0] = plantCell2vec(plantFunc, u0);
        
        e = plant0 - conv0;
        
        % initialise reactors
        base.C = c';
        base.J = plant0(1);
        base.G = permute(repmat(plant0(2:end),1,1,2),[1 3 2]);
        
        other.C = zeros(0, 2, 4);
        other.J = zeros(0, 2);
        other.G = zeros(0, 2, numel(plant0(2:end)));
        
        % get plant gradient
        [plantGradJ, plantGradG] = plantCellGrad(plantFunc, u0);
        
        % get modified grad
        mGrad = zeros(1, 2, 3);
        
        % second order modifier
        l = zeros(1, 2, 3);
        l(1, :, 1) = plantGradJ;
        l(1, :, 2:end) = plantGradG;
        
        modifiedFunc{1} = @(v)(convFunc{1}(v)+e(1,1)+(l(end,:,1))*(v-u0(1:2)));
        modifiedFunc{2} = @(v)(convFunc{2}(v)+e(1,2)+(l(end,:,2))*(v-u0(1:2)));
        modifiedFunc{3} = @(v)(convFunc{3}(v)+e(1,3)+(l(end,:,3))*(v-u0(1:2)));
        
        modifiedGrad{1} = @(v)(cGrad{1}(v)'+l(end,:,1));
        modifiedGrad{2} = @(v)(cGrad{2}(v)'+l(end,:,2));
        modifiedGrad{3} = @(v)(cGrad{3}(v)'+l(end,:,3));
        
        %run 
        unsolved = 1;
        i = 2;
        while unsolved
            % get new u
            newU = fmincon(modifiedFunc{1}, u0(1:2,i-1), [], [], [], [], [0, 0],...
                [50, 50], @(x)deal(convCell2vec({modifiedFunc{2:end}},x),[]), options_fmincon);
            
            u0(:,i) = zeros(4,1);
            u0(1:2,i) = newU;
            
            % first order modifier
            conv0 = convCell2vec(convFunc, u0(1:2,i));
            [c,plant0] = plantCell2vec(plantFunc, u0(:,i));
            
            e(i,:) = (1-gain).*e(i-1,:) + gain.*(plant0 - conv0);
            
            c0(:,i) = c;
            base.C(i,:) = c';
            base.J(i, 1) = plant0(1);
            base.G(i,:,:) = repmat(plant0(2:end),2,1);
            
            % get plant gradient
            [plantGradJ, plantGradG] = plantCellGrad(plantFunc, u0(:,i));
            
            % get model gradient
            mGrad(i,:,1) = cGrad{1}(u0(1:2,i));
            mGrad(i,:,2) = cGrad{2}(u0(1:2,i));
            mGrad(i,:,3) = cGrad{3}(u0(1:2,i));
            
            % second order modifier
            l(i, :, 1) = (1-gain).*l(i-1,:,1)+gain.*(plantGradJ - mGrad(i,:,1));
            l(i, :, 2:end) = (1-gain).*l(i-1,:,2:end)-gain.*(plantGradG - mGrad(i,:,2:end));
            
            modifiedFunc{1} = @(v)(convFunc{1}(v)+e(i,1)+(l(i,:,1))*(v-u0(1:2,i)));
            modifiedFunc{2} = @(v)(convFunc{2}(v)+e(i,2)+(l(i,:,2))*(v-u0(1:2,i)));
            modifiedFunc{3} = @(v)(convFunc{3}(v)+e(i,3)+(l(i,:,3))*(v-u0(1:2,i)));
            
            modifiedGrad{1} = @(v)(cGrad{1}(v)'+l(i,:,1));
            modifiedGrad{2} = @(v)(cGrad{2}(v)'+l(i,:,2));
            modifiedGrad{3} = @(v)(cGrad{3}(v)'+l(i,:,3));
            
            if i > 100
                unsolved = 0;
            else
                i=i+1;
            end
        end
        
        function [vec] = convCell2vec(funCell, x)
            % runs a cell of functions "funCell" with variable "x"
            vec = zeros(size(funCell));
            for j = 1:numel(funCell)
                vec(j) = funCell{j}(x);
            end
        end
        
        
        function [c, vec] = plantCell2vec(funCell, x, cNew)
            % runs a cell of functions "funCell" with variable "x"
            if nargin < 3
                cNew = c0(:,end);
            end
            vec = zeros(size(funCell(2:end)));
            c = funCell{1}(x, cNew, tau);
            
            for j = 2:numel(funCell)
                vec(j-1) = funCell{j}(x, cNew, tau);
            end
        end
        
        function [dJ, dG] = plantCellGrad(funCell, x)
            % calculates the plant gradients using the MU method
            % run other cases
            other.C(end+1,:,:) = zeros(1, size(other.C,2), size(other.C,3));
            other.J(end+1,:) = zeros(1, size(other.J,2));
            other.G(end+1,:,:) = zeros(1, size(other.G,2), size(other.G,3));
            for j = 1:2
                udu = x;
                udu(j) = udu(j)+du;
                if size(other.C,1) == 1
                    curC = c0;
                else
                    curC = permute(other.C(end-1, j, :),[3,2,1]);
                end
                [c,out] = plantCell2vec(funCell, udu, curC);
                other.C(end, j, :) = c;
                other.J(end, j) = out(1);
                other.G(end, j, :) = out(2:end);
            end
            
            dJ = (other.J(end,:)-base.J(end,:))./du;
            dG = (other.G(end,:,:)-base.G(end,:,:))./du;
            
        end
        
    end
end