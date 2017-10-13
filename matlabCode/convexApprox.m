function [] = convexApprox(fun, umin, umax, deg, opt, reactOrder, kVal, c0, V)

nP1 = 21;
nP2 = 21;

[u1mesh, u2mesh] = meshgrid(linspace(umin(1),umax(1),nP1),linspace(umin(2),umax(2),nP2));
funmesh = zeros(size(u1mesh));
du1mesh = zeros(size(u1mesh));
du2mesh = zeros(size(u2mesh));

for i = 1:numel(u1mesh)  
    cSol = CSTR(reactOrder, kVal, c0, [u1mesh(i),u2mesh(i), 0, 0], V);
    funmesh(i) = fun(c0,cSol,[u1mesh(i), u2mesh(i)]);
end

du1mesh = u1mesh - opt(1);
du2mesh = u2mesh - opt(2);

cOpt = CSTR(reactOrder, kVal, c0, [opt(1), opt(2), 0, 0], V);
funOpt = fun(c0,cOpt,[opt(1), opt(2), 0, 0]);

convPara = fminsearch(@lsr,[0,0,0])

    function [S] = lsr(x)
        
        first = zeros(size(u1mesh));
        second = zeros(size(u1mesh));
        
        for i = 1:numel(du1mesh(:,1))
            for j = 1:numel(du1mesh(1,:))
                first(i,j) = [x(1), x(2)]*[du1mesh(i,j);du2mesh(i,j)];
                second(i,j) = 0.5*[du1mesh(i,j);du2mesh(i,j)]'*(x(3)*eye(2))*...
                    [du1mesh(i,j);du2mesh(i,j)];
            end
        end
        
        funC = funOpt + first + second;
        residual = funmesh - funC;
        
        S = sum(sum(residual.^2));
    end



end