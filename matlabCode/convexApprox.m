function [convPara] = convexApprox(funArray, umin, umax, deg, opt, funOpt, reactOrder, kVal, cIn, V)
% convex approximation of a 2 variable input problem (later expand for n
% variable)

% number of points
nP1 = 31;
nP2 = 31;

% create meshes
[u1mesh, u2mesh] = meshgrid(linspace(umin(1),umax(1),nP1),linspace(umin(2),umax(2),nP2));
funmesh = zeros(nP1,nP2,length(funArray));
du1mesh = zeros(size(u1mesh));
du2mesh = zeros(size(u2mesh));
uBase = zeros(2,4); 

% solve CSTR and fun_val for all points
for i = 1:numel(u1mesh)  
    cSol = CSTR(reactOrder, kVal, cIn, [u1mesh(i),u2mesh(i), 0, 0], V);
    for j = 1:length(funArray)
        u = uBase;
        u(2,:) = cSol;
        u(1,[1,2]) = [u1mesh(i), u2mesh(i)];
        funmesh(i+numel(u1mesh)*(j-1)) = funArray{j}(u);
    end
end

% calc diff from opt point
du1mesh = u1mesh - opt(1);
du2mesh = u2mesh - opt(2);

% run solver for chosen degree
convPara = zeros(1,length(deg)*3);

for funN = 1:length(funArray)
    N = funN*3-2;
    convPara(N:N+deg(funN)) = fminsearch(@lsr,zeros(1,deg(funN)+1));
end

    function [S] = lsr(x)
        % funtion to find lsr for polyfit
        first = zeros(size(u1mesh));
        second = zeros(size(u1mesh));
        
        for i = 1:numel(du1mesh(:,1))
            for j = 1:numel(du1mesh(1,:))
                first(i,j) = [x(1), x(2)]*[du1mesh(i,j);du2mesh(i,j)];
                if deg(funN) == 2
                    second(i,j) = 0.5*[du1mesh(i,j);du2mesh(i,j)]'*(x(3)*eye(2))*...
                        [du1mesh(i,j);du2mesh(i,j)];
                else
                    second(i,j) = 0;
                end
            end
        end
        
        funC = funOpt(funN) + first + second;
        residual = funmesh(:,:,funN) - funC;
        
        S = sum(sum(residual.^2));
    end



end