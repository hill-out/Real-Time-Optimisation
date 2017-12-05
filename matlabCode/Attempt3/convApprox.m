function [cost, cons] = convApprox(stru)
% -------------------------------------------------------------------------
% gets a convex approximation of the cost and constraints of the model
%
% stru      - struct    - structure of parameters
% 
% cost      - fun_hand  - function of cost (2nd order)
% cons      - fun_hand  - function of constraints (1st order)
%
% -------------------------------------------------------------------------

% mesh u
np = 8;
[u1mesh, u2mesh] = meshgrid(linspace(1,50,np));
costmesh = zeros(size(u1mesh));
consmesh = zeros([size(u1mesh),2]);

% run for all u's
for i = 1:numel(u1mesh)
    [~,costmesh(i),consmesh(mod(i-1,np)+1,floor((i-1)/np)+1,:)] =...
        CSTRste([u1mesh(i),u2mesh(i)],stru);
end

% run lsr
du1mesh = u1mesh - stru.uOpt(1);
du2mesh = u2mesh - stru.uOpt(2);

% for cost
costPara = fminsearch(@(x)lsr(x, 2, costmesh, stru.costOpt),[0,0,0]);
costPara = [-0.8305, -0.9121, 0.08];
cost = @(x, ~)(convModel(costPara, x, stru.uOpt, stru.costOpt));

% for cons
cons1Para = zeros(1,3);
cons2Para = zeros(1,3);

cons1Para(1:2) = fminsearch(@(x)lsr(x, 1, consmesh(:,:,1), stru.consOpt(1)),[0,0]);
cons2Para(1:2) = fminsearch(@(x)lsr(x, 1, consmesh(:,:,2), stru.consOpt(2)),[0,0]);


cons = @(x, ~)([convModel(cons1Para, x, stru.uOpt, stru.consOpt(1)),...
    convModel(cons2Para, x, stru.uOpt, stru.consOpt(2))]);

    function [S] = lsr(x,deg, funmesh, funOpt)
        % funtion to find lsr for polyfit
        first = zeros(size(u1mesh));
        second = zeros(size(u1mesh));
        
        for i = 1:np
            for j = 1:np
                first(i,j) = [x(1), x(2)]*[du1mesh(i,j);du2mesh(i,j)];
                if deg == 2
                    second(i,j) = 0.5*[du1mesh(i,j);du2mesh(i,j)]'*(x(3)*eye(2))*...
                        [du1mesh(i,j);du2mesh(i,j)];
                else
                    second(i,j) = 0;
                end
            end
        end
        
        funC = funOpt + first + second;
        residual = funmesh - funC;
        
        S = sum(sum(residual.^2));
    end



end