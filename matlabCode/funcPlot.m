function [] = funcPlot(func, umin, umax, nP)
% takes in a function of u1 and u2 and plots

[u1mesh,u2mesh] = meshgrid(linspace(umin(1),umax(1),nP),...
    linspace(umin(2),umax(2),nP));

f= zeros(nP);

for i = 1:numel(u1mesh)
    f(i) = func([u1mesh(i),u2mesh(i)]);
end

surf(f)

end