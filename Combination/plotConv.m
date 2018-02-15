function plotConv(x,y,c,o)

f = convFun(c,(x-o)');
figure
hold on
scatter3(x(:,1),x(:,2),y)
scatter3(x(:,1),x(:,2),f)

end