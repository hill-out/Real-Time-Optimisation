function plotConv(x,y,c)

f = convFun(c,x');
figure
hold on
scatter3(x(:,1),x(:,2),y)
scatter3(x(:,1),x(:,2),f)

end