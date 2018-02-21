function stairs2(x,y)
    for i=1:length(x)-1
        A(:,i) = [x(i) x(i+1)];
        B(:,i) = [y(i) y(i)];
    end
    plot(A,B,'LineWidth',5,'color',[0,0.45,0.67]);
    hold on
    stairs(x,y,'LineWidth',2,'color',[0,0.45,0.67],'linestyle',':')
    %[0,0.45,0.67]
    %[193/255,0,67/255]
end