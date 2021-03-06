function [fig] = contPlot(X_As, F_Bs, fig)
% run the plant of choice until steady for a map of points to get contours

if nargin < 2
    F_Bs = linspace(3,12,41);
    if nargin < 1
        X_As = linspace(0.08,0.12,21);
    end
end

[a,b] = meshgrid(X_As, F_Bs);

rRange = [a(:),b(:)];

Kp = -1000;
T0 = 120;
xGuess = [0.09, 0.36, 0.1, 0.25, 0.1, 0.1];

for i = 1:size(a,1)
    for j = 1:size(a,2)
        r = [a(i,j),b(i,j)];
        up = plantController2(r,xGuess,Kp,T0)';
        [~,c] = ode15s(@(t,y)closedPlantODE(t,y,Kp), [0 10000],[up, xGuess]);
        up0 = c(end,1:3);
        Xp0 = c(end,4:end);
        
        phip(i,j) = phiFun(up0,Xp0);
        g1p(i,j) = g1Fun(up0,Xp0);
        g2p(i,j) = g2Fun(up0,Xp0);
        
        xGuess = Xp0;
    end
end

try figure(fig.rr) % r1 vs. r2
    hold on
catch
    fig.rr = figure('name','r1 vs. r2');
end

C = contour(X_As,F_Bs,-phip,[floor(min(min(-phip))/10)*10:10:ceil(max(max(-phip))/10)*10],'color','black');
run = 1;
i = 1;
j = 1;
ax = gca;
b1 = ax.DataAspectRatio;
b2 = ax.PlotBoxAspectRatio;
p = 9;
while run
    d = C(:,i);
    e = C(:,i+[1:d(2)]);
    diffChange = [1,find(sign(diff(e(2,1:end-1))) ~= sign(diff(e(2,2:end)))),size(e,2)];
    increasing = ~any(diffChange);
    if ~increasing
        % set e equal to the first range with p in it
        for loop = 1:numel(diffChange)-1
            if min(e(2,diffChange(loop):diffChange(loop+1)))<p && max(e(2,diffChange(loop):diffChange(loop+1)))>p
                e = e(:,diffChange(loop):diffChange(loop+1));
                break
            end
        end
    end
    
    if min(e(2,:))<p && max(e(2,:))>p
        
        f = spline(e(2,:),e(1,:),p);
        g = diff(e')';
        h = e(2,1:end-1) + g(2,:)/2;
        m(1) = spline(h,g(1,:),p);
        m(2) = spline(h,g(2,:),p);
        
        n = rad2deg(atan((b2(2)*m(2))/(b1(2)*m(1))));
        if n < 0
            n = n +180;
        end
        text(f,p,sprintf('%d',d(1)),'rotation',n,'horizontalalignment','center',...
            'verticalalignment','middle','backgroundcolor','w');
        
    end
    i = i+d(2)+1;
    if i > size(C,2)
        run = 0;
    end
    
end


figure('name','temp');
[a,~] = contourf(X_As,F_Bs,g1p,[0 0]);
close temp

figure(fig.rr)
c = patch([a(1,2:end),a(1,end)+0.01,a(1,end)+0.01],[a(2,2:end),a(2,end),a(2,2)],'red','facealpha',0.1,'linestyle','none');

end