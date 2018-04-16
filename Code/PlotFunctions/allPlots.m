function [fig] = allPlots(rpi, rp0, base, algo,tau, K, fig)
% Plots rr, phit, ut and gt
% -----------------------------------
% rpi       Iterative plant setpoints
% rp0       First setpoint from model
% base      Struct of plant time data
% algo      Algorithm for closed loop
% tau, K    Freq and gain of solution
% fig       Struct of figures to plot
% -----------------------------------
set(0,'DefaultFigureWindowStyle','docked')
dataName = sprintf('%s / \x03c4 = %4.0f / K = %5.2f',algo,tau,K);
try figure(fig.rr) % r1 vs. r2
    hold on
catch
    fig.rr = figure('name','r1 vs. r2');
    fig.h_rr = [];
end
fig.h_rr(end+1) = plot([rp0(1), rpi(:,1)'],[rp0(2), rpi(:,2)'],'-.','DisplayName',dataName);
legend('-DynamicLegend');
legend show

try figure(fig.phit) % phi vs t
    hold on
catch
    fig.phit = figure('name','phi vs t');
    fig.h_phit = [];
end
% make discontinuity at time step change
% method of discontinuity: 3 points with same time step, make middle point
% have a y-values of NaN.

i = any(base.t == base.ti,2);
i2 = (diff(i)==0 & i(1:end-1) == 1); % all places of double i
i3 = [1;find(i2);numel(i)];

base.t2 = [];
base.phip2 = [];
base.u2 = [];

for i4 = 1:numel(i3)-1
    base.t2 = [base.t2; NaN; base.t(i3(i4)+1:i3(i4+1))];
    base.phip2 = [base.phip2; NaN; base.phip(i3(i4)+1:i3(i4+1))];
    base.u2 = [base.u2; [NaN, NaN, NaN]; base.u(i3(i4)+1:i3(i4+1),:)];
    
    i3(i4) = i3(i4) + i4 - 1;
end

hold on
fig.h_phit(end+1) = plot(base.t,-base.phip,'DisplayName',dataName);
fig.h_phit(end+1) = plot(base.t2,-base.phip2,'DisplayName',dataName);
legend('-DynamicLegend');
legend show

try figure(fig.ut) % u vs t
catch
    fig.ut = figure('name','u vs t');
    fig.h_ut = [];
end

subplot(1,3,1)
hold on
fig.h_ut(end+1) = plot(base.t,base.u(:,1),'DisplayName',dataName);
fig.h_ut(end+1) = plot(base.t2,base.u2(:,1),'DisplayName',dataName);
legend('-DynamicLegend');
legend show

subplot(1,3,2)
hold on
fig.h_ut(end+1) = plot(base.t,base.u(:,2),'DisplayName',dataName);
fig.h_ut(end+1) = plot(base.t2,base.u2(:,2),'DisplayName',dataName);
legend('-DynamicLegend');
legend show

subplot(1,3,3)
hold on
fig.h_ut(end+1) = plot(base.t,base.u(:,3),'DisplayName',dataName);
fig.h_ut(end+1) = plot(base.t2,base.u2(:,3),'DisplayName',dataName);
legend('-DynamicLegend');
legend show

try figure(fig.gt) % g vs t
    hold on
catch
    fig.gt = figure('name','g vs t');
    fig.h_gt = [];
end
fig.h_gt = plot(base.t,base.g1p,'DisplayName',dataName);
legend('-DynamicLegend');
legend show

% try figure(fig.rrt) % r1 vs r2 vs t
%     hold on
% catch
%     fig.rrt = figure('name','r1 vs r2 vs t');
%     fig.h_rrt = [];
% end
% 
% t = base.t(i);
% t = t(2:2:numel(t));
% 
% fig.h_rrt = plot3([rp0(1), rpi(:,1)'],[rp0(2), rpi(:,2)'],t,'-o','DisplayName',dataName);
% legend('-DynamicLegend');
% legend show

end