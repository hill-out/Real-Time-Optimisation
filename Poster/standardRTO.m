% generate standard RTO poster graphics

%% unmodded graph
syms a b c

b = 1.2*a+10./a - 5;
%c = 0.3*a.^2-3*a+9.5;
c = 0.25*((a-5).^2)-0.3*25+9;

%optimum points
db = diff(b);
b_opt = eval(solve(db == 0));
b_opt = b_opt(b_opt>0);

dc = diff(c);
c_opt = eval(solve(dc == 0));
c_opt = c_opt(c_opt>0);

figure
hold on
a = 1:0.1:9;
plot(a,eval(b),'color',[0,50/255,95/255],'LineWidth',5)
plot(a,eval(c),'color',[193/255,0,67/255],'LineWidth',5)

a = c_opt;
mC = eval(c);
plot([c_opt c_opt],[0 mC],'--','color',[193/255,0,67/255],'LineWidth',2)

a = c_opt;
m0 = eval(b);
m1 = eval(db) - eval(dc);

clear a
syms a
bGrad = m0 + m1*(a-c_opt);
a = 1:0.1:9;
plot(a,eval(bGrad),'k--','LineWidth',2)

ylim([1,7.9])
xlim([0,9.9])

set(gca, 'Xtick', [], 'Ytick', [], 'box', 'off', 'LineWidth', 2)

xO = 0.2;  
yO = 0.1;
patch(...
    [10-xO -yO; 10-xO +yO; 10 0.0], ...
    [1+yO 8-xO; 1-yO 8-xO; 1 8], 'k', 'clipping', 'off')

patch(...
    [c_opt, c_opt; c_opt+0.1, c_opt+0.1; c_opt-0.1, c_opt-0.1], ...
    [m0, mC; m0-0.2 mC+0.2; m0-0.2 mC+0.2], 'k', 'clipping', 'off', 'EdgeColor', 'k')

plot([c_opt c_opt],[m0 mC], 'k', 'LineWidth', 2)

annotation('arrow',[0.67 0.80],[0.73 0.69],'Color',[0.50 0.50 0.50],...
    'LineWidth',2);

text(    7, 0.65, 'Conditions (u)', 'fontsize', 32)
text(  -0.4,  6, 'Cost (\phi)', 'rotation', 90, 'fontsize', 32)
text( c_opt-0.1, 0.6, 'u_0^{\ast}', 'fontsize', 32)
a = 9;
text( 9.2, eval(b), '\phi_p', 'fontsize', 40)
text( 9.2, eval(c), '\phi_0', 'fontsize', 40)
a = 3;
text( 3.2, 6.3, '\lambda_0 = \nabla_{u}\phi_p({u^{\ast}_0})', 'fontsize', 32)
text( c_opt+0.2, (m0+mC)/2+0.1, '$\varepsilon_0$', 'Interpreter', 'latex', 'fontsize', 48)

%% modified graph
figure
hold on

a = c_opt;
m0 = eval(b) - eval(c);
m1 = eval(db) - eval(dc);

clear a
syms a
mod = c + m0 + m1*(a-c_opt);
a = 1:0.1:9;
plot(a,eval(b),'color',[0,50/255,95/255],'LineWidth',5)
plot(a,eval(mod),'color',[193/255,0,67/255],'LineWidth',5)

dmod = diff(mod);
mod_opt = eval(solve(dmod == 0));
mod_opt = mod_opt(mod_opt>0);

ylim([1,7.9])
xlim([0,9.9])

set(gca, 'Xtick', [], 'Ytick', [], 'box', 'off', 'LineWidth', 2)

xO = 0.2;  
yO = 0.1;
patch(...
    [10-xO -yO; 10-xO +yO; 10 0.0], ...
    [1+yO 8-xO; 1-yO 8-xO; 1 8], 'k', 'clipping', 'off')


a = mod_opt;
plot([mod_opt mod_opt],[0 eval(mod)],'--','color',[193/255,0,67/255],'LineWidth',2)
a = c_opt;
plot([c_opt c_opt],[0 eval(mod)],'k--','LineWidth',2)

text(    7, 0.65, 'Conditions (u)', 'fontsize', 32)
text(  -0.4,  6, 'Cost (\phi)', 'rotation', 90, 'fontsize', 32)
text( c_opt-0.1, 0.6, 'u_0^{\ast}', 'fontsize', 32)
text( mod_opt-0.1, 0.6, 'u_1^{\ast}', 'fontsize', 32)
a = 9;
text( 9.2, eval(b), '\phi_p', 'fontsize', 40)
text( 6.5, eval(b), '\phi_1', 'fontsize', 40)

%% full optimisation
unsolved = 1;
u_opt = [c_opt, mod_opt];
while unsolved
    a = mod_opt;
    m0 = eval(b) - eval(c);
    m1 = eval(db) - eval(dc);
    clear a
    syms a
    mod = c + m0 + m1*(a-mod_opt);
    dmod = diff(mod);
    mod_opt = eval(solve(dmod == 0));
    mod_opt = mod_opt(mod_opt>0);
    
    u_opt(end+1) = mod_opt;
    if abs(mod_opt-b_opt) <1e-3
        unsolved = 0;
    end
end

figure
hold on
iter = 0:(numel(u_opt)-1);
plot([0,(numel(u_opt)-1)],[b_opt,b_opt],':','color',[0,50/255,95/255],'LineWidth',5)
plot(iter,u_opt,'Color',[193/255,0,67/255],'LineWidth',5)

set(gca, 'Xtick', [0 1 2 3 4 5 6 7 8 9], 'Ytick', [], 'box', 'off', 'LineWidth', 2, ...
    'FontSize', 16)

ya = 2;
yb = 6;
xa = 0;
xb = 10;
ylim([ya,yb-0.1])
xlim([xa,xb-0.1])

xO = 0.2;  
yO = 0.1*(yb-ya)/7;
x1 = 0.2*(yb-ya)/7;
y1 = 0.1;
patch(...
    [10-xO 0-y1; 10-xO 0+y1; 10 0.0], ...
    [2+yO 6-x1; 2-yO 6-x1; 2 6], 'k', 'clipping', 'off')

text(  xb/2-1.2, 1.6, 'Iteration (k)', 'fontsize', 32)
text(  -0.4, 3.5, 'Optimum Condition', 'rotation', 90, 'fontsize', 30)
text(  -0.65, b_opt, 'u^{\ast}_p', 'fontsize', 32)
text(  1, 4, 'u^{\ast}_k', 'fontsize', 32)


