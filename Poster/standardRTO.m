% generate standard RTO poster graphics

syms a b c

b = a+10./a - 4;
c = 0.3*a.^2-3*a+9.5;

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
plot(a,eval(b),'b','LineWidth',3)
plot(a,eval(c),'r','LineWidth',3)

a = c_opt;
mC = eval(c);
plot([c_opt c_opt],[0 mC],'r--')

a = c_opt;
m0 = eval(b);
m1 = eval(db) - eval(dc);

clear a
syms a
bGrad = m0 + m1*(a-c_opt);
a = 1:0.1:9;
plot(a,eval(bGrad),'k--')

ylim([1,9])
xlim([0,10])

set(gca, 'Xtick', [], 'Ytick', [], 'box', 'off')

xO = 0.3;  
yO = 0.15;
patch(...
    [10-xO -yO; 10-xO +yO; 10 0.0], ...
    [1+yO 9-xO; 1-yO 9-xO; 1 9], 'k', 'clipping', 'off')

patch(...
    [c_opt, c_opt; c_opt+0.1, c_opt+0.1; c_opt-0.1, c_opt-0.1], ...
    [m0, mC; m0-0.2 mC+0.2; m0-0.2 mC+0.2], 'k', 'clipping', 'off', 'EdgeColor', 'k')

plot([c_opt c_opt],[m0 mC], 'k')

text(    6.5, 0.65, 'Constraints, u', 'fontsize', 18)
text(  -0.35,  6.5, 'Cost, \phi', 'rotation', 90, 'fontsize', 18)
text( c_opt-0.1, 0.6, 'u_0^{\ast}', 'fontsize', 18)
a = 9.1;
text( 9.2, eval(b), '\phi_p', 'fontsize', 18)
text( 9.2, eval(c), '\phi_0', 'fontsize', 18)
text( 9.2, eval(bGrad), '\lambda_0', 'fontsize', 18)
text( c_opt+0.2, (m0+mC)/2+0.1, '\epsilon_0', 'fontsize', 14)

figure
hold on

a = c_opt;
m0 = eval(b) - eval(c);
m1 = eval(db) - eval(dc);

clear a
syms a
mod = c + m0 + m1*(a-c_opt);
a = 1:0.1:9;
plot(a,eval(b),'b','LineWidth',3)
plot(a,eval(mod),'r','LineWidth',3)

dmod = diff(mod);
mod_opt = eval(solve(dmod == 0));
mod_opt = mod_opt(mod_opt>0);

ylim([1,9])
xlim([0,10])

set(gca, 'Xtick', [], 'Ytick', [], 'box', 'off')

xO = 0.3;  
yO = 0.15;
patch(...
    [10-xO -yO; 10-xO +yO; 10 0.0], ...
    [1+yO 9-xO; 1-yO 9-xO; 1 9], 'k', 'clipping', 'off')


a = mod_opt;
plot([mod_opt mod_opt],[0 eval(mod)],'r--')

text(    6.5, 0.65, 'Constraints, u', 'fontsize', 18)
text(  -0.35,  6.5, 'Cost, \phi', 'rotation', 90, 'fontsize', 18)
text( mod_opt-0.1, 0.6, 'u_1^{\ast}', 'fontsize', 18)
a = 9.1;
text( 9.2, eval(b), '\phi_p', 'fontsize', 18)
text( 8.5, 8, '\phi_1', 'fontsize', 18)
