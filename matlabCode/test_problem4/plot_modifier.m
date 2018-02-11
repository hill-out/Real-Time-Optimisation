
% To generate plots for psfrag 
%axis([-70 70 0 50])
axis tight
set(gcf,'color','w')
set(gca,'fontsize',14)
%% 
myfig = gcf; 
grid off
set(myfig, 'WindowStyle', 'normal')
set(myfig, 'Position', [1200 500 450 300]) % [x(bottom left) y(bottom left) width height] for Windows!

%%
xlabel('x')

ylabel('y')

zlabel('z')
%%
title('')
%%
break

set(myleg,'FontSize',11);
%%
legend('d','e','f','g','h')
%% if plot won't fit in figure
myaxes = gca;
set(myaxes,'Position',[0.17 0 0.65 1]);
%%
%--- position legend first!
legend boxoff

%% 
legend( 'tracking nominal path', 'optimal')
%%
axis([0 130 800 1200])

%% 
axis fit