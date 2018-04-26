%% Read in and Prep Results

 
load 'Data Files/Covariance.mat';


%% Plot Covariance 3D

x = 1:1:2701;
y = x;
[X,Y] = meshgrid(x,y);

Z = P_bar;

h = pcolor(Y,X,Z');
set(h,'EdgeColor','none')
%clabel(C, h)
colorbar
caxis([0,.3]); %force colorbar to be the same for each plot
ylabel('')
xlabel('')
zlabel('')
title(['Covariance'])