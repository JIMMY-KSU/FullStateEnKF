%% Read in and Prep Results

 
load 'Data Files/Covariance.mat';


%% Plot Covariance 3D

x = 1:1:2701;
y = x;
[X,Y] = meshgrid(x,y);

Z = P_bar(7:2707,7:2707);

h = pcolor(Y,X,Z');
set(h,'EdgeColor','none')
%clabel(C, h)
colorbar
caxis([.5e-9,1e-9]); %force colorbar to be the same for each plot
ylabel('')
xlabel('')
zlabel('')
title(['Covariance'])