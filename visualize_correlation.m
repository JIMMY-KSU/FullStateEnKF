
load 'Data Files/ensemble_cor_smooth.mat'; %size: 2701, 2701



% x = 1:1:2701;
% y = 1:1:2701;
% [X,Y] = meshgrid(x,y);
% 
% Z = cor_mat; 
% 
% figure;
% mesh(X,Y,Z')
% xlabel('X')
% ylabel('Y')
% zlabel('correlations')
% title('Correlations')



cor_mat_reshaped = reshape(cor_mat', 2701, 37, 73);
cor_mat_reshaped = reshape(cor_mat_reshaped, 37, 73, 37, 73);

%correlations with lat/lst = 0
cor_mat_oneRow = cor_mat_reshaped(19,1,:,:);

x = -90:5:90;
y = 0:5:360;
[X,Y] = meshgrid(x,y);

Z = reshape(cor_mat_oneRow, 37, 73); 

figure;
%contourf(X,Y,Z')
pcolor(X,Y,Z')
xlabel('latitude')
ylabel('longitude')
zlabel('correlations')
title('Correlations w/ 0 lat/LST')



% %correlations with lat = 0, lst = 180
% cor_mat_oneRow = cor_mat_reshaped(19,37,:,:);
% 
% x = -90:5:90;
% y = 0:5:360;
% [X,Y] = meshgrid(x,y);
% 
% Z = reshape(cor_mat_oneRow, 37, 73); 
% 
% figure;
% mesh(X,Y,Z')
% xlabel('latitude')
% ylabel('longitude')
% zlabel('correlations')
% title('Correlations w/ 0 lat, 180 LST')
% 
% 
% Z = reshape(density(1,:,:), 37, 73);
% 
% figure;
% mesh(X,Y,Z')
% xlabel('latitude')
% ylabel('longitude')
% zlabel('density')
% title('Density Member')


%%
%investigate different kinds of contour plots as function of LST
%close all;

x = rad2deg(latitude_grid);
y = rad2deg(longitude_grid);
[X,Y] = meshgrid(x,y);

counter = 1;
figure
lat_index = 28;
%loop through LSTs at the equator
for ii = 1:5:73
    
    cor_mat_oneRow = cor_mat_reshaped(lat_index,ii,:,:); %19=0lat, 28=45lat, 37=90lat
    Z = reshape(cor_mat_oneRow, 37, 73); 
    
    %mean normalize correlation values (Z)
    Z = (Z-mean(mean(Z)))/(max(max(Z))-min(min(Z)));

    
    subplot(3,5,counter)       % add first plot in 2 x 2 grid
    %pcolor(Y,X,Z')
    h = pcolor(Y,X,Z');
    set(h,'EdgeColor','none')
    %clabel(C, h)
    colorbar
    caxis([-.5,.5]); %force colorbar to be the same for each plot
    ylabel('Lat')
    xlabel('LST')
    zlabel('correlations')
    title([num2str(x(lat_index)),' Lat, ', num2str(y(ii)), ' LST'])%', longitude_grid(ii),'

    counter = counter + 1;
    
end


%%
%investigate different kinds of contour plots as function of latitude
%close all;

x = rad2deg(latitude_grid);
y = rad2deg(longitude_grid);
[X,Y] = meshgrid(x,y);

counter = 1;
figure
lst_index = 37;
%loop through LSTs at the equator
for ii = sort(1:3:37, 'descend')
    
    cor_mat_oneRow = cor_mat_reshaped(ii,lst_index,:,:); %37=180lst, 19=90lst
    Z = reshape(cor_mat_oneRow, 37, 73);
    
    %mean normalize correlation values (Z)
    Z = (Z-mean(mean(Z)))/(max(max(Z))-min(min(Z)));

    
    subplot(3,5,counter)       % add first plot in 2 x 2 grid
    pcolor(Y,X,Z')
    %clabel(C, h)
    colorbar
    caxis([-.5,.5]); %force colorbar to be the same for each plot
    ylabel('Lat')
    xlabel('LST')
    zlabel('correlations')
    title([num2str(x(ii)),' Lat, ', num2str(y(lst_index)), ' LST'])%', longitude_grid(ii),'

    counter = counter + 1;
    
end



