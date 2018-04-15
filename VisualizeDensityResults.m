
%% Read in and Prep Results

load 'Data Files/7Period_Results.mat'; 

% X_mean_updated_list_EnKF = mydict['X_mean_updated_list_EnKF'] 
% P_list_EnKF = mydict['P_list_EnKF'] 
% post_fit_list_EnKF = mydict['post_fit_list_EnKF'] 
% density_MSIS_array = mydict['density_MSIS_array'] 
% est_density_array = mydict['est_density_array'] 
% X_distribution = mydict['X_distribution'] 
% density_distribution = mydict['density_distribution'] 
% lat_lst_array = mydict['lat_lst_array'] 
% final_density_ensemble_est = mydict['final_density_ensemble_est']
% final_X_ensemble = mydict['final_X_ensemble']
% 'true_density_array': true_density_array*1e9, 
% 'final_density_grid_truth': final_density_grid_truth

%final filter ensemble of density grids estimate 
%average to get single grid
final_density_grid_est = mean(final_density_ensemble_est, 1);
final_density_grid_est = reshape(final_density_grid_est, 37, 73);

%final true density grid
final_density_grid_truth = reshape(final_density_grid_truth, 37, 73);


error_sum = sum(sum(Z))


%% Plot Final Time Results

x = -90:5:90;
y = 0:5:360;
[X,Y] = meshgrid(x,y);

%difference between true and estimated density grid
Z = final_density_grid_truth - final_density_grid_est;


h = pcolor(Y,X,Z');
set(h,'EdgeColor','none')
%clabel(C, h)
colorbar
caxis([-2e-4,4e-4]); %force colorbar to be the same for each plot
ylabel('Lat')
xlabel('LST')
zlabel('Density')
title(['Density Estimate Error, Error Integral = ', num2str(error_sum)])%', longitude_grid(ii),'


%% Process Time Series Results

error_integral_array = zeros(length(est_density_grid_array),1);
counter = 1;

for ii = 1:25:length(est_density_grid_array)
    
    %difference between true and estimated density grid
    Z = final_density_grid_truth_timeSeries(ii,:,:) - est_density_grid_array(ii,:,:);
    Z = reshape(Z, 37, 73);
    Z = abs(Z);

    error_integral_array(counter,1) = sum(sum(Z));  
    
    figure;
    h = pcolor(Y,X,Z');
    set(h,'EdgeColor','none')
    %clabel(C, h)
    colorbar
    caxis([0,4e-4]); %force colorbar to be the same for each plot
    ylabel('Lat')
    xlabel('LST')
    zlabel('Density')
    title(['Density Estimate Error, Error Integral = ', num2str(error_integral_array(ii,1))])
    
    counter = counter + 1;
    
end

figure;  
plot(error_integral_array)








