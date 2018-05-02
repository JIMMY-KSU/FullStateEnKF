
%% Read in and Prep Results

 
load 'Data Files/RR_EnSRF_q39_75obs.mat';

%load 'Data Files/7Period_Results_180RAAN.mat';


% X_mean_updated_list_EnKF = mydict['X_mean_updated_list_EnKF'] 
% P_list_EnKF = mydict['P_list_EnKF'] 
% post_fit_list_EnKF = mydict['post_fit_list_EnKF'] 
% density_MSIS_array = mydict['density_MSIS_array'] 
% est_density_array = mydict['est_density_array'] 
% X_distribution = mydict['X_distribution'] 
% density_distribution = mydict['density_distribution'] 
% lat_lst_array = mydict['lat_lst_array'] 
% final_density_ensem+ble_est = mydict['final_density_ensemble_est']
% final_X_ensemble = mydict['final_X_ensemble']
% 'true_density_array': true_density_array, 
% 'final_density_grid_truth_timeSeries': final_density_grid_truth_timeSeries

%final filter ensemble of density grids estimate 
%average to get single grid
final_density_grid_est = mean(final_density_ensemble_est, 1);
final_density_grid_est = reshape(final_density_grid_est, 37, 73);

stop_index = length(final_density_grid_truth_timeSeries);

%final true density grid
final_density_grid_truth = reshape(final_density_grid_truth_timeSeries(stop_index,:,:), 37, 73);


%indices by station
indices_1 = measurement_array(1:stop_index,2) == 1;
indices_2 = measurement_array(1:stop_index,2) == 2;
indices_3 = measurement_array(1:stop_index,2) == 3;
indices_4 = measurement_array(1:stop_index,2) == 4;

%% percent



perc_error = 100 * abs(est_density_array - true_density_array(1:stop_index))./true_density_array(1:stop_index);

figure;
plot(measurement_array(1:stop_index,1), perc_error)




%% Plot Final Time Results

x = -90:5:90;
y = 0:5:360;
[X,Y] = meshgrid(x,y);
% 
% %difference between true and estimated density grid
% Z = final_density_grid_truth - final_density_grid_est;
% Z = Z.^2;
% error_sum = sum(sum(Z));
% 
% h = pcolor(Y,X,Z');
% set(h,'EdgeColor','none')
% %clabel(C, h)
% colorbar
% caxis([0,1e-6]); %force colorbar to be the same for each plot
% ylabel('Lat')
% xlabel('LST')
% zlabel('Density')
% title(['Density Estimate Error, Error Integral = ', num2str(error_sum)])%', longitude_grid(ii),'


%% Process Time Series Results


counter = 1;
figure;
title(['Density Estimate Error'])
    
for ii = 1:11:length(est_density_grid_array)
    
    ii
    
    %difference between true and estimated density grid
    Z = final_density_grid_truth_timeSeries(ii,:,:) - est_density_grid_array(ii,:,:);
    Z = reshape(Z, 37, 73);
    Z = Z.^2;


    subplot(1,7,counter)
    h = pcolor(Y,X,Z');
    set(h,'EdgeColor','none')
    %clabel(C, h)
    colorbar
    caxis([0,1e-7]); %force colorbar to be the same for each plot
    ylabel('Lat')
    xlabel('LST')
    zlabel('Density')
    title(['Time Step: ', num2str(ii)])
    
    counter = counter + 1;
    
end

%figure;  
%plot(error_integral_array(2:end))



%%
error_integral_array = 1;% = zeros(length(est_density_grid_array),1);
error_integral_array_sqrd = 1;

for ii = 1:1:length(est_density_grid_array)
    
    %difference between true and estimated density grid
    Z = final_density_grid_truth_timeSeries(ii,:,:) - est_density_grid_array(ii,:,:);
    
    error_integral_array = [error_integral_array; sum(sum(Z))];
    
    Z = Z.^2;
    error_integral_array_sqrd = [error_integral_array_sqrd; sum(sum(Z))];
    
    
end

% figure;  
% plot(error_integral_array(2:end))
% ylabel('Error Integral')
% xlabel('Time')
% title(['Density Error Integral Over Grid Time Series (not abs or sqrd)'])

error_integral_array_sqrd = error_integral_array_sqrd(2:end);

figure;  
hold;
scatter(measurement_array(indices_1,1), error_integral_array_sqrd(indices_1),'filled', 's', 'MarkerFaceColor','m')
scatter(measurement_array(indices_2,1), error_integral_array_sqrd(indices_2),'filled', '^', 'MarkerFaceColor',[0 0.5 0])
scatter(measurement_array(indices_3,1), error_integral_array_sqrd(indices_3),'filled', 'd', 'MarkerFaceColor','r')
scatter(measurement_array(indices_4,1), error_integral_array_sqrd(indices_4),'filled', 'o', 'MarkerFaceColor','k')
ylabel('Error Integral')
xlabel('Time')
%legend({'Station 1','Station 2', 'Station 3', 'Station 4'})
title(['Density Error Integral Over Grid Time Series '])



%% Covariance Plots

counter = 1;
figure;

x = 1:1:7;
y = 1:1:7;
[X,Y] = meshgrid(x,y);


for ii = 1:10:length(est_density_grid_array)
    
    Z = P_list_EnKF(ii,:,:);
    Z = reshape(Z, 7, 7);
    
    
    subplot(1,7,counter)
    h = pcolor(Y,X,Z');
    set(h,'EdgeColor','none')
    %clabel(C, h)
    colorbar
    caxis([0,1e-7]); %force colorbar to be the same for each plot
    ylabel('State: XYZ Position, XYZ Velocity, Density')
    xlabel('State')
    zlabel('State Covariance')
    title(['Time Step: ', num2str(ii)])
    
    counter = counter + 1;
    
end




