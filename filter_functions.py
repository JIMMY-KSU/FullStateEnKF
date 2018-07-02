import numpy as np
import math
import scipy
import operator
import copy
import matplotlib.pyplot as plt

#MSIS: https://github.com/DeepHorizons/Python-NRLMSISE-00
#import time
from nrlmsise_00_header import *
from nrlmsise_00 import *
#SUBROUTINE GTD7D -- d[5] is the "effective total mass density
#for drag" and is the sum of the mass densities of all species
#in this model, INCLUDING anomalous oxygen.




def calc_MSIS_density(t, X_vector, day_of_year_init, day_of_month_init, hour_init, month_init, year_init, omega_const, r_earth_const):
    
    state = X_vector[0:3] 
    
    
    (latitude, longitude, altitude, day_of_month, hour) = calc_lat_lon_from_t_R(t, state, day_of_month_init, hour_init, \
                                                                           month_init, year_init, omega_const, r_earth_const)
    
    day_of_year = math.floor(t/86400) + day_of_year_init
    t = t - math.floor(t/86400) * 86400
    
    
    if longitude < 0:
            longitude = longitude + 2*math.pi
    lst = calc_LST(hour, longitude) #lst in units of hours
    
    #convert to radians
    lst_rad = (lst/24) * np.radians(360) #radians, not used, for debugging only
    
    lon = 0 #dependent only on lst


    Output = nrlmsise_output()
    Input = nrlmsise_input()
    flags = nrlmsise_flags()
    aph = ap_array()

    for i in range(7):
        aph.a[i]=100
    flags.switches[0] = 1
    for i in range(1, 24):
        flags.switches[i]=1

    Input.doy = day_of_year
    Input.year = 0 #/* without effect */
    Input.sec = t
    Input.alt = altitude #km
    Input.g_lat = math.degrees(latitude)
    Input.g_long = math.degrees(lon)
    Input.lst = lst
    Input.f107A = 180 #I believe this is a "nominal" value
    Input.f107 = 180
    Input.ap = 10 

    gtd7d(Input, flags, Output)

    density = Output.d[5] #total mass density (grams/m^3, m^3 b/c switches[0] = 1)
    
    return density


#altitude in km
def gen_one_ensemble(latitudes, longitudes, alt, day_of_year, t):

    num_lat = len(latitudes)
    num_lon = len(longitudes)
    density_grid = np.zeros((num_lat, num_lon))
    
    day_of_year = math.floor(t/86400) + day_of_year
    t = t - math.floor(t/86400) * 86400
    
    
    #loop through spatial latitude/longitude grid

    for ii in range(num_lat):

        lat = latitudes[ii]

        for jj in range(num_lon):

            lon = longitudes[jj]
            
            #check if lon is negative or wrapped around after 2 pi
            if lon > 2*math.pi:
                lon = lon - 2*math.pi
            elif lon < 0:
                lon = lon + 2*math.pi

            lst = lon/(2*math.pi) * 24 #fraction of 360 degrees converted to hours

            lon = 0

            Output = nrlmsise_output()
            Input = nrlmsise_input()
            flags = nrlmsise_flags()
            aph = ap_array()

            for i in range(7):
                aph.a[i]=100
            flags.switches[0] = 1
            for i in range(1, 24):
                flags.switches[i]=1

            Input.doy = day_of_year
            Input.year = 0 #/* without effect */
            Input.sec = t
            Input.alt = alt
            Input.g_lat = math.degrees(lat)
            Input.g_long = math.degrees(lon)
            Input.lst = lst #I believe this needs to be lst b/c will populate for all local sidereal times
            #as though greenwich is midnight, but this shell is not actually dependent on location on Earth,
            #it is only dependent on its angle/location WRT the sun, will convert this grid to sun fixed 
            #coordinate system post generation
            Input.f107A = 180 #I believe this is a "nominal" value
            Input.f107 = 180
            Input.ap = 10 

            gtd7d(Input, flags, Output)

            density = Output.d[5] #total mass density (grams/m^3, m^3 b/c switches[0] = 1)
            density = density 
            #print(rho)

            density_grid[ii, jj] = density
            
            
    return density_grid



#t_UT : current universal time
#current longitude of space object
#returns LST in hours
def calc_LST(t_UT, longitude):
    
    #t_UT = t_UT%24 this was causing an issue in filter (greater error) not sure why..
    
    #calculate change in longitude from UT and convert to change in hours/time
    #note: negative lon/lon above 180 is a negative difference in time & opposite for less than 180
    delta_lon = 360 - np.degrees(longitude)
    delta_lst = (delta_lon/360) * 24 #24 hours, unit = hours
    
    t_lst = t_UT - delta_lst
    if t_lst < 0:
        t_lst = t_lst + 24
    
    #return now if want in hours

    
    return t_lst


#depending on the local time of the orbit (its pos & time) calculate the corresponding lat/lon
#indices of the ensemble in order to extract corresponding density
#t : time since beginning of orbit simulation (seconds)
#year, month, day, hour_init : the UT time at the beginning of the simulation
#R : current position of orbital object, ECI
def calc_lat_lst_indices(t, R, day_of_month_init, hour_init, month_init, year_init, omega_const, r_earth_const,\
                         lat_res, lon_res):
    
    (latitude, longitude, altitude, day_of_month, hour) = calc_lat_lon_from_t_R(t, R, \
                                    day_of_month_init, hour_init, month_init, year_init, omega_const, r_earth_const)
    

    if longitude < 0:
        longitude = longitude + 2*math.pi
    lst = calc_LST(hour, longitude) #lst in units of hours
    
    if lst > 24:
        lst = 24
    
    #convert to radians
    lst = (lst/24) * np.radians(360) #radians


    #in ensembles: longitude = 0 -> LST = 12, lon = 0 -> LST = 0
    #determine closest ensemble lat/lon to the actual lat/lon of the orbital object ------------------------------

    #calculate grid tick for latitude
    remainder = latitude%lat_res
    lat_grid_ticks = (latitude - remainder)/lat_res
    #check if the grid tick should be that one or the next one, if remainder > half of the resolution,
    #should be next grid tick
    if remainder > lat_res/2:
        lat_grid_ticks = lat_grid_ticks + 1
    lat_grid_ticks = lat_grid_ticks + math.floor(np.radians(90)/lat_res) #grid tick = 0 -> lat = -90 deg

    #calculate grid tick for latitude
    remainder = lst%lon_res
    lst_grid_ticks = (lst - remainder)/lon_res
    #check if the grid tick should be that one or the next one, if remainder > half of the resolution,
    #should be next grid tick
    if remainder > lon_res/2:
        lst = lst + 1
    
    lat_grid_ticks = int(lat_grid_ticks)
    lst_grid_ticks = int(lst_grid_ticks)
        

    return (lat_grid_ticks, lst_grid_ticks)


#calculate the latitude and longitude of an object based UCT time (t) and object ECI position (R)
def calc_lat_lon_from_t_R(t, R, day_of_month_init, hour_init, month_init, year_init, omega_const, r_earth_const):
    
    day_of_month = math.floor(t/86400) + day_of_month_init
    t = t - math.floor(t/86400) * 86400
        
    #Calculations for theta_gmst (rotation btwn ECI & ECEF for the current time of simulation/orbit)--------------
    hour = hour_init + t/(60*60) #hours (float) since midnight UT
    jd = calc_julian_date(year_init, month_init, day_of_month, hour) #date of interest in UT   
    
    T_UT = (jd - 2451545)/36525 #calc T_UT at this delta_t/date 
    #calculate theta gmst @ 0 hr using the T_UT and the eq. from Vallado 
    theta_gmst_0_hr = math.radians(100.4606184 + 36000.77005361*T_UT + \
                                           .00038793*T_UT**2 - 2.6e-8*T_UT**3)
    seconds = hour * 60 * 60 #seconds since midnight UT 
    #calculate theta gmst using theta gmst @ 0 hr and seconds 
    #since midnight*the rotation rate of earth
    theta_gmst = theta_gmst_0_hr + omega_const * seconds 
    theta_gmst = theta_gmst % (2*math.pi)
    
    #calculate position in ECEF & then the lat/lon/alt-----------------------------------------------------------
    r_ecef = eci2ecef(R, theta_gmst)
    (latitude, longitude, altitude) = ecef2geo_lat_lon_alt(r_ecef, r_earth_const)
    
    
    return (latitude, longitude, altitude, day_of_month, hour)


#input the year, month, day (of month, not of year), and hour in order to calculate the 
#corresponding julian date
def calc_julian_date(year, month, day, hour):
    julian_day = day - 32075 + int(1461 * (year + 4800 + int((month - 14) / 12)) / 4) + \
            int(367 * (month - 2 - int((month - 14) / 12) * 12) / 12) - \
            int(3 * (int((year + 4900 + int((month - 14) / 12)) / 100)) / 4)
    julian_day = julian_day - 0.5 + hour / 24.0
    return julian_day

def RIC_2_ECI(pos, vel, Q_ric):
    
    u_R = pos/np.linalg.norm(pos)
    u_N = np.cross(pos, vel)/np.linalg.norm(np.cross(pos, vel))
    u_T = np.cross(u_N, u_R)
    
    gamma = np.array([u_R, u_T, u_N]).reshape(3,3)
    
    Q_eci = np.dot(gamma.T,np.dot(Q_ric, gamma))
    return Q_eci



def calc_display_results_pre(pre_fit_list, prefit_bounds_list, measurement_array, R, meas_type, stop_index, saveFig_bool, time_str):
    
    
    rms_1 = 'Range ='
    unit_1 = 'km'
    ylabel_1 = 'Range Residuals (km)'
    title_1 = 'Range pre-fit Residuals'
    save_fig_1 = 'prefit_range.png'
    rms_2 = 'Range Rate ='
    unit_2 = 'km/s'
    ylabel_2 = 'Range Rate Residuals (km/s)'
    title_2 = 'Range Rate pre-fit Residuals'
    save_fig_2 = 'prefit_rangeRate.png'
    pre_fit_list_new = copy.deepcopy(pre_fit_list)
         
    if (meas_type == 2) or (meas_type == 3):
        rms_1 = 'Azimuth ='
        unit_1 = 'degrees'
        rms_2 = 'Elevation ='
        unit_2 = 'degrees'
        pre_fit_list_new[:, 0] = np.degrees(pre_fit_list[:, 0])
        pre_fit_list_new[:, 1] = np.degrees(pre_fit_list[:, 1])
        ylabel_1 = 'Azimuth Residuals (degrees)'
        title_1 = 'Azimuth pre-fit Residuals'
        save_fig_1 = 'prefit_az.png'
        ylabel_2 = 'Elevation Residuals (degrees)'
        title_2 = 'Elevation pre-fit Residuals'
        save_fig_2 = 'prefit_el_rate.png'

        
    

    times = measurement_array[:stop_index,0]/(60)
    
    indices_1 = np.where(measurement_array[:stop_index, 1] == 1)[0]
    indices_2 = np.where(measurement_array[:stop_index, 1] == 2)[0]
    indices_3 = np.where(measurement_array[:stop_index, 1] == 3)[0]
    indices_4 = np.where(measurement_array[:stop_index, 1] == 4)[0]

    
    
    #pre-fit
    print('pre-fit RMS:')
    pre_fit_1_list_4RMS = pre_fit_list_new[:, 0]
    prefit_1_rms = np.sqrt(np.mean(np.square(pre_fit_1_list_4RMS)))
    print(rms_1, "%.4f" % prefit_1_rms, unit_1)

    pre_fit_2_list_4RMS = pre_fit_list_new[:, 1]
    prefit_2_rms = np.sqrt(np.mean(np.square(pre_fit_2_list_4RMS)))
    print(rms_2, "%.4f" % prefit_2_rms, unit_2)
    
    if meas_type == 3:
        pre_fit_3_list_4RMS = pre_fit_list_new[:, 2]
        prefit_3_rms = np.sqrt(np.mean(np.square(pre_fit_3_list_4RMS)))
        print('Range =', "%.3f" % prefit_3_rms, 'km')
    
    
    covar_env_upper1 = np.degrees(abs(prefit_bounds_list[:stop_index, 0]))*3
    covar_env_upper2 = np.degrees(abs(prefit_bounds_list[:stop_index, 1]))*3

   
    #pre-fit Residuals
    fig_preFit_az = plt.figure()
    plt.plot(times, covar_env_upper1, label='_nolegend_', c='g')
    plt.plot(times, -covar_env_upper1, label='_nolegend_', c='g')
    plt.scatter(times[indices_1], pre_fit_list_new[indices_1, 0], s=50, c='m', marker='s')
    plt.scatter(times[indices_2], pre_fit_list_new[indices_2, 0], s=50, c='g', marker='^')
    plt.scatter(times[indices_3], pre_fit_list_new[indices_3, 0], s=50, c='r', marker='D')
    plt.scatter(times[indices_4], pre_fit_list_new[indices_4, 0], s=50, c='k', marker='o')
    plt.ylabel(ylabel_1, fontsize=18)
    plt.xlabel(time_str, fontsize=18)
    plt.title(title_1, fontsize=18)
    legend_names = ['Station 1', 'Station 2']
    plt.legend(legend_names, fontsize=10)
    plt.xlim([times[0] - 5, times[-1] + 5])
    plt.tight_layout()
    plt.show()
    #fig.savefig(save_fig_1)

    fig_preFit_el = plt.figure()
    plt.plot(times, covar_env_upper2, label='_nolegend_', c='g')
    plt.plot(times, -covar_env_upper2, label='_nolegend_', c='g')
    plt.scatter(times[indices_1], pre_fit_list_new[indices_1, 1], s=50, c='m', marker='s')
    plt.scatter(times[indices_2], pre_fit_list_new[indices_2, 1], s=50, c='g', marker='^')
    plt.scatter(times[indices_3], pre_fit_list_new[indices_3, 1], s=50, c='r', marker='D')
    plt.scatter(times[indices_4], pre_fit_list_new[indices_4, 1], s=50, c='k', marker='o')
    plt.ylabel(ylabel_2, fontsize=18)
    plt.xlabel(time_str, fontsize=18)
    plt.title(title_2, fontsize=18)
    plt.legend(legend_names, fontsize=10)
    #plt.ylim([-.01,.01])
    plt.xlim([times[0] - 5, times[-1] + 5])
    plt.tight_layout()
    plt.show()
    #fig.savefig(save_fig_2)
    
    if meas_type == 3:
        
        covar_env_upper3 = abs(prefit_bounds_list[:stop_index, 2])*3
        
        fig_preFit_range = plt.figure()
        plt.plot(times, covar_env_upper3, label='_nolegend_', c='g')
        plt.plot(times, -covar_env_upper3, label='_nolegend_', c='g')
        plt.scatter(times[indices_1], pre_fit_list_new[indices_1, 2], s=50, c='m', marker='s')
        plt.scatter(times[indices_2], pre_fit_list_new[indices_2, 2], s=50, c='g', marker='^')
        plt.scatter(times[indices_3], pre_fit_list_new[indices_3, 2], s=50, c='r', marker='D')
        plt.scatter(times[indices_4], pre_fit_list_new[indices_4, 2], s=50, c='k', marker='o')
        plt.ylabel('Range Residuals (km)', fontsize=18)
        plt.xlabel(time_str, fontsize=18)
        plt.title('Range pre-fit Residuals', fontsize=18)
        plt.legend(legend_names, fontsize=10)
        #plt.ylim([-.01,.01])
        plt.xlim([times[0] - 5, times[-1] + 5])
        plt.tight_layout()
        plt.show()
        #fig.savefig('prefit_range.png')
    
    if saveFig_bool:
        fig_preFit_az.savefig('Figures/preFit_az.png')
        fig_preFit_el.savefig('Figures/preFit_el.png')
        fig_preFit_range.savefig('Figures/preFit_range.png')
        
        
def calc_display_results(post_fit_list, measurement_array, R, meas_type, stop_index, saveFig_bool, time_str):
    
    
    rms_1 = 'Range ='
    unit_1 = 'km'
    ylabel_1 = 'Range Residuals (km)'
    title_1 = 'Range Post-fit Residuals'
    save_fig_1 = 'postfit_range.png'
    rms_2 = 'Range Rate ='
    unit_2 = 'km/s'
    ylabel_2 = 'Range Rate Residuals (km/s)'
    title_2 = 'Range Rate Post-fit Residuals'
    save_fig_2 = 'postfit_rangeRate.png'
    post_fit_list_new = copy.deepcopy(post_fit_list)
         
    if (meas_type == 2) or (meas_type == 3):
        rms_1 = 'Azimuth ='
        unit_1 = 'degrees'
        rms_2 = 'Elevation ='
        unit_2 = 'degrees'
        post_fit_list_new[:, 0] = np.degrees(post_fit_list[:, 0])
        post_fit_list_new[:, 1] = np.degrees(post_fit_list[:, 1])
        ylabel_1 = 'Azimuth Residuals (degrees)'
        title_1 = 'Azimuth Post-fit Residuals'
        save_fig_1 = 'postfit_az.png'
        ylabel_2 = 'Elevation Residuals (degrees)'
        title_2 = 'Elevation Post-fit Residuals'
        save_fig_2 = 'postfit_el_rate.png'

        
    

    times = measurement_array[:stop_index,0]/(60)
    
    indices_1 = np.where(measurement_array[:stop_index, 1] == 1)[0]
    indices_2 = np.where(measurement_array[:stop_index, 1] == 2)[0]
    indices_3 = np.where(measurement_array[:stop_index, 1] == 3)[0]
    indices_4 = np.where(measurement_array[:stop_index, 1] == 4)[0]

    
    
    #Post-fit
    print('Post-fit RMS:')
    post_fit_1_list_4RMS = post_fit_list_new[:, 0]
    postfit_1_rms = np.sqrt(np.mean(np.square(post_fit_1_list_4RMS)))
    print(rms_1, "%.4f" % postfit_1_rms, unit_1)

    post_fit_2_list_4RMS = post_fit_list_new[:, 1]
    postfit_2_rms = np.sqrt(np.mean(np.square(post_fit_2_list_4RMS)))
    print(rms_2, "%.4f" % postfit_2_rms, unit_2)
    
    if meas_type == 3:
        post_fit_3_list_4RMS = post_fit_list_new[:, 2]
        postfit_3_rms = np.sqrt(np.mean(np.square(post_fit_3_list_4RMS)))
        print('Range =', "%.5f" % postfit_3_rms, 'km')
    
    
    covar_env_upper1 = np.ones((stop_index)) * np.degrees(np.sqrt(abs(R[0, 0])))*3
    covar_env_upper2 = np.ones((stop_index)) * np.degrees(np.sqrt(abs(R[1, 1])))*3

   
    #Post-fit Residuals
    fig_postFit_az = plt.figure()
    plt.plot(times, covar_env_upper1, label='_nolegend_', c='g')
    plt.plot(times, -covar_env_upper1, label='_nolegend_', c='g')
    plt.scatter(times[indices_1], post_fit_list_new[indices_1, 0], s=50, c='m', marker='s')
    plt.scatter(times[indices_2], post_fit_list_new[indices_2, 0], s=50, c='g', marker='^')
    plt.scatter(times[indices_3], post_fit_list_new[indices_3, 0], s=50, c='r', marker='D')
    plt.scatter(times[indices_4], post_fit_list_new[indices_4, 0], s=50, c='k', marker='o')
    plt.ylabel(ylabel_1, fontsize=18)
    plt.xlabel(time_str, fontsize=18)
    plt.title(title_1, fontsize=18)
    legend_names = ['Station 1', 'Station 2']
    plt.legend(legend_names, fontsize=10)
    plt.xlim([times[0] - 5, times[-1] + 5])
    plt.tight_layout()
    plt.show()
    #fig.savefig(save_fig_1)

    fig_postFit_el = plt.figure()
    plt.plot(times, covar_env_upper2, label='_nolegend_', c='g')
    plt.plot(times, -covar_env_upper2, label='_nolegend_', c='g')
    plt.scatter(times[indices_1], post_fit_list_new[indices_1, 1], s=50, c='m', marker='s')
    plt.scatter(times[indices_2], post_fit_list_new[indices_2, 1], s=50, c='g', marker='^')
    plt.scatter(times[indices_3], post_fit_list_new[indices_3, 1], s=50, c='r', marker='D')
    plt.scatter(times[indices_4], post_fit_list_new[indices_4, 1], s=50, c='k', marker='o')
    plt.ylabel(ylabel_2, fontsize=18)
    plt.xlabel(time_str, fontsize=18)
    plt.title(title_2, fontsize=18)
    plt.legend(legend_names, fontsize=10)
    #plt.ylim([-.01,.01])
    plt.xlim([times[0] - 5, times[-1] + 5])
    plt.tight_layout()
    plt.show()
    #fig.savefig(save_fig_2)
    
    if meas_type == 3:
        
        covar_env_upper3 = np.ones((stop_index)) * np.sqrt(abs(R[2, 2]))*3
        
        fig_postFit_range = plt.figure()
        plt.plot(times, covar_env_upper3, label='_nolegend_', c='g')
        plt.plot(times, -covar_env_upper3, label='_nolegend_', c='g')
        plt.scatter(times[indices_1], post_fit_list_new[indices_1, 2], s=50, c='m', marker='s')
        plt.scatter(times[indices_2], post_fit_list_new[indices_2, 2], s=50, c='g', marker='^')
        plt.scatter(times[indices_3], post_fit_list_new[indices_3, 2], s=50, c='r', marker='D')
        plt.scatter(times[indices_4], post_fit_list_new[indices_4, 2], s=50, c='k', marker='o')
        plt.ylabel('Range Residuals (km)', fontsize=18)
        plt.xlabel(time_str, fontsize=18)
        plt.title('Range Post-fit Residuals', fontsize=18)
        legend_names = ['Station 1', 'Station 2']
        plt.legend(legend_names, fontsize=10)
        #plt.ylim([-.01,.01])
        plt.xlim([times[0] - 5, times[-1] + 5])
        plt.tight_layout()
        plt.show()
        #fig.savefig('postfit_range.png')
    
    if saveFig_bool:
        fig_postFit_az.savefig('Figures/postFit_az.png')
        fig_postFit_el.savefig('Figures/postFit_el.png')
        fig_postFit_range.savefig('Figures/postFit_range.png')
        
    
    
def plot_error_covar_xref(P_list, x_ref_updated_list, obs_data_truth, density_truth, x_range, y_range, \
                          z_range, xv_range, yv_range, zv_range, measurement_array, time, stop_index, \
                          saveFig_bool, time_str, title_str):
    
    #Compare to the Truth Data : Estimation Errors------

    times = time[:stop_index]/(60)
    
    indices_1 = np.where(measurement_array[:stop_index, 1] == 1)[0]
    indices_2 = np.where(measurement_array[:stop_index, 1] == 2)[0]
    indices_3 = np.where(measurement_array[:stop_index, 1] == 3)[0]
    indices_4 = np.where(measurement_array[:stop_index, 1] == 4)[0]
    
    density_covar_env_upper = np.sqrt(abs(P_list[:stop_index, -1, -1]))*3
    density_covar_env_lower = -density_covar_env_upper
    density_error = x_ref_updated_list[:stop_index,-1] - density_truth[:stop_index]
    
    x_covar_env_upper = np.sqrt(abs(P_list[:stop_index, 0, 0]))*3
    x_covar_env_lower = -x_covar_env_upper
    x_error = x_ref_updated_list[:stop_index,0] - obs_data_truth[:stop_index, 0]
    
    y_covar_env_upper = np.sqrt(abs(P_list[:stop_index, 1, 1]))*3
    y_covar_env_lower = -y_covar_env_upper
    y_error = x_ref_updated_list[:stop_index,1] - obs_data_truth[:stop_index, 1]
    
    z_covar_env_upper = np.sqrt(abs(P_list[:stop_index, 2, 2]))*3
    z_covar_env_lower = -z_covar_env_upper
    z_error = x_ref_updated_list[:stop_index,2] - obs_data_truth[:stop_index, 2]
    
    #error_pos_norm = np.sqrt(x_error**2 + y_error**2 + z_error**2)
    error_density_rms_3D = np.sqrt(np.mean(np.square(density_error)))
    print('Density RMS =', "%.5f" % error_density_rms_3D, r'$(kg/km^3)$')
    
    print('Position RMS:')
    error_x_pos_rms_3D = np.sqrt(np.mean(np.square(x_error)))
    print('X =', "%.4f" % error_x_pos_rms_3D, 'km')
    
    error_y_pos_rms_3D = np.sqrt(np.mean(np.square(y_error)))
    print('Y =', "%.4f" % error_y_pos_rms_3D, 'km')
    
    error_z_pos_rms_3D = np.sqrt(np.mean(np.square(z_error)))
    print('Z =', "%.4f" % error_z_pos_rms_3D, 'km')
    
    pos_rms = np.sqrt(error_x_pos_rms_3D**2 + error_y_pos_rms_3D**2 + error_z_pos_rms_3D**2)
    print('Overall =', "%.4f" % pos_rms, 'km')
    
    
    #Density
    fig_dens = plt.figure()
    plt.plot(times, density_covar_env_upper, label='_nolegend_', c='g')
    plt.plot(times, density_covar_env_lower, label='_nolegend_', c='g')
    plt.scatter(times[indices_1], density_error[indices_1], s=50, c='m', marker='s')
    plt.scatter(times[indices_2], density_error[indices_2], s=50, c='g', marker='^')
    plt.scatter(times[indices_3], density_error[indices_3], s=50, c='r', marker='D')
    plt.scatter(times[indices_4], density_error[indices_4], s=50, c='k', marker='o')
    plt.ylabel(r'$(kg/km^3)$', fontsize=18)
    plt.xlabel(time_str, fontsize=18)
    legend_names = ['Station 1', 'Station 2']
    plt.legend(legend_names, fontsize=10)
    plt.title(title_str + ' Density Error & Covariance Envelope', fontsize=18)
    #plt.ylim([-x_range,x_range])
    plt.xlim([times[0] - 5, times[-1] + 5])
    plt.show()
    

    #x Position
    fig_xpos = plt.figure()
    plt.plot(times, x_covar_env_upper, label='_nolegend_', c='g')
    plt.plot(times, x_covar_env_lower, label='_nolegend_', c='g')
    
    plt.scatter(times[indices_1], x_error[indices_1], s=50, c='m', marker='s')
    plt.scatter(times[indices_2], x_error[indices_2], s=50, c='g', marker='^')
    plt.scatter(times[indices_3], x_error[indices_3], s=50, c='r', marker='D')
    plt.scatter(times[indices_4], x_error[indices_4], s=50, c='k', marker='o')
    plt.ylabel('km', fontsize=18)
    plt.xlabel(time_str, fontsize=18)
    plt.legend(legend_names, fontsize=10)
    plt.title(title_str + ' X Position Error & Covariance Envelope', fontsize=18)
    plt.ylim([-x_range,x_range])
    plt.xlim([times[0] - 5, times[-1] + 5])
    plt.show()

    #y Position 
    fig_ypos = plt.figure()
    plt.plot(times, y_covar_env_upper, label='_nolegend_', c='g')
    plt.plot(times, y_covar_env_lower, label='_nolegend_', c='g')
    plt.scatter(times[indices_1], y_error[indices_1], s=50, c='m', marker='s')
    plt.scatter(times[indices_2], y_error[indices_2], s=50, c='g', marker='^')
    plt.scatter(times[indices_3], y_error[indices_3], s=50, c='r', marker='D')
    plt.scatter(times[indices_4], y_error[indices_4], s=50, c='k', marker='o')
    plt.ylabel('km', fontsize=18)
    plt.xlabel(time_str, fontsize=18)
    plt.legend(legend_names, fontsize=10)
    plt.title(title_str + ' Y Position Error & Covariance Envelope', fontsize=18)
    plt.ylim([-y_range,y_range])
    plt.xlim([times[0] - 5, times[-1] + 5])
    plt.show()
    #fig.savefig('y_pos_error.png')

    #z Position
    fig_zpos = plt.figure()
    plt.plot(times, z_covar_env_upper, label='_nolegend_', c='g')
    plt.plot(times, z_covar_env_lower, label='_nolegend_', c='g')
    plt.scatter(times[indices_1], z_error[indices_1], s=50, c='m', marker='s')
    plt.scatter(times[indices_2], z_error[indices_2], s=50, c='g', marker='^')
    plt.scatter(times[indices_3], z_error[indices_3], s=50, c='r', marker='D')
    plt.scatter(times[indices_4], z_error[indices_4], s=50, c='k', marker='o')
    plt.ylabel('km', fontsize=18)
    plt.xlabel(time_str, fontsize=18)
    plt.legend(legend_names, fontsize=10)
    plt.title(title_str + ' Z Position Error & Covariance Envelope', fontsize=18)
    plt.ylim([-z_range,z_range])
    plt.xlim([times[0] - 5, times[-1] + 5])
    plt.show()
    
    #x Velocity
    x_dot_covar_env_upper = np.sqrt(abs(P_list[:stop_index, 3, 3]))*3
    x_dot_covar_env_lower = -np.sqrt(abs(P_list[:stop_index, 3, 3]))*3
    x_vel_error = x_ref_updated_list[:stop_index,3] - obs_data_truth[:stop_index, 3]

    fig_xvel = plt.figure()
    plt.plot(times, x_dot_covar_env_upper, label='_nolegend_', c='g')
    plt.plot(times, x_dot_covar_env_lower, label='_nolegend_', c='g')
    plt.scatter(times[indices_1], x_vel_error[indices_1], s=50, c='m', marker='s')
    plt.scatter(times[indices_2], x_vel_error[indices_2], s=50, c='g', marker='^')
    plt.scatter(times[indices_3], x_vel_error[indices_3], s=50, c='r', marker='D')
    plt.scatter(times[indices_4], x_vel_error[indices_4], s=50, c='k', marker='o')
    plt.ylabel('km/second', fontsize=18)
    plt.xlabel(time_str, fontsize=18)
    plt.legend(legend_names, fontsize=10)
    plt.title(title_str + ' X Velocity Error & Covariance Envelope', fontsize=18)
    plt.ylim([-xv_range,xv_range])
    plt.xlim([times[0] - 5, times[-1] + 5])
    plt.show()
    #fig.savefig('x_vel_error.png')

    #y Velocity
    y_dot_covar_env_upper = np.sqrt(abs(P_list[:stop_index, 4, 4]))*3
    y_dot_covar_env_lower = -np.sqrt(abs(P_list[:stop_index, 4, 4]))*3
    y_vel_error = x_ref_updated_list[:stop_index,4] - obs_data_truth[:stop_index, 4]

    fig_yvel = plt.figure()
    plt.plot(times, y_dot_covar_env_upper, label='_nolegend_', c='g')
    plt.plot(times, y_dot_covar_env_lower, label='_nolegend_', c='g')
    plt.scatter(times[indices_1], y_vel_error[indices_1], s=50, c='m', marker='s')
    plt.scatter(times[indices_2], y_vel_error[indices_2], s=50, c='g', marker='^')
    plt.scatter(times[indices_3], y_vel_error[indices_3], s=50, c='r', marker='D')
    plt.scatter(times[indices_4], y_vel_error[indices_4], s=50, c='k', marker='o')
    plt.ylabel('km/second', fontsize=18)
    plt.xlabel(time_str, fontsize=18)
    plt.legend(legend_names, fontsize=10)
    plt.title(title_str + ' Y Velocity Error & Covariance Envelope', fontsize=18)
    plt.ylim([-yv_range,yv_range])
    plt.xlim([times[0] - 5, times[-1] + 5])
    plt.show()

    #z Velocity
    z_dot_covar_env_upper = np.sqrt(abs(P_list[:stop_index, 5, 5]))*3
    z_dot_covar_env_lower = -np.sqrt(abs(P_list[:stop_index, 5, 5]))*3
    z_vel_error = x_ref_updated_list[:stop_index,5] - obs_data_truth[:stop_index, 5]

    fig_zvel = plt.figure()
    plt.plot(times, z_dot_covar_env_upper, label='_nolegend_', c='g')
    plt.plot(times, z_dot_covar_env_lower, label='_nolegend_', c='g') 
    plt.scatter(times[indices_1], z_vel_error[indices_1], s=50, c='m', marker='s')
    plt.scatter(times[indices_2], z_vel_error[indices_2], s=50, c='g', marker='^')
    plt.scatter(times[indices_3], z_vel_error[indices_3], s=50, c='r', marker='D')
    plt.scatter(times[indices_4], z_vel_error[indices_4], s=50, c='k', marker='o')
    plt.ylabel('km/second', fontsize=18)
    plt.xlabel(time_str, fontsize=18)
    plt.legend(legend_names, fontsize=10)
    plt.title(title_str + ' Z Velocity Error & Covariance Envelope', fontsize=18)
    plt.ylim([-zv_range,zv_range])
    plt.xlim([times[0] - 5, times[-1] + 5])
    plt.show()
    
    
    print('Velocity RMS:')
    error_x_vel_rms_3D = np.sqrt(np.mean(np.square(x_vel_error)))
    print('X =', "%.6f" % error_x_vel_rms_3D, 'km/second')
    
    error_y_vel_rms_3D = np.sqrt(np.mean(np.square(y_vel_error)))
    print('Y =', "%.6f" % error_y_vel_rms_3D, 'km/second')
    
    error_z_vel_rms_3D = np.sqrt(np.mean(np.square(z_vel_error)))
    print('Z =', "%.6f" % error_z_vel_rms_3D, 'km/second')
    
    vel_rms = np.sqrt(error_x_vel_rms_3D**2 + error_y_vel_rms_3D**2 + error_z_vel_rms_3D**2)
    print('Overall =', "%.6f" % vel_rms, 'km/s')

    
    if saveFig_bool:
        fig_dens.savefig('Figures/dens_error.png')
        fig_xpos.savefig('Figures/x_pos_error.png')
        fig_ypos.savefig('Figures/y_pos_error.png')
        fig_zpos.savefig('Figures/z_pos_error.png')
        fig_xvel.savefig('Figures/x_vel_error.png')
        fig_yvel.savefig('Figures/y_vel_error.png')
        fig_zvel.savefig('Figures/z_vel_error.png')
    
        
    

def plot_error_covar_xref_noDensity(P_list, x_ref_updated_list, obs_data_truth, x_range, y_range, \
                          z_range, xv_range, yv_range, zv_range, measurement_array, time, stop_index, saveFig_bool, time_str):
    
    #Compare to the Truth Data : Estimation Errors------

    times = time[:stop_index]/(60)
    
    indices_1 = np.where(measurement_array[:stop_index, 1] == 1)[0]
    indices_2 = np.where(measurement_array[:stop_index, 1] == 2)[0]
    indices_3 = np.where(measurement_array[:stop_index, 1] == 3)[0]
    indices_4 = np.where(measurement_array[:stop_index, 1] == 4)[0]
    
    
    x_covar_env_upper = np.sqrt(abs(P_list[:stop_index, 0, 0]))*3
    x_covar_env_lower = -x_covar_env_upper
    x_error = x_ref_updated_list[:stop_index,0] - obs_data_truth[:stop_index, 0]
    
    y_covar_env_upper = np.sqrt(abs(P_list[:stop_index, 1, 1]))*3
    y_covar_env_lower = -y_covar_env_upper
    y_error = x_ref_updated_list[:stop_index,1] - obs_data_truth[:stop_index, 1]
    
    z_covar_env_upper = np.sqrt(abs(P_list[:stop_index, 2, 2]))*3
    z_covar_env_lower = -z_covar_env_upper
    z_error = x_ref_updated_list[:stop_index,2] - obs_data_truth[:stop_index, 2]
    
    
    print('Position RMS:')
    error_x_pos_rms_3D = np.sqrt(np.mean(np.square(x_error)))
    print('X =', "%.4f" % error_x_pos_rms_3D, 'km')
    
    error_y_pos_rms_3D = np.sqrt(np.mean(np.square(y_error)))
    print('Y =', "%.4f" % error_y_pos_rms_3D, 'km')
    
    error_z_pos_rms_3D = np.sqrt(np.mean(np.square(z_error)))
    print('Z =', "%.4f" % error_z_pos_rms_3D, 'km')
    
    pos_rms = np.sqrt(error_x_pos_rms_3D**2 + error_y_pos_rms_3D**2 + error_z_pos_rms_3D**2)
    print('Overall =', "%.4f" % pos_rms, 'km')
    
   
    

    #x Position
    fig_xpos = plt.figure()
    plt.plot(times, x_covar_env_upper, label='_nolegend_', c='g')
    plt.plot(times, x_covar_env_lower, label='_nolegend_', c='g')
    
    plt.scatter(times[indices_1], x_error[indices_1], s=50, c='m', marker='s')
    plt.scatter(times[indices_2], x_error[indices_2], s=50, c='g', marker='^')
    plt.scatter(times[indices_3], x_error[indices_3], s=50, c='r', marker='D')
    plt.scatter(times[indices_4], x_error[indices_4], s=50, c='k', marker='o')
    plt.ylabel('km', fontsize=18)
    plt.xlabel(time_str, fontsize=18)
    legend_names = ['Station 1', 'Station 2', 'Station 3', 'Station 4']
    plt.legend(legend_names, fontsize=10)
    plt.title('EnKF X Position Error & Covariance Envelope', fontsize=18)
    plt.ylim([-x_range,x_range])
    plt.xlim([0,times[-1]])
    plt.show()

    #y Position 
    fig_ypos = plt.figure()
    plt.plot(times, y_covar_env_upper, label='_nolegend_', c='g')
    plt.plot(times, y_covar_env_lower, label='_nolegend_', c='g')
    plt.scatter(times[indices_1], y_error[indices_1], s=50, c='m', marker='s')
    plt.scatter(times[indices_2], y_error[indices_2], s=50, c='g', marker='^')
    plt.scatter(times[indices_3], y_error[indices_3], s=50, c='r', marker='D')
    plt.scatter(times[indices_4], y_error[indices_4], s=50, c='k', marker='o')
    plt.ylabel('km', fontsize=18)
    plt.xlabel(time_str, fontsize=18)
    plt.legend(legend_names, fontsize=10)
    plt.title('EnKF Y Position Error & Covariance Envelope', fontsize=18)
    plt.ylim([-y_range,y_range])
    plt.xlim([0,times[-1]])
    plt.show()
    #fig.savefig('y_pos_error.png')

    #z Position
    fig_zpos = plt.figure()
    plt.plot(times, z_covar_env_upper, label='_nolegend_', c='g')
    plt.plot(times, z_covar_env_lower, label='_nolegend_', c='g')
    plt.scatter(times[indices_1], z_error[indices_1], s=50, c='m', marker='s')
    plt.scatter(times[indices_2], z_error[indices_2], s=50, c='g', marker='^')
    plt.scatter(times[indices_3], z_error[indices_3], s=50, c='r', marker='D')
    plt.scatter(times[indices_4], z_error[indices_4], s=50, c='k', marker='o')
    plt.ylabel('km', fontsize=18)
    plt.xlabel(time_str, fontsize=18)
    plt.legend(legend_names, fontsize=10)
    plt.title('EnKF Z Position Error & Covariance Envelope', fontsize=18)
    plt.ylim([-z_range,z_range])
    plt.xlim([0,times[-1]])
    plt.show()
    
    #x Velocity
    x_dot_covar_env_upper = np.sqrt(abs(P_list[:stop_index, 3, 3]))*3
    x_dot_covar_env_lower = -np.sqrt(abs(P_list[:stop_index, 3, 3]))*3
    x_vel_error = x_ref_updated_list[:stop_index,3] - obs_data_truth[:stop_index, 3]

    fig_xvel = plt.figure()
    plt.plot(times, x_dot_covar_env_upper, label='_nolegend_', c='g')
    plt.plot(times, x_dot_covar_env_lower, label='_nolegend_', c='g')
    plt.scatter(times[indices_1], x_vel_error[indices_1], s=50, c='m', marker='s')
    plt.scatter(times[indices_2], x_vel_error[indices_2], s=50, c='g', marker='^')
    plt.scatter(times[indices_3], x_vel_error[indices_3], s=50, c='r', marker='D')
    plt.scatter(times[indices_4], x_vel_error[indices_4], s=50, c='k', marker='o')
    plt.ylabel('km/second', fontsize=18)
    plt.xlabel(time_str, fontsize=18)
    plt.legend(legend_names, fontsize=10)
    plt.title('EnKF X Velocity Error & Covariance Envelope', fontsize=18)
    plt.ylim([-xv_range,xv_range])
    plt.xlim([0,times[-1]])
    plt.show()
    #fig.savefig('x_vel_error.png')

    #y Velocity
    y_dot_covar_env_upper = np.sqrt(abs(P_list[:stop_index, 4, 4]))*3
    y_dot_covar_env_lower = -np.sqrt(abs(P_list[:stop_index, 4, 4]))*3
    y_vel_error = x_ref_updated_list[:stop_index,4] - obs_data_truth[:stop_index, 4]

    fig_yvel = plt.figure()
    plt.plot(times, y_dot_covar_env_upper, label='_nolegend_', c='g')
    plt.plot(times, y_dot_covar_env_lower, label='_nolegend_', c='g')
    plt.scatter(times[indices_1], y_vel_error[indices_1], s=50, c='m', marker='s')
    plt.scatter(times[indices_2], y_vel_error[indices_2], s=50, c='g', marker='^')
    plt.scatter(times[indices_3], y_vel_error[indices_3], s=50, c='r', marker='D')
    plt.scatter(times[indices_4], y_vel_error[indices_4], s=50, c='k', marker='o')
    plt.ylabel('km/second', fontsize=18)
    plt.xlabel(time_str, fontsize=18)
    plt.legend(legend_names, fontsize=10)
    plt.title('EnKF Y Velocity Error & Covariance Envelope', fontsize=18)
    plt.ylim([-yv_range,yv_range])
    plt.xlim([0,times[-1]])
    plt.show()

    #z Velocity
    z_dot_covar_env_upper = np.sqrt(abs(P_list[:stop_index, 5, 5]))*3
    z_dot_covar_env_lower = -np.sqrt(abs(P_list[:stop_index, 5, 5]))*3
    z_vel_error = x_ref_updated_list[:stop_index,5] - obs_data_truth[:stop_index, 5]

    fig_zvel = plt.figure()
    plt.plot(times, z_dot_covar_env_upper, label='_nolegend_', c='g')
    plt.plot(times, z_dot_covar_env_lower, label='_nolegend_', c='g') 
    plt.scatter(times[indices_1], z_vel_error[indices_1], s=50, c='m', marker='s')
    plt.scatter(times[indices_2], z_vel_error[indices_2], s=50, c='g', marker='^')
    plt.scatter(times[indices_3], z_vel_error[indices_3], s=50, c='r', marker='D')
    plt.scatter(times[indices_4], z_vel_error[indices_4], s=50, c='k', marker='o')
    plt.ylabel('km/second', fontsize=18)
    plt.xlabel(time_str, fontsize=18)
    plt.legend(legend_names, fontsize=10)
    plt.title('EnKF Z Velocity Error & Covariance Envelope', fontsize=18)
    plt.ylim([-zv_range,zv_range])
    plt.xlim([0,times[-1]])
    plt.show()
    
    
    print('Velocity RMS:')
    error_x_vel_rms_3D = np.sqrt(np.mean(np.square(x_vel_error)))
    print('X =', "%.6f" % error_x_vel_rms_3D, 'km/second')
    
    error_y_vel_rms_3D = np.sqrt(np.mean(np.square(y_vel_error)))
    print('Y =', "%.6f" % error_y_vel_rms_3D, 'km/second')
    
    error_z_vel_rms_3D = np.sqrt(np.mean(np.square(z_vel_error)))
    print('Z =', "%.6f" % error_z_vel_rms_3D, 'km/second')
    
    vel_rms = np.sqrt(error_x_vel_rms_3D**2 + error_y_vel_rms_3D**2 + error_z_vel_rms_3D**2)
    print('Overall =', "%.6f" % vel_rms, 'km/s')

    
    if saveFig_bool:
        fig_xpos.savefig('Figures/x_pos_error.png')
        fig_ypos.savefig('Figures/y_pos_error.png')
        fig_zpos.savefig('Figures/z_pos_error.png')
        fig_xvel.savefig('Figures/x_vel_error.png')
        fig_yvel.savefig('Figures/y_vel_error.png')
        fig_zvel.savefig('Figures/z_vel_error.png')
    
    
    

#phys_model is a dictionary
def calculate_area(phys_model, V, Q):
    
    #Tracer() ()
    C = functions.Quaternion_2_DCM(Q)
    
    num_facets = phys_model['num_facets']
    facet_areas = phys_model['facet_areas'] #meter^2
    normal_vecs = phys_model['facet_orientations']['normal_vecs']
    #in_plane_vecs = phys_model['facet_orientations']['in_plane_vecs']
    
    projected_area = 0
    
    for ii in range(num_facets):
        
        #rotate vector/facet using quaternion
        normal_vec_rotated = np.dot(C.T, normal_vecs[ii])
        #in_plane_vec_rotated = np.dot(C.T, in_plane_vecs[ii])
        
        #determine angle between velocity and normal vector
        val = np.dot(normal_vec_rotated, -V) / \
            (np.linalg.norm(normal_vec_rotated) * np.linalg.norm(-V))
        
        angle = functions.safe_acos(val) #because if directly aligned, sometimes val>1 or < -1
        
        if angle > math.radians(275):
            print('angle > 275')
        
        if ii == 0:
            global facet_v_angle
            facet_v_angle = math.degrees(angle)
        
        #if less than 90 degrees determine projected area
        if angle < np.radians(90):
            projected_area = projected_area + facet_areas[ii] * math.cos(angle)
    #print(projected_area)
    return projected_area

#area = calculate_area(GRACE_phys_model, v_eci)


#auto opposite dir (neg sign on V)
def calc_attitude_from_vec(X, V, phys_model):
    
    #Tracer() ()
    #normalize vectors
    V_unit = (-V/np.linalg.norm(V)).reshape(3,1)
    X_unit = (-X/np.linalg.norm(X)).reshape(3,1)
    
    (R, I, C) = ECI_2_RIC(X, V)
    I = I.reshape(3,1)
    
    radial_unit = X/np.linalg.norm(X)
    facets = phys_model['facet_orientations']['normal_vecs']
    front_facet = facets[0].reshape(3,1)
    front_facet_unit = front_facet/np.linalg.norm(front_facet)
    radial_facet = facets[-1].reshape(3,1)
    radial_facet_unit = radial_facet/np.linalg.norm(radial_facet)
    
    #allign front facet with projection of velocity vecto in orbital plane &
    #orthogonal with the radial direction
    B = np.dot(front_facet_unit, I.T) + np.dot(radial_facet_unit, X_unit.T)
    
    U, s, V = np.linalg.svd(B, full_matrices=True)
    
    M = np.diag([1, 1, np.linalg.det(U)*np.linalg.det(V)])
    #print(np.shape(U), np.shape(s), np.shape(V), np.shape(M))
    
    C = np.dot(U, np.dot(M, V))
    #print(C)
    
    
    Q = functions.DCM_2_Quaternion(C)
    
    return Q


#Q = calc_attitude_from_vec(np.array([[1], [0], [0]]), np.array([[0], [1], [0]]), phys_model)
#print(Q)


def calc_torque(X, V, phys_model, Q):
    
    C = functions.Quaternion_2_DCM(Q)
    
    num_facets = phys_model['num_facets']
    facet_areas = phys_model['facet_areas'] #meter^2
    normal_vecs = phys_model['facet_orientations']['normal_vecs']
    dist_CoM_vecs = phys_model['facet_orientations']['dist_CoM_vecs']
    
    r = np.linalg.norm(X)
    atm_density = rho_0_const * math.exp( -(r-r_0_const) / H_const )
    
    
    torque = np.zeros((3,1))
    
    
    for ii in range(num_facets):
        
        #rotate vector/facet using quaternion
        normal_vec_rotated = np.dot(C.T, normal_vecs[ii])
        
        #determine angle between velocity and normal vector
        val = np.dot(normal_vec_rotated, -V) / \
            (np.linalg.norm(normal_vec_rotated) * np.linalg.norm(-V))
        
        angle = functions.safe_acos(val) #because if directly aligned, sometimes val>1 or < -1
        #Tracer() ()
        if angle > math.radians(275):
            print('angle > 275')
        
        #if less than 90 degrees determine projected area
        if angle < np.radians(90):
            
            
            projected_area = facet_areas[ii] * math.cos(angle)
            
            force = .5 * atm_density * np.dot(V, V.T) * C_D_est * projected_area * (V/np.linalg.norm(V))
            
            
            dist_vec_rotated = np.dot(C.T, dist_CoM_vecs[ii])
            dist_norm_rotated = dist_vec_rotated/np.linalg.norm(dist_vec_rotated)
            
            #print('normal vec rotated:', normal_vec_rotated)
            #print('dist vec rotated:', dist_norm_rotated)
            #print('vel. norm:', (V/np.linalg.norm(V)))
            
            
            dum = np.cross(dist_norm_rotated, force)
            
            dum2 = np.cross(dist_CoM_vecs[ii]/np.linalg.norm(dist_CoM_vecs[ii]), np.dot(C, force))
            
            torque = torque + np.cross(dist_norm_rotated, force).reshape(3,1)
    
    return torque



def Quaternion_2_DCM(Q):
    
    q0, q1, q2, q3 = Q.reshape(4, 1)
    
    C = np.array([[q0**2 + q1**2 - q2**2 - q3**2, 2*(q1*q2 + q0*q3), 2*(q1*q3 - q0*q2)],\
                  [2*(q1*q2 - q0*q3), q0**2 - q1**2 + q2**2 - q3**2, 2*(q2*q3 + q0*q1)],\
                  [2*(q1*q3 + q0*q2), 2*(q2*q3 - q0*q1), q0**2 - q1**2 - q2**2 + q3**2]]).reshape(3,3)
                  
    return C


def safe_acos(value):
    
    if value > 1:
        value = 1
    elif value < -1:
        value = -1
    
    return math.acos(value)


#Stanley's Method to find Eulers params of DCM
#Eq. 3.94 & 3.95, Analytical Mechanics of Space Systems
def DCM_2_Quaternion(C):
    
    beta0_sqrd = .25 * (1 + np.trace(C))
    beta1_sqrd = .25 * (1 + 2*C[0, 0] - np.trace(C))
    beta2_sqrd = .25 * (1 + 2*C[1, 1] - np.trace(C))
    beta3_sqrd = .25 * (1 + 2*C[2, 2] - np.trace(C))
    
    beta_sqrd_list = [beta0_sqrd, beta1_sqrd, beta2_sqrd, beta3_sqrd]
    
    index, value = max(enumerate(beta_sqrd_list), key=operator.itemgetter(1))
    #print(beta_sqrd_list)
    
    if index == 0:
        beta0 = math.sqrt(beta0_sqrd)
        beta1 = (C[1, 2] - C[2, 1])/(4*beta0)
        beta2 = (C[2, 0] - C[0, 2])/(4*beta0)
        beta3 = (C[0, 1] - C[1, 0])/(4*beta0)
    
    elif index == 1:
        beta1 = math.sqrt(beta1_sqrd)
        beta0 = (C[1, 2] - C[2, 1])/(4*beta1)
        beta2 = (C[0, 1] + C[1, 0])/(4*beta1)
        beta3 = (C[2, 0] + C[0, 2])/(4*beta1)
    
    elif index == 2:
        beta2 = math.sqrt(beta2_sqrd)
        beta0 = (C[2, 0] - C[0, 2])/(4*beta2)
        beta1 = (C[0, 1] + C[1, 0])/(4*beta2)
        beta3 = (C[1, 2] + C[2, 1])/(4*beta2)
    
    elif index == 3:
        beta3 = math.sqrt(beta3_sqrd)
        beta0 = (C[0, 1] - C[1, 0])/(4*beta3)
        beta1 = (C[2, 0] + C[0, 2])/(4*beta3)
        beta2 = (C[1, 2] + C[2, 1])/(4*beta3)
    
    Q = np.array([beta0, beta1, beta2, beta3])
    
    return Q



def load_satellite_struct(filename, sat_name):
    
    sat_data = scipy.io.loadmat(filename)
    
    physical_model = {}
    brdf = {}
    facet_orientations = {}
    
    physical_model = sat_data[sat_name]['physical_model'][0,0][0,0]
    #print(physical_model)


    #print(physical_model.dtype)
    physical_model['num_facets'] = physical_model[0][0][0]
    physical_model['facet_areas'] = physical_model[1][0]
    
    
    brdf2 = physical_model[2]
    #print(brdf2.dtype)
    brdf['spec'] = physical_model[2][0][0][0][0]
    brdf['diff'] = physical_model[2][0][0][1][0]
    brdf['nu'] = physical_model[2][0][0][2][0]
    brdf['nv'] = physical_model[2][0][0][3][0]
    
    
    facet_orientations2 = physical_model[3]
    #print(facet_orientations2.dtype)
    facet_orientations['normal_vecs'] = physical_model[3][0][0][0]
    facet_orientations['in_plane_vecs'] = physical_model[3][0][0][1]
    facet_orientations['dist_CoM_vecs'] = physical_model[3][0][0][2]
    facet_orientations['dimensions'] = physical_model[3][0][0][3]
    
    physical_model['mass'] = physical_model[4][0][0]
    
    
    physical_model['brdf'] = brdf
    physical_model['facet_orientations'] = facet_orientations
    #print(physical_model)
    
    return physical_model


#time in seconds
def calculate_elevation(station_ecef, object_eci, time, omega, r_earth):
    station_ecef = station_ecef[0:3].reshape(3,1)
    object_eci = object_eci[0:3].reshape(3,1)

    #begin by converting object eci to ecef
    theta_GST = omega * time
    object_ecef = eci2ecef(object_eci, theta_GST)

    #calculate lat and lon of sensor
    r_earth = 6378136.3 #meters
    (lat, lon, alt) = ecef2geo_lat_lon_alt(station_ecef, r_earth)

    #transform to topo wrt to station
    #r_sez = ecef2topo(object_ecef, station_lat, station_lon)

    #calculate elevation from topo
    #(az, el, r) = topo2az_el_range(r_sez)

    rot_mat = np.array([[-math.sin(lon), math.cos(lon), 0], \
                [-math.sin(lat)*math.cos(lon), -math.sin(lat)*math.sin(lon), math.cos(lat)],\
                [math.cos(lat)*math.cos(lon), math.cos(lat)*math.sin(lon), math.sin(lat)]])
    X_L = np.dot(rot_mat, (object_ecef - station_ecef))
    az = math.atan2(X_L[0], X_L[1])
    if(az < 0):
        az = az + 2*math.pi

    el = math.asin(X_L[2]/np.linalg.norm(X_L))

    range = np.linalg.norm(object_ecef - station_ecef)

    #print(math.degrees(az), math.degrees(el), range)

    #return in degrees!
    return math.degrees(el)

#function to execute a rotation about the z axis of magnitude theta (radians)
#of the position vector r
def z_rotation(theta, r):
    rotation_mat = np.array([[math.cos(theta), math.sin(theta), 0], \
                [-math.sin(theta), math.cos(theta), 0], \
                [0, 0, 1]])
    result = np.dot(rotation_mat, r)
    return result

#perform a rotation of magniture theta GST to convert from eci to ecef
#for an obj at position pos_eci
def eci2ecef(pos_eci, theta_GST):
    r_ecef = z_rotation(theta_GST, pos_eci)
    return r_ecef

def ecef2eci(pos_ecef, theta_GST):
    r_eci = z_rotation(-theta_GST, pos_ecef)
    return r_eci

#input the ecef x, y, z position of a space obj in order to calculate
#the lat, lon, (radians) and alt of a space obj
def ecef2geo_lat_lon_alt(r_ecef, main_body_radius):
    r = np.linalg.norm(r_ecef)
    altitude = r - main_body_radius

    x, y, z = r_ecef.reshape(3, 1)
    latitude = math.asin(z/r)

    longitude = math.atan2(y, x)
    #check by using eq. for x
    #x_calc = r*math.cos(latitude)*math.cos(longitude)
    #print(x_calc - x)
    return (latitude, longitude, altitude)

#function to execute a rotation about the y axis of magnitude theta (radians)
def y_rotation_mat(theta):
    rotation_mat = np.array([[math.cos(theta), 0, -math.sin(theta)], \
                [0, 1, 0], \
                [math.sin(theta), 0, math.cos(theta)]])
    return rotation_mat

#input the ecef position of a space obj relative to a ground station
#the station is located at input params lat and lon (radians)
#to calculate the position of the space obj. in topo/sez frame
def ecef2topo(pos_ecef, lat, lon):
    r_sez = np.dot(y_rotation_mat((math.pi/2)-lat), z_rotation(lon, pos_ecef))
    return r_sez

#input the lat, lon, alt, (radians) and radius of main body
#in order to calculate the ecef x, y, z position of a ground sensor/position
def topo2ecef(latitude, longitude, altitude, main_body_radius):
    r = altitude + main_body_radius
    x_ecef = r*math.cos(latitude)*math.cos(longitude)
    y_ecef = r*math.cos(latitude)*math.sin(longitude)
    z_ecef = r*math.sin(latitude)
    r_ecef = np.array([x_ecef, y_ecef, z_ecef])
    return r_ecef

#input the position of a satellite in topo/sez frame
#in order to calculate the azimuth, elevation, and range of the satellite/object
def topo2az_el_range(pos_topo):
    s, e, z = pos_topo.reshape(3,1)
    r = np.linalg.norm(pos_topo)
    se_norm = (s**2+e**2)**(1/2)

    el_sin = z/r
    el_cos = (se_norm)/r
    el = math.atan2(el_sin, el_cos)

    az_sin = e/se_norm
    az_cos = -s/se_norm
    az = math.atan2(az_sin, az_cos)
    if(az < 0):
        az = az + 2*math.pi
    return (az, el, r)


def skew(v):
    if len(v) == 4:
        v = v[:3]/v[3]
    skv = np.roll(np.roll(np.diag(np.ndarray.flatten(v)), 1, 1), -1, 0)
    #v.flatten()

    return skv - skv.T

#print(skew(np.array([1, 2, 3])))
