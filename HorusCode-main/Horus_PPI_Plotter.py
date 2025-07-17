#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 12 09:31:51 2024

Horus Plotter 
Note: in order to run this code you need HorusFunctions.py

Created by: Laura Shedd (lshedd123@ou.edu)
"""

import HorusFunctions
import os, re,sys
import numpy as np
import matplotlib.pyplot as plt
import pyart
from pyart.config import get_metadata
import cartopy.crs as ccrs

## TODO
# Add AzShear (and other MRMS calcs)
# consider adding option for only plotting the gridded data on comparisons (limits figure burden)
# potentially add some capability for cross-sections wherever instead of along azimuth
# set up numerical comparisons?? need to determine what 
# finalize off plotting the texture fields (started but not done)
# could consider exploring the cfad option for Horus
# potential for cappi's over elevation angles for ungridded data (if possible)
    # might be possible using WDSSII and can just implement plotting

########################################
### Script Overview and Requirements ###
########################################
# This script (at present) plots Horus data and Comparisons to WSR-88D data
# Requirements: data must be in 107WDSS-II formatting
# Required Algorithms Run: w2hail, w2dualpol, w2echotop

# NOTE: the WSR-88D data have been reconfigured such that the 0.5 tilt is NOT included 
    # SAILS mode made some of the indexing more challenging so it was negated for the moment
    # But, Horus likely deals with a decent amount of GC at low tilts so...
################################
####                        #### 
#### User Defined Variables ####
####                        ####
################################10

# In Case there are time problems (goal: eventually remove need for this)
epoch = 0

# In case you need to adjust the range or something (note actually doing this will require re-running WDSSI but can use for a test)
adjust_range = False
# This needs to be in meters, not km
range_adj = -5000 

#Specify if you want to match Horus times or WSR-88D times
# False = match Horus; True = match WSR-88D
# Note: This should only be true when comparison is set to true
match_wsr = False

# If you want to plot everything or just a single PPI
plot_ppi_Horus = True
plot_all = True
plot_MESH = False
SidebySide = False
FanPlot = False ## note this one will work if plot all is False and RHI is True

# Do you want the mping for hail? 
# Note: can be modified to include other reports as well 
mping = False
# Save Figures? 
# Will do plt.show() otherwise
save_fig = True

# the date you want to use (these should match) 
#Horus_path = 'output_corrected'
# Horus_path = 'Scan13_247m'
Horus_path = 'Horus_27April2024'
date_Horus = '20240428'
vcp = '102'
mesh_ele = '00.75'
echo_thresh = ['18','30','50']
Horus_lat = [35.1864]
Horus_lon = [-97.4459]

# This is where the figures will be going (make sure its where you want!)
# figpath = '/Volumes/BRANDON/27April_Figs_4March'
figpath = '/Volumes/BRANDON/27April_Figs_27March_2'
# figpath = '/Volumes/BRANDON/Seminar_Things'
# figpath = '/Users/brandoncohen/Desktop/Horus_hold/likely_use/'
#figpath = '/Users/brandoncohen/Desktop/Horus_Annual'#+date_Horus
PPI = 0
# For April Case144
xlims = [-40,-0]
ylims = [-35,5]
# April Zoomed - BRANDON
# xlims = [-18,-10]
# ylims = [-12,-4]
# For May Case
#xlims = [-80,-20]
#ylims = [40,100]

# This one decides if we want the code to do comparisons with another radar or to just do base PPI plots
comparison = False
comp_rad = 'KTLX'
sweep_number = 1
sweep_Horus = 15
ktlx_mesh_path = '/Volumes/MiniHorus/Data/'+comp_rad+'_'+date_Horus+'/output/MESH/at00.50'
# For April Case
#lon_lim = [-98.1,-97.6]
#lat_lim = [34.7,35.1]
# April Scan 12
lon_lim = [-97.9,-97.5]
lat_lim = [34.9,35.2]
#test test on scan12
#lon_lim = [-97.9,-97.5]
#lat_lim = [34.9,35.3]
# For May Case
#lon_lim = [-98.5,-97.3]
#lat_lim = [35.4,36.3]
# Colorbar Placement for the pyart maps
cbar_ax = [-97.8, 34.8]

# Below Determines if we want to Grid and If So What Configurations to Follow
gridding = False
# GRID LIMS ARE Z, Y, X 
# For April
grid_lims = ((200, 14000), (-50000, -000), (-55000, 000))
# For May
#grid_lims = ((500,30000), (30000,110000), (-110000,-30000))
grid_x_y = 230
grid_z = 20
grid_level = 3
weighting = 'barnes2'
# Options are map_to_grid or map_gates_to_grid
grid_algo = 'map_to_grid'
roi_func = 'dist_beam'
min_radius = 100

## Below Determines if you want to plot RHIs and if so, what azimuth to complete them
# For now, these RHIs will be along a fixed azimuth, but configuration for anywhere is possible
RHI = False
az_plot = 251
zlim = [0,10] #[0,17.5]
xlim = [5,35] #[5,20]
cbar_ax_RHI = [-97.9, 34.2]
save_rhi = False
# Process radar data for each elevation angle
radar_fields = {'Reflectivity': [], 'RhoHV': [], 'ZDR': [], 'velocity': [], 'corrected_velocity': [], 'PhiDP': []}


######################################
###                                ###
###       Code Begins Here         ###
###                                ###
######################################
# First, we want to read in the elevation angles that we need
lvls1, eles = HorusFunctions.Horus_elevations(date_Horus,vcp,PPI)

# Next we want to determine what time to start at for plotting
if match_wsr is not True:
    #base_folder = '/Volumes/BRANDON/Horus_27April2024'
    # base_folder = '/Users/brandoncohen/Dropbox/Grogu/Horus_Figures_From_LauraNew/output_corrected'
    base_folder = '/Volumes/BRANDON/Horus_27April2024'
    #base_folder = '/Volumes/GROGU/Horus_Figures_From_Laura/'+Horus_path
    # base_folder = '/Volumes/GROGU/Horus/'+Horus_path
    os.chdir(base_folder+'/Reflectivity/'+lvls1[0])
    l_Hrus_index = np.sort([f for f in os.listdir('.') if re.search('.netcdf', f)])
    for j in range(0,len(l_Hrus_index)):
        print(str(j)+'. '+l_Hrus_index[j])
    timet = int(input('What Time Do You Want To Start With? '))
else:
    #base_folder = '/Users/brandoncohen/Downloads/Data_27April/'#+comp_rad+date_Horus#+'_Raw'
    base_folder = '/Users/brandoncohen/Dropbox/Grogu/Horus_Figures_From_LauraNew/output_corrected'
    match_folder = '/Volumes/GROGU/Horus_Figures_From_Laura/'+Horus_path
    os.chdir(base_folder)
    l_index = np.sort([f for f in os.listdir('.') if re.search('V06', f)])
    l_Hrus_index = l_index
    for j in range(0,len(l_index)):
        print(str(j)+'. '+l_index[j])
    timet = int(input('What Time Do You Want To Start With? '))
    
if comparison is not True:
    if plot_all is True:
        for i in range(timet,len(l_Hrus_index)):
            for j in lvls1:
                os.chdir(base_folder+'/Reflectivity/'+j)
                radar = HorusFunctions.read_Horus(date_Horus,i,j,base_folder,i,lvls1[0])
                #os.chdir(base_folder+'/MESH/at'+lvls1[0])
                # os.chdir(base_folder+'/MESH/at'+mesh_ele)
                # mesh_idx = np.sort([f for f in os.listdir('.') if re.search('.netcdf', f)])
                # rad_MESH = HorusFunctions.read_Horus_MESH(date_Horus,mesh_idx[i],'Horus')
                # os.chdir(base_folder+'/TraditionalEchoTop_'+echo_thresh[0]+'/'+mesh_ele)
                # rad_echo = HorusFunctions.read_echotop(mesh_idx[i],'Horus',echo_thresh,base_folder,mesh_ele)
                
                # Append the data for each field (Reflectivity, RhoHV, etc.) 
                radar_fields['Reflectivity'].append(radar.fields['Reflectivity']['data'])
                radar_fields['RhoHV'].append(radar.fields['RhoHV']['data'])
                radar_fields['ZDR'].append(radar.fields['ZDR']['data'])
                radar_fields['velocity'].append(radar.fields['velocity']['data'])
                radar_fields['corrected_velocity'].append(radar.fields['corrected_velocity']['data'])
                radar_fields['PhiDP'].append(radar.fields['PhiDP']['data'])
                
                if Horus_path == '20231004_PPI':
                    radar.time['data'] = np.array([1696459808])
                    rad_MESH.time['data'] = np.array([1696459808])
                if adjust_range is True:
                    radar.range['data'] = radar.range['data'] + range_adj
                # Plot the PPIs 
                if mping is True and i == timet: 
                    ping_lat,ping_lon,ping_time,ping_des = HorusFunctions.mping_read()
                    #sys.exit()
                if plot_ppi_Horus is True:
                    HorusFunctions.plot_Horus_PPI(radar,xlims,ylims,figpath,date_Horus,l_Hrus_index[i][9:15],j,echo_thresh,save_fig,plot_MESH)
                    #HorusFunctions.plot_Horus_PPI(radar,rad_MESH,rad_echo,xlims,ylims,figpath,date_Horus,l_Hrus_index[i][9:15],j,echo_thresh,save_fig,plot_MESH)
                if SidebySide is True:
                    HorusFunctions.SidebySidePlot(radar,date_Horus,l_Hrus_index[i][9:15],j,xlims,ylims,figpath,save_fig)
                
            if RHI is True:
                # First, we need to get a volume radar object for Horus
                Horus_radar_volume = HorusFunctions.Horus_Volume(date_Horus,base_folder,lvls1,i)
                if adjust_range is True:
                    Horus_radar_volume.range['data'] = Horus_radar_volume.range['data'] + range_adj
                # Now Plot the RHIs
                xsect = HorusFunctions.Horus_RHI(Horus_radar_volume,az_plot,date_Horus,l_Hrus_index[i][9:15],zlim,xlim,save_fig,figpath,xlims,ylims)
                
            fixed_angles = np.array([0.75, 2.25, 3.75, 5.25, 6.75, 8.25, 9.75, 11.25, 
                                      12.75, 14.25, 15.75, 17.25, 18.75, 20.25, 21.75, 
                                      23.25, 24.75, 26.25, 27.75, 29.25], dtype=np.float32)

            # Assign these values to the radar object
            # radar.fixed_angle = {
            #     'long_name': 'Target angle for sweep',
            #     'units': 'degrees',
            #     'standard_name': 'target_fixed_angle',
            #     'data': fixed_angles
            # }
            # radar.fixed_angle["data"] = fixed_angles.astype(np.float32)

            # # Now combine the data across all sweeps (20 sweeps in total)
            # reshaped_reflectivity = np.stack(radar_fields['Reflectivity'], axis=0)  # Shape (20, 44, 1822)
            # reshaped_rhohv = np.stack(radar_fields['RhoHV'], axis=0)  # Shape (20, 44, 1822)
            # reshaped_zdr = np.stack(radar_fields['ZDR'], axis=0)  # Shape (20, 44, 1822)
            # reshaped_velocity = np.stack(radar_fields['velocity'], axis=0)  # Shape (20, 44, 1822)
            # reshaped_corrected_velocity = np.stack(radar_fields['corrected_velocity'], axis=0)  # Shape (20, 44, 1822)
            # reshaped_phidp = np.stack(radar_fields['PhiDP'], axis=0)  # Shape (20, 44, 1822)
            
            # # Now combine the data across all sweeps (20 sweeps in total)
            # reshaped_reflectivity = np.concatenate([radar_fields['Reflectivity']], axis=0)  # Shape (20, 44, 1822)
            # reshaped_rhohv = np.concatenate([radar_fields['RhoHV']], axis=0)  # Shape (20, 44, 1822)
            # reshaped_zdr = np.concatenate([radar_fields['ZDR']], axis=0)  # Shape (20, 44, 1822)
            # reshaped_velocity = np.concatenate([radar_fields['velocity']], axis=0)  # Shape (20, 44, 1822)
            # reshaped_corrected_velocity = np.concatenate([radar_fields['corrected_velocity']], axis=0)  # Shape (20, 44, 1822)
            # reshaped_phidp = np.concatenate([radar_fields['PhiDP']], axis=0)  # Shape (20, 44, 1822)
            
            # # Ensure the reshaped data has the desired shape (20, 44, 1822)
            # reshaped_reflectivity = reshaped_reflectivity.reshape(20, 44, 1822)
            # reshaped_rhohv = reshaped_rhohv.reshape(20, 44, 1822)
            # reshaped_zdr = reshaped_zdr.reshape(20, 44, 1822)
            # reshaped_velocity = reshaped_velocity.reshape(20, 44, 1822)
            # reshaped_corrected_velocity = reshaped_corrected_velocity.reshape(20, 44, 1822)
            # reshaped_phidp = reshaped_phidp.reshape(20, 44, 1822)

            # Assuming radar_fields is a list of arrays for each elevation sweep
            reshaped_reflectivity = np.reshape(radar_fields['Reflectivity'], (20, 44, 1822))
            reshaped_rhohv = np.reshape(radar_fields['RhoHV'], (20, 44, 1822))
            reshaped_zdr = np.reshape(radar_fields['ZDR'], (20, 44, 1822))
            reshaped_velocity = np.reshape(radar_fields['velocity'], (20, 44, 1822))
            reshaped_corrected_velocity = np.reshape(radar_fields['corrected_velocity'], (20, 44, 1822))
            reshaped_phidp = np.reshape(radar_fields['PhiDP'], (20, 44, 1822))
            
            # Ensure reshaped data has the correct shape
            assert reshaped_reflectivity.shape == (20, 44, 1822), f"Unexpected shape for Reflectivity: {reshaped_reflectivity.shape}"
            assert reshaped_rhohv.shape == (20, 44, 1822), f"Unexpected shape for RhoHV: {reshaped_rhohv.shape}"
            assert reshaped_zdr.shape == (20, 44, 1822), f"Unexpected shape for ZDR: {reshaped_zdr.shape}"
            assert reshaped_velocity.shape == (20, 44, 1822), f"Unexpected shape for velocity: {reshaped_velocity.shape}"
            assert reshaped_corrected_velocity.shape == (20, 44, 1822), f"Unexpected shape for corrected_velocity: {reshaped_corrected_velocity.shape}"
            assert reshaped_phidp.shape == (20, 44, 1822), f"Unexpected shape for PhiDP: {reshaped_phidp.shape}"

            # Create the final radar object (assuming the radar has all necessary metadata)
            final_radar = radar  # Copy the original radar object structure
            finalradar = radar
            
            # Update the fields with the reshaped data
            final_radar.fields['Reflectivity']['data'] = reshaped_reflectivity
            final_radar.fields['RhoHV']['data'] = reshaped_rhohv
            final_radar.fields['ZDR']['data'] = reshaped_zdr
            final_radar.fields['velocity']['data'] = reshaped_velocity
            final_radar.fields['corrected_velocity']['data'] = reshaped_corrected_velocity
            final_radar.fields['PhiDP']['data'] = reshaped_phidp
            
            # Set the fixed_angle (elevation angles)
            final_radar.fixed_angle = {
                'long_name': 'Target angle for sweep',
                'units': 'degrees',
                'standard_name': 'target_fixed_angle',
                'data': fixed_angles
            }

            # Set sweep_start_ray_index and sweep_end_ray_index for each sweep
            n_rays_per_sweep = 44  # 44 rays per sweep
            sweep_start_indices = np.arange(0, 20 * n_rays_per_sweep, n_rays_per_sweep)
            sweep_end_indices = np.arange(n_rays_per_sweep - 1, 20 * n_rays_per_sweep, n_rays_per_sweep)
            
            final_radar.sweep_start_ray_index['data'] = sweep_start_indices
            final_radar.sweep_end_ray_index['data'] = sweep_end_indices

            radar.instrument_parameters["frequency"]["data"] = np.array([radar.instrument_parameters["frequency"]["data"]])

            # Assuming fixed_angles is already defined
            fixed_angles = np.array([0.75, 2.25, 3.75, 5.25, 6.75, 8.25, 9.75, 11.25, 
                                     12.75, 14.25, 15.75, 17.25, 18.75, 20.25, 21.75, 
                                     23.25, 24.75, 26.25, 27.75, 29.25], dtype=np.float32)
            
            # Assuming radar is your radar object and reflectivity is one of the fields
            reflectivity_data = radar.fields['Reflectivity']['data']  # shape (880, 1822)
            
            # Reshape the data to 20 sweeps, 44 rays per sweep (20, 44, 1822)
            reshaped_reflectivity = reflectivity_data.reshape(20, 44, 1822)
            
            # Now assign the reshaped reflectivity back to the radar object
            radar.fields['Reflectivity']['data'] = reshaped_reflectivity
            
            # Set the fixed angle (elevation) for each sweep
            radar.fixed_angle['data'] = fixed_angles
            
            # Number of sweeps (elevation angles)
            n_sweeps = 20
            n_rays_per_sweep = 44  # 44 rays per sweep
            
            # Assign sweep_number for each sweep
            radar.sweep_number['data'] = np.arange(1, n_sweeps + 1)
            
            # Assign sweep_start_ray_index and sweep_end_ray_index for each sweep
            sweep_start_indices = np.arange(0, n_sweeps * n_rays_per_sweep, n_rays_per_sweep)
            sweep_end_indices = np.arange(n_rays_per_sweep - 1, n_sweeps * n_rays_per_sweep, n_rays_per_sweep)
            
            radar.sweep_start_ray_index['data'] = sweep_start_indices
            radar.sweep_end_ray_index['data'] = sweep_end_indices


            # Define fixed angles (20 angles)
            # fixed_angles = np.array([0.75, 2.25, 3.75, 5.25, 6.75, 8.25, 9.75, 11.25, 
                                      # 12.75, 14.25, 15.75, 17.25, 18.75, 20.25, 21.75, 
                                      # 23.25, 24.75, 26.25, 27.75, 29.25], dtype=np.float32)
            
            # Repeat each angle 44 times to match 880 sweeps
            # fixed_angles_expanded = np.repeat(fixed_angles, 44)  # Shape (880,)
            
            # # Assign properly formatted fixed_angle to radar object
            # radar.fixed_angle = {
            #     'long_name': 'Target angle for sweep',
            #     'units': 'degrees',
            #     'standard_name': 'target_fixed_angle',
            #     'data': fixed_angles_expanded  # Ensure it matches nsweeps=880
            # }

            # print("Fixed angle data:", Horusradar.fixed_angle["data"])
            # print("Sweep numbers:", Horusradar.sweep_number["data"])


            #os.chdir(base_folder+'/NCfiles')

            # pyart.io.write_cfradial('/Users/brandoncohen/Desktop/BRANDON/Horus_27April2024/NCfiles/' + 'Horus_'+date_Horus+'_'+l_Hrus_index[i][9:15]+".nc", radar)
            
    else:
        for j in lvls1:
            os.chdir(base_folder+'/Reflectivity/'+j)
            radar = HorusFunctions.read_Horus(date_Horus,timet,j,base_folder,timet,lvls1[0])
            #os.chdir(base_folder+'/MESH/at'+lvls1[0])
            #os.chdir(base_folder+'/MESH/at'+mesh_ele)
            #mesh_idx = np.sort([f for f in os.listdir('.') if re.search('.netcdf', f)])
            #rad_MESH = HorusFunctions.read_Horus_MESH(date_Horus,mesh_idx[timet],'Horus')
            # Read in Echo Tops
            # os.chdir(base_folder+'/TraditionalEchoTop_'+echo_thresh[0]+'/'+mesh_ele)
            # rad_echo = HorusFunctions.read_echotop(mesh_idx[timet],'Horus',echo_thresh,base_folder,mesh_ele)
            if mping is True:
                    ping_lat,ping_lon,ping_time,ping_des = HorusFunctions.mping_read()
            if Horus_path == '20231004_PPI':
                radar.time['data'] = np.array([1696459808])
                rad_MESH.time['data'] = np.array([1696459808])
            # Plot the PPIs
            if plot_ppi_Horus is True:
                HorusFunctions.plot_Horus_PPI(radar,xlims,ylims,figpath,date_Horus,l_Hrus_index[timet][9:15],j,echo_thresh,save_fig,plot_MESH)

                # HorusFunctions.plot_Horus_PPI(radar,rad_MESH,rad_echo,xlims,ylims,figpath,date_Horus,l_Hrus_index[timet][9:15],j,echo_thresh,save_fig,plot_MESH)
            
            # Calculate the texture field for Horus
            #HorusFunctions.Texture_Calcs(radar,[],comparison,winsize=91)
            
        if RHI is True:
            # First, we need to get a volume radar object for Horus
            Horus_radar_volume = HorusFunctions.Horus_Volume(date_Horus,base_folder,lvls1,timet)
            # Now Plot the RHIs
            xsect = HorusFunctions.Horus_RHI(Horus_radar_volume,az_plot,date_Horus,l_Hrus_index[timet][9:15],zlim,xlim,save_fig,figpath,xlims,ylims)
            if FanPlot is True:
                HorusFunctions.fan_plots(Horus_radar_volume, radar, date_Horus, l_Hrus_index[timet][9:15], xlims, ylims, xlim, zlim, figpath, save_fig)
            
        
if comparison is True:
    if match_wsr is not True:
        if plot_all is True:
            for i in range(timet,len(l_Hrus_index)):
                # Horus data
                os.chdir(base_folder+'/Reflectivity')
                Horus_radar_volume = HorusFunctions.Horus_Volume(date_Horus,base_folder,lvls1,i)
                os.chdir(base_folder+'/MESH/at'+mesh_ele)
                mesh_idx = np.sort([f for f in os.listdir('.') if re.search('.netcdf', f)])
                Horus_rad_MESH = HorusFunctions.read_Horus_MESH(date_Horus,mesh_idx[i],'Horus')
                
                # Read KTLX data in
                #print(date_Horus)
                KTLX_base_path = '/Users/brandoncohen/Downloads/Data_27April/'+comp_rad+'_'+date_Horus#+'_Raw'
                os.chdir(KTLX_base_path)
                l_KTLX = np.sort([f for f in os.listdir('.') if re.search('V06', f)])
                time_idx = HorusFunctions.Time_Finder(l_KTLX,l_Hrus_index[i],match_wsr)
                time_ktlx_name = l_KTLX[time_idx][13:19]
                print(time_ktlx_name)
                KTLX_rad = pyart.io.read(l_KTLX[time_idx])
                # Read and Add the Dual-Pol in for KTLX
                KTLX_rad_algs = HorusFunctions.WSR_volume(date_Horus,'/Volumes/GROGU/KTLX_27April/KTLX_'+date_Horus+'/output',time_idx,comp_rad)
                #sys.exit()
                os.chdir(ktlx_mesh_path)
                mesh_idx_ktlx = np.sort([f for f in os.listdir('.') if re.search('.netcdf', f)])
                #KTLX_rad_MESH = HorusFunctions.read_Horus_MESH(mesh_idx_ktlx[time_idx],comp_rad)
                if time_ktlx_name[0] == '0' and date_Horus == '20240506':
                    date_Horus = '20240507'
                KTLX_rad_MESH = HorusFunctions.read_Horus_MESH(date_Horus,date_Horus+'-'+time_ktlx_name+'.netcdf',comp_rad)
                #KTLX_rad_MESH = HorusFunctions.read_HOrus_MESH(mesh_idk_ktlx[i])
                if date_Horus == '20240428':
                    print('Adjusting Azimuth...')
                    Horus_radar_volume.azimuth['data'] = Horus_radar_volume.azimuth['data'] + 5
                    print('Adjusting Reflectivity...')
                    Horus_radar_volume.fields['Reflectivity']['data'] = Horus_radar_volume.fields['Reflectivity']['data'] - 5
                # Plot
                HorusFunctions.Comparison_Plots(Horus_radar_volume,KTLX_rad,Horus_rad_MESH,KTLX_rad_MESH,KTLX_rad_algs,
                                    lon_lim,lat_lim,figpath,date_Horus,l_Hrus_index[i][9:15],cbar_ax,sweep_number,sweep_Horus,save_fig)
                if date_Horus == '20240507':
                    print('Adjusting Date Back')
                    date_Horus = '20240506'
                
        else:
                    
            KTLX_base_path = '/Users/brandoncohen/Downloads/Data_27April/'+comp_rad+'_'+date_Horus+'V06'#+'_Raw'
            
            os.chdir(base_folder+'/Reflectivity')
            Horus_radar_volume = HorusFunctions.Horus_Volume(date_Horus,base_folder,lvls1,timet)
        
        #os.chdir(base_folder+'/MESH/at'+lvls1[0])
            os.chdir(base_folder+'/MESH/at'+mesh_ele)
            mesh_idx = np.sort([f for f in os.listdir('.') if re.search('.netcdf', f)])
            Horus_rad_MESH = HorusFunctions.read_Horus_MESH(date_Horus,mesh_idx[timet],'Horus')
        
        
            if Horus_path == '20231004_PPI':
                Horus_radar_volume.time['data'] = np.array([1696459808])
                Horus_rad_MESH.time['data'] = np.array([1696459808])
            #if date_Horus == '20240428':
            ##        print('Adjusting Azimuth...')
            #        Horus_radar_volume.azimuth['data'] = Horus_radar_volume.azimuth['data'] + 5
            #        print('Adjustng Range...')
            #        Horus_radar_volume.range['data'] = Horus_radar_volume.range['data'] - 5000
            #        print('Adjusting Reflectivity...')
            #        Horus_radar_volume.fields['Reflectivity']['data'] = Horus_radar_volume.fields['Reflectivity']['data'] - 5
            #if date_Horus == '20240428':
                #print('Adjusting Azimuth...')
                #print(np.nanmin(Horus_radar_volume.range['data']))
                #Horus_radar_volume.azimuth['data'] = Horus_radar_volume.azimuth['data'] + 30
                #Horus_rad_MESH.azimuth['data'] = Horus_rad_MESH.azimuth['data'] + 30
            #    Horus_radar_volume.range['data'] = Horus_radar_volume.range['data'] - 5000
            #    Horus_radar_volume.init_gate_x_y_z()
                #print(np.nanmin(Horus_radar_volume.range['data']))
    
        # Now that the Horus radar object is setup, we need to find the closest corresponding time to KTLX
        # First, we need to swap into the KTLX (or KCRI) Directory
            os.chdir(KTLX_base_path)
            l_KTLX = np.sort([f for f in os.listdir('.') if re.search('V06', f)])
            
            # Now run our time matching function to find the closest time and output the KTLX radar object
            time_idx = HorusFunctions.Time_Finder(l_KTLX,l_Hrus_index[timet],match_wsr)
            time_ktlx_name = l_KTLX[time_idx][13:19]
            if time_ktlx_name[0] == '0' and date_Horus == '20240506':
                date_Horus = str(int(date_Horus) + 1)
            KTLX_rad = pyart.io.read(l_KTLX[time_idx])
        
        # Read MESH in
            os.chdir(ktlx_mesh_path)
            mesh_idx_ktlx = np.sort([f for f in os.listdir('.') if re.search('.netcdf', f)])
            #KTLX_rad_MESH = HorusFunctions.read_Horus_MESH(mesh_idx_ktlx[time_idx],comp_rad)
            KTLX_rad_MESH = HorusFunctions.read_Horus_MESH(date_Horus,date_Horus+'-'+time_ktlx_name+'.netcdf',comp_rad)
            KTLX_rad_algs = HorusFunctions.WSR_volume(date_Horus,'/Volumes/MiniHorus/Data/KTLX_'+date_Horus+'/output',time_idx,comp_rad)
    else:
        
        #Read in the WSR Radar data and MESH
        KTLX_rad = pyart.io.read(l_index[timet])
        
        # #TO DO
        # # need to set it up such that it matches the time for the specific sweep, not the volume
        # os.chdir(ktlx_mesh_path)
        # mesh_idx_ktlx = np.sort([f for f in os.listdir('.') if re.search('.netcdf', f)])
        # KTLX_rad_MESH = HorusFunctions.read_Horus_MESH(date_Horus,mesh_idx_ktlx[timet],comp_rad)
        
        #Now we want to find the closest Horus time
        os.chdir(match_folder+'/Reflectivity/'+lvls1[0])
        l_Horus = np.sort([f for f in os.listdir('.') if re.search('.netcdf', f)])
        time_idx = HorusFunctions.Time_Finder(l_Horus,l_index[timet],match_wsr)
        #print(time_idx)
        print('Horus Time',l_Horus[time_idx])
        # Read the data in
        Horus_radar_volume = HorusFunctions.Horus_Volume(date_Horus,match_folder,lvls1,time_idx)
        if adjust_range is True:
            Horus_radar_volume.range['data'] = Horus_radar_volume.range['data'] + range_adj
        os.chdir(match_folder+'/MESH/at'+mesh_ele)
        mesh_idx = np.sort([f for f in os.listdir('.') if re.search('.netcdf', f)])
        Horus_rad_MESH = HorusFunctions.read_Horus_MESH(date_Horus,mesh_idx[time_idx],'Horus')
        
        #sys.exit()
        
    # Plot the Baseline Comparisons
    if plot_all is not True: #this ensures that we don't replot what was made above
        HorusFunctions.Comparison_Plots(Horus_radar_volume,KTLX_rad,KTLX_rad_algs,
                                    lon_lim,lat_lim,figpath,date_Horus,l_Hrus_index[timet][9:15],cbar_ax,sweep_number,sweep_Horus,save_fig)
    
    if gridding:
        grid_KTLX, grid_Horus = HorusFunctions.Gridding(Horus_radar_volume,KTLX_rad,
                                                        eles,grid_lims,lon_lim,lat_lim,
                                                        grid_x_y,grid_z,grid_level,weighting,roi_func,min_radius,grid_algo,cbar_ax,save_fig)
        
        refl1 = grid_KTLX.fields['reflectivity']['data']
        refl2 = grid_Horus.fields['Reflectivity']['data']
        
        zdr1 = grid_KTLX.fields['differential_reflectivity']['data']
        zdr2 = grid_Horus.fields['ZDR']['data']
        zdr1[zdr1<-100] = np.nan
        zdr2[zdr2<-100] = np.nan
        
        # Find the Reflectivity Difference with the Gridding Data and add it into one of the radars
        refl_dif = refl1 - refl2
        ref_diff_dict = get_metadata('Reflectivity Difference')
        ref_diff_dict['data'] = refl_dif
        ref_diff_dict['units'] = 'dBZ'
        grid_Horus.add_field('Reflectivity Difference',ref_diff_dict)
        
        #Plot the Differences
        fig = plt.figure(figsize=(18, 6))
        ax1 = plt.subplot(111, projection=ccrs.PlateCarree())
        display = pyart.graph.GridMapDisplay(grid_Horus)
        #display_grid2 = pyart.graph.GridMapDisplay(grid2)
        #display_grid3 = pyart.graph.GridMapDisplay(grid3)
        pcm = display.plot_grid(
            "Reflectivity Difference",
            level = grid_level,
            ax=ax1,
            cmap="pyart_NWSRef",
            vmin=-10,
            vmax=20,
            colorbar_flag=True
            )
        plt.xlim(lon_lim)
        plt.ylim(lat_lim)
        
        # Filter Out the Nan
        refl1_lvl1 = refl1[:,:,:]
        refl2_lvl1 = refl2[:,:,:]
        zdr1_lvl1 = zdr1[:,:,:]
        zdr2_lvl1 = zdr2[:,:,:]
        
        bad_indices = np.isnan(refl1_lvl1) | np.isnan(refl2_lvl1)
        bad_indices_zdr = np.isnan(zdr1_lvl1) | np.isnan(zdr2_lvl1) 
        good_indices = ~bad_indices
        good_indices_zdr = ~bad_indices_zdr
        good_refl1 = refl1_lvl1[good_indices]
        good_refl2 = refl2_lvl1[good_indices]
        good_zdr1 = zdr1_lvl1[good_indices_zdr]
        good_zdr2 = zdr2_lvl1[good_indices_zdr]
        
        # Plot a 2D Histogram to Compare the Reflectivity 
        fig = plt.figure()
        plt.hist2d(good_refl1.flatten(),good_refl2.flatten(),bins = 30,cmap='jet',density=True,cmin=0.00001)
        plt.plot([1,70],[1,70],color='k')
        plt.colorbar()
        plt.xlim([10,70])
        plt.ylim([10,70])
        plt.xlabel('KTLX Gridded Reflectivity')
        plt.ylabel('Horus Gridded Reflectivity')
        
        fig = plt.figure()
        plt.hist2d(good_zdr1.flatten(),good_zdr2.flatten(),bins = 300,cmap='jet',density=True,cmin=0.05)
        plt.plot([-1,3],[-1,3],color='k')
        plt.colorbar()
        plt.xlim([-1,3])
        plt.ylim([-1,3])
        plt.xlabel('KTLX Gridded ZDR')
        plt.ylabel('Horus Gridded ZDR')
        
    if RHI == True:
        # IMPORTANT: for RHI comparisons, the data MUST be the gridded datasets 
        # So for the RHI comparisons, comparisons AND gridding must be set to true
        if gridding is not True:
            print('Gridding Must be on, closing program')
            sys.exit()

        #xsect = HorusFunctions.Horus_RHI(Horus_radar_volume,az_plot,date_Horus,l_Hrus_index[timet][9:15],zlim,save_fig)
        HorusFunctions.RHI_Comparisons(Horus_radar_volume, grid_Horus, grid_KTLX, comp_rad, az_plot,date_Horus,l_Hrus_index[timet][9:15],cbar_ax_RHI,figpath,save_fig)
       