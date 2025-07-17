#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 12 10:08:44 2024

Horus Function Folder

@author: lshedd123
"""
import os ,re
import pyart
from datetime import datetime,timezone,timedelta
import netCDF4
import numpy as np
from pyart.config import get_metadata
import matplotlib.pyplot as plt
from matplotlib import colors
import cartopy.crs as ccrs
import matplotlib.cm as cm
import requests
import geopy.distance
import numpy.matlib

# Useful for plotting correlation coefficient 
cdict = {'red':   [[0.0,  0.0, 0.0],
                   [0.5,  0.0, 0.0],
                   [0.85,  1.0, 1.0],
                   [1.0,  1.0, 1.0]],
         'green': [[0.0,  0.0, 0.0],
                   [0.25, 0.0, 0.0],
                   [0.5,  0.0, 0.0],
                   [0.85,  1.0, 1.0],
                   [1.0,  0.0, 0.0]],
         'blue':  [[0.0,  1.0, 1.0],
                   [0.5,  1.0, 1.0],
                   [1.0,  0.1, 0.1]]}
cc_cmap = colors.LinearSegmentedColormap('testCmap', segmentdata=cdict, N=256)

ticks = np.arange(1,15)
labs = ['L/MR','HR','R/Ha','BD','AP','BI','UK','NE','DS','WS','CY','','GR','+']
hca_cmap = colors.ListedColormap(['#AAAAAA','#0000FE','#0057FF','#00A9FF','#00FFFF','#56FFAA','#ACFD56','#FFFF00','#FEAB02','#FF5700','#FF00FF','#B300B3','#DC143C','#FF00FF'])
hca_cmap = colors.ListedColormap(['limegreen','green','#DC143C','#FFFF00','dimgray','lightgray','#B300B3','white','#00FFFF','#0000FE','mediumturquoise','#000000','#FF00FF','#000000'])

def Horus_elevations(dateHorus,vcp,PPI):
    if dateHorus == '20230511' or dateHorus == '20231004' or dateHorus == '20230808':
        lvls1 = ['00.50','01.00','01.50','02.00','02.50','03.00','03.50','04.00','04.50','05.00','05.50','06.00',
                   '06.50','07.00','07.50','08.00','08.50','09.00','09.50','10.00','10.50','11.00','11.50','12.00',
                   '12.50','13.00','13.50','14.00', '14.50','15.00','15.50','16.00','16.50','17.00','17.50','18.00',
                   '18.50','19.00','19.50','20.00','20.50','21.00','21.50','22.00','22.50','23.00','23.50','24.00',
                   '24.50','25.00','25.50','26.00','26.50','27.00','27.50','28.00','28.50','29.00','29.50','30.00',
                   '30.50']
        eles = [00.50,01.00,01.50,02.00,02.50,03.00,03.50,04.00,04.50,05.00,05.50,06.00,
                   06.50,07.00,07.50,08.00,08.50,09.00,09.50,10.00,10.50,11.00,11.50,12.00,
                   12.50,13.00,13.50,14.00, 14.50,15.00,15.50,16.00,16.50,17.00,17.50,18.00,
                   18.50,19.00,19.50,20.00,20.50,21.00,21.50,22.00,22.50,23.00,23.50,24.00,
                   24.50,25.00,25.50,26.00,26.50,27.00,27.50,28.00,28.50,29.00,29.50,30.00,
                   30.50]
    elif dateHorus == '20230216':
        lvls1 = ['00.00','01.33','02.67','04.00','05.33','06.67','08.00','09.33','10.67',
                 '12.00','13.33','14.67','16.00','17.33','18.67']
        eles = [00.00,01.33,02.67,04.00,05.33,06.67,08.00,09.33,10.67,
                 12.00,13.33,14.67,16.00,17.33,18.67]
    elif dateHorus == '20221213':
        lvls1 = ['00.00','01.00','02.00','03.00','04.00','05.00','06.00'
                   ,'07.00','08.00','09.00','10.00','11.00','12.00'
                   ,'13.00','14.00','15.00','16.00','17.00','18.00'
                   ,'19.00','20.00','21.00','22.00','23.00','24.00'
                   ,'25.00','26.00','27.00','28.00','29.00','30.00']
        eles = [00.00,01.00,02.00,03.00,04.00,05.00,06.00,
                   07.00,08.00,09.00,10.00,11.00,12.00,
                   13.00,14.00,15.00,16.00,17.00,18.00,
                   19.00,20.00,21.00,22.00,23.00,24.00,
                   25.00,26.00,27.00,28.00,29.00,30.00]
    elif vcp == '101':
        # Still having problems with the 00.62 elevation angle that may need to be resovled but removing for now
        lvls1 = ['00.62','01.88', '03.12', '04.38', '05.62', '06.88', '08.12', '09.38', '10.62', '11.88',
                '13.12', '14.38', '15.62', '16.88', '18.12', '19.38']
        eles = [ 00.62,01.88, 03.12, 04.38, 05.62, 06.88, 08.12, 09.38, 10.62, 11.88,
                13.12, 14.38, 15.62, 16.88, 18.12, 19.38]
    elif vcp == '102':
        lvls1 = ['00.75','02.25','03.75','05.25','06.75','08.25','09.75','11.25','12.75','14.25','15.75',
                 '17.25','18.75','20.25','21.75','23.25','24.75','26.25','27.75','29.25']
        eles = [00.75,02.25,03.75,05.25,06.75,08.25,09.75,11.25,12.75,14.25,15.75,
                 17.25,18.75,20.25,21.75,23.25,24.75,26.25,27.75,29.25]
    if PPI == 1:
            lvls1 = ['02.50','10.00','12.50']
            eles = [02.50,10.00,12.50]
    return lvls1,eles

def read_Horus(dateHorus,idz,elevations,base,HCAtime,low_ele):
    os.chdir(base+'/Reflectivity/'+elevations)
    l2=np.sort([f for f in os.listdir('.') if re.search('.netcdf', f)])
    file_path=l2[idz]
    #print(file_path)
    with netCDF4.Dataset(file_path,mode='r') as data:
            try:
                Ref=np.array(data['Reflectivity'])
            except: 
                Ref = np.array(data['Reflectivity_EXTRA']) 
            I,J=Ref.shape
            #print(I,J)
            #Az=np.array(nc['Azimuth'])

            # Fill the Radar Object with important fields
            rad = pyart.testing.make_empty_ppi_radar(J, I, 1)
            rad.ngates = J
            rad.nrays = I
            time_dict=get_metadata('Time')
            #if epoch == 1: 
            #    time_1 = datetime.fromtimestamp(int(np.array([data.Time])))
            #    time_1 = time_1 - timedelta(hours=5, minutes=0)
            #    time_dict['data'] = np.array([int(time_1.timestamp())])
            #else:
            time_dict['data']=np.array(np.repeat(data.Time,I))
            #print(data.Time)
            #sys.exit()
            time_dict['units']='seconds since 1970-01-01T00:00:00'
            time_dict['calendar']='gregorian'
            rad.time = time_dict
            if dateHorus == '20240428':
                rad.latitude['data'] = np.array([35.1864])
                rad.longitude['data'] = np.array([-97.4459])
            else:
                rad.latitude['data'] = np.array([data.Latitude])
                rad.longitude['data'] = np.array([data.Longitude])
            rad.sweep_number['data'] = np.array([1])
            rad.azimuth['data'] = np.array(data['Azimuth'])
            #rad.azimuth['data'] = np.array(60.25)
            #azs = rad.azimuth['data'][0]
            #rad.azimuth['data'][0:-1] = rad.azimuth['data'][1:]
            #rad.azimuth['data'][-1]=azs
            rad.fixed_angle['data']=np.array([elevations]).astype(float)
            #print(np.ones(I)*np.array([elevations]).astype(float))
            rad.elevation['data'] = np.ones(I) * np.array([elevations]).astype(float)
            #rad.sweep_mode['data'] = np.array(data['sweep_mode'])
            rad.instrument_name='Horus'
            rad.metadata['instrument_name']='Horus'
            parameters = {}
            frequency = {'data': int(2705000000), 'units': 'Hertz'}
            nyquistvel = {'data': np.array(data['NyquistVelocity'])[0],'units':'m/s'}
            parameters['frequency'] = frequency
            parameters['nyquist_velocity'] = nyquistvel
            rad.instrument_parameters=parameters
            
            
            
            # Lets Have Fun With Ranges
            numgates=(len(data.dimensions['Gate']))
            gatewidth=np.array(data['GateWidth'])[0]
            start=data.RangeToFirstGate
            ranges=np.ones(numgates)+start
            for i in range(1,len(ranges)):
                ranges[i]=ranges[i-1]+gatewidth
            
            range_dict=get_metadata('Range')
            range_dict['data']=ranges
            range_dict['Units'] = 'Meters'
            rad.range=range_dict
            
            rad.init_gate_altitude()
            rad.init_gate_longitude_latitude()
            
            # Now get other data - ie Radar Fields
            # Reflectivity
            ref_dict = get_metadata('Reflectivity')
            try:
                ref_dict['data'] = np.array(data['Reflectivity'])
                ref_dict['data'][ref_dict['data']<-10000] = np.nan
            except:
                ref_dict['data'] = np.array(data['Reflectivity_EXTRA'])
            ref_dict['units']='dBZ'
            rad.add_field('Reflectivity', ref_dict)
            
            #RhoHV
            os.chdir(base+'/RhoHV/'+elevations)
            l2=np.sort([f for f in os.listdir('.') if re.search('.netcdf', f)])
            file_path=l2[idz]
            with netCDF4.Dataset(file_path,mode='r') as data:
                rho_dict = get_metadata('RhoHV')
                try:
                    rho_dict['data'] = np.array(data['RhoHV'])
                    rho_dict['data'][rho_dict['data']==-99900] = np.nan
                except:
                    rho_dict['data'] = np.array(data['RhoHV_EXTRA'])
                rho_dict['units'] = 'unitless'
                rad.add_field('RhoHV',rho_dict)
            
            #ZDR
            print(base+'/Zdr/'+elevations)
            os.chdir(base+'/Zdr/'+elevations)
            l2=np.sort([f for f in os.listdir('.') if re.search('.netcdf', f)])
            file_path=l2[idz]
            with netCDF4.Dataset(file_path,mode='r') as data:
                zdr_dict = get_metadata('Zdr')
                try:
                    zdr_dict['data'] = np.array(data['Zdr'])
                    zdr_dict['data'][zdr_dict['data']==-99900] = np.nan
                except:
                    zdr_dict['data'] = np.array(data['Zdr_EXTRA'])  
                zdr_dict['units'] = 'dB'
                rad.add_field('ZDR',zdr_dict)
                
            # Velocity
            os.chdir(base+'/AliasedVelocity/'+elevations)
            l2=np.sort([f for f in os.listdir('.') if re.search('.netcdf', f)])
            file_path=l2[idz]
            with netCDF4.Dataset(file_path,mode='r') as data:
                zdr_dict = get_metadata('velocity')
                try:
                    zdr_dict['data'] = np.array(data['AliasedVelocity'])
                    zdr_dict['data'][zdr_dict['data']==-99900] = np.nan
                except:
                    zdr_dict['data'] = np.array(data['AliasedVelocity_EXTRA'])
                zdr_dict['units'] = 'm/s'
                rad.add_field('velocity',zdr_dict)
            # create a gate filter which specifies gates to exclude from dealiasing
            rad.instrument_parameters["nyquist_velocity"]={"data":np.resize([24.41],rad.fields['Reflectivity']['data'].shape[0])}
            gatefilter = pyart.filters.GateFilter(rad)
            gatefilter.exclude_transition()
            gatefilter.exclude_invalid("velocity")
            gatefilter.exclude_invalid("Reflectivity")
            gatefilter.exclude_outside("Reflectivity", -10, 90)

            # perform dealiasing
            dealias_data = pyart.correct.dealias_region_based(rad, gatefilter=gatefilter)
            rad.add_field("corrected_velocity", dealias_data)
            
            # PhiDP
            os.chdir(base+'/PhiDP/'+elevations)
            l2=np.sort([f for f in os.listdir('.') if re.search('.netcdf', f)])
            file_path=l2[idz]
            with netCDF4.Dataset(file_path,mode='r') as data:
                phi_dict = get_metadata('PhiDP')
                phi_dict['data'] = np.array(data['PhiDP'])
                phi_dict['data'][phi_dict['data']==-99900] = np.nan
                phi_dict['units'] = 'degrees'
                rad.add_field('PhiDP',phi_dict)
            # MESH
            '''
            os.chdir(base+'/MESH/at'+low_ele)
            file_path=idz
            with netCDF4.Dataset(file_path,mode='r') as data:
                zdr_dict = get_metadata('MESH')
                zdr_dict['data'] = np.array(data['MESH'])
                print(np.array(data['MESH']).shape)
                zdr_dict['units'] = 'mm'
                rad.add_field('MESH',zdr_dict)'''
            # # HCA
            # os.chdir(base+'/DHCA/'+elevations)
            # l2=np.sort([f for f in os.listdir('.') if re.search('.netcdf', f)])
            # file_path=l2[HCAtime]
            # #print(file_path)
            # with netCDF4.Dataset(file_path,mode='r') as data:
            #     zdr_dict = get_metadata('HCA')
            #     zdr_dict['data'] = np.array(data['DHCA'])
            #     zdr_dict['units'] = 'unitless'
            #     rad.add_field('HCA',zdr_dict)
                
            # # KDP
            # os.chdir(base+'/DKDP/'+elevations)
            # l2=np.sort([f for f in os.listdir('.') if re.search('.netcdf', f)])
            # file_path=l2[HCAtime]
            # #print(file_path)
            # with netCDF4.Dataset(file_path,mode='r') as data:
            #     zdr_dict = get_metadata('KDP')
            #     zdr_dict['data'] = np.array(data['DKDP'])
            #     zdr_dict['units'] = 'deg/km'
            #     rad.add_field('KDP',zdr_dict)
            
    return rad

def read_Horus_MESH(dateHorus,idz,rad_name):
    file_path=idz
    with netCDF4.Dataset(file_path,mode='r') as data:
            Ref=np.array(data['MESH'])
            I,J=Ref.shape
            #Az=np.array(nc['Azimuth'])

            # Fill the Radar Object with important fields
            rad = pyart.testing.make_empty_ppi_radar(J, I, 1)
            time_dict=get_metadata('Time')
            #if epoch == 1: 
            #    time_1 = datetime.fromtimestamp(int(np.array([data.Time])))
            #    time_1 = time_1 - timedelta(hours=5, minutes=0)
            #    time_dict['data'] = np.array([int(time_1.timestamp())])
            #else:
            time_dict['data']=np.array([data.Time])
            #print(data.Time)
            #sys.exit()
            time_dict['units']='seconds since 1970-01-01T00:00:00'
            time_dict['calendar']='gregorian'
            rad.time = time_dict
            if dateHorus == '20240428' and rad_name=='Horus':
                rad.latitude['data'] = np.array([35.1864])
                rad.longitude['data'] = np.array([-97.4459])
            else:
                rad.latitude['data'] = np.array([data.Latitude])
                rad.longitude['data'] = np.array([data.Longitude])
            rad.sweep_number['data'] = np.array([1])
            rad.azimuth['data'] = np.array(data['Azimuth'])
            #rad.azimuth['data'] = np.array(60.25)
            #azs = rad.azimuth['data'][0]
            #rad.azimuth['data'][0:-1] = rad.azimuth['data'][1:]
            #rad.azimuth['data'][-1]=azs
            rad.fixed_angle['data']=np.array([data.Elevation])
            rad.elevation['data'] = np.array([data.Elevation])
            #rad.sweep_mode['data'] = np.array(data['sweep_mode'])
            rad.instrument_name=rad_name
            rad.metadata['instrument_name']=rad_name
            parameters = {}
            frequency = {'data': int(2705000000), 'units': 'Hertz'}
            parameters['frequency'] = frequency
            rad.instrument_parameters=parameters
            
            # Lets Have Fun With Ranges
            numgates=(len(data.dimensions['Gate']))
            gatewidth=np.array(data['GateWidth'])[0]
            start=data.RangeToFirstGate
            ranges=np.ones(numgates)+start
            for i in range(1,len(ranges)):
                ranges[i]=ranges[i-1]+gatewidth
            
            range_dict=get_metadata('Range')
            range_dict['data']=ranges
            range_dict['Units'] = 'Meters'
            rad.range=range_dict
            
            rad.init_gate_altitude()
            rad.init_gate_longitude_latitude()
            
            # Now get other data - ie Radar Fields
            # Reflectivity
            ref_dict = get_metadata('MESH')
            try:
                ref_dict['data'] = np.array(data['MESH'])
                ref_dict['data'][ref_dict['data']==0] = np.nan
            except:
                ref_dict['data'] = np.array(data['MESH_EXTRA'])
            ref_dict['units']='mm'
            rad.add_field('MESH', ref_dict)
            
            #if rad_name == 'KTLX' or rad_name == 'KCRI':
            #    os.chdir(base + '/DHCA/'+)
            return rad

def Horus_Volume(dateHorus,base_folder,lvls1,Hrus_file):
    Hrus_rads = []
    sweep_start = []
    sweep_end = []
    for j in range(0,len(lvls1)):
        os.chdir(base_folder+'/Reflectivity/'+lvls1[j])
        Hrusrad = read_Horus(dateHorus,Hrus_file,lvls1[j],base_folder,Hrus_file,lvls1[0])
        nrays = Hrusrad.nrays
        sweep_start.append(j*nrays)
        sweep_end.append(((j+1)*nrays)-1)
        Hrus_rads.append(Hrusrad)
    while len(Hrus_rads) > 1:
            bob = Hrus_rads.pop(1)
        #print(bob.nrays)
            Hrus_rads[0] = pyart.util.join_radar(Hrus_rads[0], bob)
    #print(sweep_start)
    #print(sweep_end)
    Horusradar=Hrus_rads[0]
    Horusradar.sweep_number['data'] = np.array([np.linspace(1,len(lvls1),num=len(lvls1))])
    Horusradar.nrays = len(Horusradar.azimuth['data'])
    Horusradar.sweep_start_ray_index['data'] = sweep_start
    Horusradar.sweep_end_ray_index['data'] = sweep_end
    return Horusradar

# def Horus_Volume(dateHorus, base_folder, lvls1, Hrus_file):
#     Hrus_rads = []  # List to hold radar objects for each elevation angle
#     sweep_start = []  # Track starting index of each sweep
#     sweep_end = []  # Track ending index of each sweep
    
#     for j in range(0, len(lvls1)):
#         os.chdir(base_folder+'/Reflectivity/'+lvls1[j])
#         Hrusrad = read_Horus(dateHorus, Hrus_file, lvls1[j], base_folder, Hrus_file, lvls1[0])
        
#         # Check if Hrusrad is valid (i.e., not None or empty)
#         if Hrusrad is None:
#             print(f"Error: Radar object is None for level {lvls1[j]}")
#             continue
        
#         nrays = Hrusrad.nrays
#         sweep_start.append(j * nrays)
#         sweep_end.append(((j + 1) * nrays) - 1)
#         Hrus_rads.append(Hrusrad)
    
#     # Check if Hrus_rads has been populated correctly
#     if len(Hrus_rads) == 0:
#         print("Error: No radar objects were added to Hrus_rads.")
#         return None  # Exit if no radar objects
    
#     # Join radar objects (sweeps from different elevations)
#     while len(Hrus_rads) > 1:
#         bob = Hrus_rads.pop(1)
#         Hrus_rads[0] = pyart.util.join_radar(Hrus_rads[0], bob)
    
#     # Final radar object with all sweeps
#     Horusradar = Hrus_rads[0]
    
#     # Set the sweep number for each sweep (should be the length of lvls1)
#     Horusradar.sweep_number['data'] = np.array([np.linspace(1, len(lvls1), num=len(lvls1))]).astype(np.int32)
    
#     # Update sweep start and end ray indices
#     Horusradar.sweep_start_ray_index['data'] = np.array(sweep_start)
#     Horusradar.sweep_end_ray_index['data'] = np.array(sweep_end)

#     # Now, make sure you update the fixed_angle with the correct array for all sweeps:
#     fixed_angles = np.array([0.75, 2.25, 3.75, 5.25, 6.75, 8.25, 9.75, 11.25, 
#                              12.75, 14.25, 15.75, 17.25, 18.75, 20.25, 21.75, 
#                              23.25, 24.75, 26.25, 27.75, 29.25], dtype=np.float32)

#     # Create an expanded fixed_angles array that matches the number of sweeps
#     fixed_angles_expanded = np.repeat(fixed_angles, len(lvls1))

#     # Update the fixed_angle data in radar
#     Horusradar.fixed_angle = {'data': fixed_angles_expanded}

#     return Horusradar

# def Horus_Volume(dateHorus, base_folder, lvls1, Hrus_file):
#     Hrus_rads = []
#     sweep_start = []
#     sweep_end = []
    
#     # Load each radar data for each elevation (sweep)
#     for j in range(0, len(lvls1)):
#         os.chdir(base_folder+'/Reflectivity/'+lvls1[j])
#         Hrusrad = read_Horus(dateHorus, Hrus_file, lvls1[j], base_folder, Hrus_file, lvls1[0])
#         nrays = Hrusrad.nrays
#         sweep_start.append(j * nrays)
#         sweep_end.append(((j + 1) * nrays) - 1)
#         Hrus_rads.append(Hrusrad)

#     # Merge all sweeps (radar data) into one radar object
#     while len(Hrus_rads) > 1:
#         bob = Hrus_rads.pop(1)
#         Hrus_rads[0] = pyart.util.join_radar(Hrus_rads[0], bob)

#     # Final merged radar object
#     Horusradar = Hrus_rads[0]
    
#     # Assign sweep indices and number of sweeps
#     Horusradar.sweep_number['data'] = np.array([np.linspace(1, len(lvls1), num=len(lvls1))])
#     Horusradar.nrays = len(Horusradar.azimuth['data'])
#     Horusradar.sweep_start_ray_index['data'] = sweep_start
#     Horusradar.sweep_end_ray_index['data'] = sweep_end

#     # Now, return the fully populated radar object with all sweeps
#     return Horusradar


def WSR_volume(dateHorus,base_folder,time_idx,rad_name):
    WSR_rads = []
    os.chdir(base_folder+'/Reflectivity')
    eles = np.array(os.listdir())
    eles = eles[3:]
    for j in range(0,len(eles)):
        #print(eles)
        os.chdir(base_folder+'/Reflectivity/'+eles[j])
        radar = read_WSR_volume(dateHorus,base_folder,eles[j],time_idx,rad_name)
        WSR_rads.append(radar)
    while len(WSR_rads) > 1:
         #print(len(WSR_rads))
         bob = WSR_rads.pop(1)
         WSR_rads[0] = pyart.util.join_radar(WSR_rads[0],bob)
    WSR_rad = WSR_rads[0]
    
    return WSR_rad

def read_WSR_volume(dateHorus,base_path,ele,time_idx,rad_name):
    os.chdir(base_path+'/DHCA/'+ele)
    # First find the right elevation
    l2=np.sort([f for f in os.listdir('.') if re.search('.netcdf', f)])
    file_path=l2[time_idx]
    #print(file_path)
    with netCDF4.Dataset(file_path,mode='r') as data:
            
            Ref=np.array(data['DHCA'])
            I,J=Ref.shape
            # Fill the Radar Object with important fields
            rad = pyart.testing.make_empty_ppi_radar(J, I, 1)
            rad.ngates = J
            rad.nrays = I
            time_dict=get_metadata('Time')
            time_dict['data']=np.array(np.repeat(data.Time,I))
            #print(data.Time)
            time_dict['units']='seconds since 1970-01-01T00:00:00'
            time_dict['calendar']='gregorian'
            rad.time = time_dict
            rad.latitude['data'] = np.array([data.Latitude])
            rad.longitude['data'] = np.array([data.Longitude])
            rad.sweep_number['data'] = np.array([1])
            rad.azimuth['data'] = np.array(data['Azimuth'])
            #rad.azimuth['data'] = np.array(60.25)
            #azs = rad.azimuth['data'][0]
            #rad.azimuth['data'][0:-1] = rad.azimuth['data'][1:]
            #rad.azimuth['data'][-1]=azs
            rad.fixed_angle['data']=np.array([ele]).astype(float)
            #print(np.ones(I)*np.array([elevations]).astype(float))
            rad.elevation['data'] = np.ones(I) * np.array([ele]).astype(float)
            #rad.sweep_mode['data'] = np.array(data['sweep_mode'])
            rad.instrument_name=rad_name
            rad.metadata['instrument_name']=rad_name
            parameters = {}
            frequency = {'data': int(2705000000), 'units': 'Hertz'}
            parameters['frequency'] = frequency
            rad.instrument_parameters=parameters
            
            # Lets Have Fun With Ranges
            numgates=(len(data.dimensions['Gate']))
            gatewidth=np.array(data['GateWidth'])[0]
            start=data.RangeToFirstGate
            ranges=np.ones(numgates)+start
            for i in range(1,len(ranges)):
                ranges[i]=ranges[i-1]+gatewidth
            
            range_dict=get_metadata('Range')
            range_dict['data']=ranges
            range_dict['Units'] = 'Meters'
            rad.range=range_dict
            
            rad.init_gate_altitude()
            rad.init_gate_longitude_latitude()
            
            # Now get other data - ie Radar Fields
            # Reflectivity
            ref_dict = get_metadata('DHCA')
            ref_dict['data'] = np.array(data['DHCA'])
            #ref_dict['data'][ref_dict['data']<-10000] = np.nan
            ref_dict['units']='unitless'
            rad.add_field('DHCA', ref_dict)
                            
            # KDP
            os.chdir(base_path+'/DKDP/'+ele)
            l2=np.sort([f for f in os.listdir('.') if re.search('.netcdf', f)])
            file_path=l2[time_idx]
            #print(file_path)
            with netCDF4.Dataset(file_path,mode='r') as data:
                zdr_dict = get_metadata('KDP')
                zdr_dict['data'] = np.array(data['DKDP'])
                zdr_dict['units'] = 'deg/km'
                rad.add_field('KDP',zdr_dict)
    
    return rad
    
def read_echotop(idz,rad_name, echo_thresh,base_path,eles):
    file_path=idz
    with netCDF4.Dataset(file_path,mode='r') as data:
            #print(data)
            Ref=np.array(data['TraditionalEchoTop_'+echo_thresh[0]])
            I,J=Ref.shape
            rad = pyart.testing.make_empty_ppi_radar(J, I, 1)
            time_dict=get_metadata('Time')
            time_dict['data']=np.array([data.Time])
            time_dict['units']='seconds since 1970-01-01T00:00:00'
            time_dict['calendar']='gregorian'
            rad.time = time_dict
            rad.latitude['data'] = np.array([data.Latitude])
            rad.longitude['data'] = np.array([data.Longitude])
            rad.sweep_number['data'] = np.array([1])
            rad.azimuth['data'] = np.array(data['Azimuth'])
            rad.fixed_angle['data']=np.array([data.Elevation])
            rad.elevation['data'] = np.array([data.Elevation])
            rad.instrument_name=rad_name
            rad.metadata['instrument_name']= rad_name
            parameters = {}
            frequency = {'data': int(2705000000), 'units': 'Hertz'}
            parameters['frequency'] = frequency
            rad.instrument_parameters=parameters
            numgates=(len(data.dimensions['Gate']))
            gatewidth=np.array(data['GateWidth'])[0]
            start=data.RangeToFirstGate
            ranges=np.ones(numgates)+start
            for i in range(1,len(ranges)):
                ranges[i]=ranges[i-1]+gatewidth
            
            range_dict=get_metadata('Range')
            range_dict['data']=ranges
            range_dict['Units'] = 'Meters'
            rad.range=range_dict
            
            rad.init_gate_altitude()
            rad.init_gate_longitude_latitude()
            ref_dict = get_metadata('EchoTop')
            ref_dict['data'] = np.array(data['TraditionalEchoTop_'+echo_thresh[0]])
            ref_dict['data'][ref_dict['data']==-99900] = np.nan
            ref_dict['units']='kilometers'
            rad.add_field('EchoTop_'+echo_thresh[0], ref_dict)
            
    for i in range(1,len(echo_thresh)):
        os.chdir(base_path+'/TraditionalEchoTop_'+echo_thresh[i]+'/'+eles)
        with netCDF4.Dataset(file_path,mode='r') as data:
            ref_dict = get_metadata('EchoTop')
            ref_dict['data'] = np.array(data['TraditionalEchoTop_'+echo_thresh[i]])
            ref_dict['data'][ref_dict['data']==-99900] = np.nan
            ref_dict['units']='kilometers'
            rad.add_field('EchoTop_'+echo_thresh[i], ref_dict)
            
    return rad
def mping_read():
            reqheaders = {
                'content-type':'application/json',
                'Authorization': 'Token 526911e6a3b40b1ed9d74b1b6bd4b73d5f72ef01'
                }
            
            reqparams = {
                'category':'Hail',
                'year': '2024',
                'month': '04',
                'day': '28'
                }
            
            url = 'http://mping.ou.edu/mping/api/v2/reports'
            response = requests.get(url,params=reqparams, headers=reqheaders)
            
            if response.status_code != 200:
                print ('Request Failed with status code %i' )# response.status_code
            else:
                print ('Request Successful')
                #print (response.url)
                data = response.json()
                # Pretty print the data
                #print (json.dumps(data,indent=4))
            
            data_res = data['results']   
            
            times = []
            des = []
            geom = []
            for i in range(0,len(data_res)):
                lists = data_res[i]
                #time = lists['obtime'][11:]
                times = np.append(times,lists['obtime'][11:])
                des = np.append(des,lists['description'])
                geom = np.append(geom,lists['geom'])
            
            lat = []
            lon = []    
            #Now need to split the lat and lons
            for j in range(0,len(geom)):
                coords = geom[j]['coordinates']
                lon = np.append(lon,coords[0])
                lat = np.append(lat,coords[1])
            
            times = np.array(times)
            
            return lat, lon, times, des

def plot_Horus_PPI(radar,xlims,ylims,figpath,dateHrus,Hrus_time,ele,echo_thresh,save_fig,plot_MESH):
#def plot_Horus_PPI(radar,rad_MESH,rad_Echo,xlims,ylims,figpath,dateHrus,Hrus_time,ele,echo_thresh,save_fig,plot_MESH):
    cdict = {'red':   [[0.0,  0.0, 0.0],
                   [0.5,  0.0, 0.0],
                   [0.85,  1.0, 1.0],
                   [1.0,  1.0, 1.0]],
         'green': [[0.0,  0.0, 0.0],
                   [0.25, 0.0, 0.0],
                   [0.5,  0.0, 0.0],
                   [0.85,  1.0, 1.0],
                   [1.0,  0.0, 0.0]],
         'blue':  [[0.0,  1.0, 1.0],
                   [0.5,  1.0, 1.0],
                   [1.0,  0.1, 0.1]]}
    cc_cmap = colors.LinearSegmentedColormap('testCmap', segmentdata=cdict, N=256)
    ticks = np.arange(1,15)
    labs = ['L/MR','HR','R/Ha','BD','AP','BI','UK','NE','DS','WS','CY','','GR','+']
    hca_cmap = colors.ListedColormap(['#AAAAAA','#0000FE','#0057FF','#00A9FF','#00FFFF','#56FFAA','#ACFD56','#FFFF00','#FEAB02','#FF5700','#FF00FF','#B300B3','#DC143C','#FF00FF'])
    hca_cmap = colors.ListedColormap(['limegreen','green','#DC143C','#FFFF00','dimgray','lightgray','#B300B3','white','#00FFFF','#0000FE','mediumturquoise','#000000','#FF00FF','#000000'])

    os.chdir(figpath)
    fig = plt.figure()
    display = pyart.graph.RadarMapDisplay(radar)
    display.plot_ppi('Reflectivity',vmin=-10,vmax=70,cmap='pyart_NWSRef')
    plt.xlim(xlims)
    plt.ylim(ylims)
    if save_fig is True:
        plt.savefig(dateHrus+'_'+str(Hrus_time)+'_'+ele[0:2]+ele[3:]+'_Reflectivity',bbox_inches = 'tight')
    else:
        plt.show(fig)
    plt.close(fig)
    #plt.show(fig)
    
    fig = plt.figure()
    display = pyart.graph.RadarMapDisplay(radar)
    display.plot_ppi('corrected_velocity',vmin=-50, vmax=50,cmap = 'pyart_BuDRd18')
    plt.xlim(xlims)
    plt.ylim(ylims)
    if save_fig is True:
        plt.savefig(dateHrus+'_'+str(Hrus_time)+'_'+ele[0:2]+ele[3:]+'_Velocity',bbox_inches = 'tight')
    else:
        plt.show(fig)
    plt.close(fig)
    #plt.show(fig)
    
    fig = plt.figure()
    display.plot_ppi('ZDR',vmin=-1,vmax=7,cmap='pyart_NWSRef')
    plt.xlim(xlims)
    plt.ylim(ylims)
    if save_fig is True:
        plt.savefig(dateHrus+'_'+str(Hrus_time)+'_'+ele[0:2]+ele[3:]+'_ZDR',bbox_inches = 'tight')
    else:
        plt.show(fig)
    plt.close(fig)
    #plt.show(fig)
    
    fig = plt.figure()
    display.plot_ppi('RhoHV',vmin=0.5,vmax=1,cmap = cc_cmap)
    plt.xlim(xlims)
    plt.ylim(ylims)
    if save_fig is True:
        plt.savefig(dateHrus+'_'+str(Hrus_time)+'_'+ele[0:2]+ele[3:]+'_RhoHV',bbox_inches = 'tight')
    else:
        plt.show(fig)
    plt.close(fig)
    #plt.show(fig)
    
    fig = plt.figure()
    display.plot_ppi('PhiDP',vmin=-180,vmax=180,cmap = 'pyart_RRate11')
    plt.xlim(xlims)
    plt.ylim(ylims)
    if save_fig is True:
        plt.savefig(dateHrus+'_'+str(Hrus_time)+'_'+ele[0:2]+ele[3:]+'_PhiDP',bbox_inches = 'tight')
    else:
        plt.show(fig)
    plt.close(fig)
    '''
    fig = plt.figure()
    display.plot_ppi('KDP',vmin=-1,vmax=10,cmap = 'pyart_NWSRef')
    plt.xlim(xlims)
    plt.ylim(ylims)
    if save_fig is True:
        plt.savefig(dateHrus+'_'+str(Hrus_time)+'_'+ele[0:2]+ele[3:]+'_KDP',bbox_inches = 'tight')
    else:
        plt.show(fig)
    plt.close(fig)'''
    
    # if plot_MESH is True:
    #     fig = plt.figure()
    #     #ax = fig.add_subplot(111)
    #     display2 = pyart.graph.RadarMapDisplay(rad_MESH)
    #     display.plot_ppi('Reflectivity',vmin=-10,vmax=70,cmap='pyart_NWSRef',alpha=0.1,colorbar_flag=False)
    #     display2.plot_ppi('MESH',cmap='pyart_NWSRef',vmin=0,vmax=100)
    #     plt.xlim(xlims)
    #     plt.ylim(ylims)
    #     if save_fig is True:
    #         plt.savefig(dateHrus+'_'+str(Hrus_time)+'_'+ele[0:2]+ele[3:]+'_MESHRef',bbox_inches = 'tight')
    #     else:
    #         plt.show(fig)
    #     plt.close(fig)
    #     #plt.show(fig)
        
    #     fig = plt.figure()
    #     #ax = fig.add_subplot(111)
    #     display2 = pyart.graph.RadarMapDisplay(rad_MESH)
    #     #display.plot_ppi('Reflectivity',vmin=-10,vmax=70,cmap='pyart_NWSRef',alpha=0.1,colorbar_flag=False)
    #     display2.plot_ppi('MESH',cmap='pyart_NWSRef',vmin=0,vmax=100)
    #     plt.xlim(xlims)
    #     plt.ylim(ylims)
    #     if save_fig is True:
    #         plt.savefig(dateHrus+'_'+str(Hrus_time)+'_MESH',bbox_inches = 'tight')
    #     else:
    #         plt.show(fig)
    #     plt.close(fig)
    
    # for j in range(0,len(echo_thresh)):
    #     fig = plt.figure()
    #     display3 = pyart.graph.RadarMapDisplay(rad_Echo)
    #     display3.plot_ppi('EchoTop_'+echo_thresh[j],cmap='pyart_NWSRef',vmin=0,vmax=20)
    #     plt.xlim(xlims)
    #     plt.ylim(ylims)
    #     if save_fig is True:
    #         plt.savefig(dateHrus+'_'+str(Hrus_time)+'_EchoTop_'+echo_thresh[j],bbox_inches = 'tight')
    #     else:
    #         plt.show(fig)
    #     plt.close(fig)
    # # Sets up the HCA colormaps to be used later
    # ticks = np.arange(1,15)
    # labs = ['L/MR','HR','R/Ha','BD','AP','BI','UK','NE','DS','WS','CY','','GR','+']
    # hca_cmap = colors.ListedColormap(['#AAAAAA','#0000FE','#0057FF','#00A9FF','#00FFFF','#56FFAA','#ACFD56','#FFFF00','#FEAB02','#FF5700','#FF00FF','#B300B3','#DC143C','#FF00FF'])
    # hca_cmap = colors.ListedColormap(['limegreen','green','#DC143C','#FFFF00','dimgray','lightgray','#B300B3','white','#00FFFF','#0000FE','mediumturquoise','#000000','#FF00FF','#000000'])
    
    # fig = plt.figure()
    # #ax = fig.add_subplot(111)
    # #display2 = pyart.graph.RadarMapDisplay(rad_MESH)
    # display.plot_ppi('HCA',vmin=0.5,vmax=14.5,cmap=hca_cmap,ticklabs=labs,ticks=ticks)
    # #display2.plot_ppi('MESH',cmap='pyart_NWSRef',vmin=0,vmax=10)
    # plt.xlim(xlims)
    # plt.ylim(ylims)
    # if save_fig is True:
    #     plt.savefig(dateHrus+'_'+str(Hrus_time)+'_'+ele[0:2]+ele[3:]+'_DHCA',bbox_inches = 'tight')
    # else:
    #     plt.show(fig)
    # plt.close(fig)
    
def Time_Finder(l_KTLX,Horusfile,match_wsr):
    if match_wsr is not True:
        Hrs_time_str = Horusfile[9:15]
        Hrs_time = int(Hrs_time_str)
    else:
        Hrs_time_str = Horusfile[13:19]
        Hrs_time = int(Hrs_time_str)
    break_flag = 0
    #print(Hrs_time_str)
    for m in range(0,len(l_KTLX)):
        #print(m)
        if break_flag == 1:
            break
        if match_wsr is not True:
            KTLX_time = l_KTLX[m][13:19]
            #print(KTLX_time)
            #print(Hrs_time_str)
            #print(int(Hrs_time_str[0]))
            if int(KTLX_time) < 100000:
                if int(Hrs_time_str[0]) == 0:
                    time_diff_idx = (int(Hrs_time)+240000) - (int(KTLX_time)+240000)
                else:
                    time_diff_idx = int(Hrs_time) - (int(KTLX_time)+240000)
            else:
                if int(Hrs_time_str[0]) == 0:
                    time_diff_idx = (int(Hrs_time)+240000) - int(KTLX_time)
                else: 
                    time_diff_idx = int(Hrs_time) - int(KTLX_time)
            #print(time_diff_idx)
        else:
            #print('Down Here')
            KTLX_time = l_KTLX[m][9:15]
            #print(KTLX_time)
            if int(KTLX_time) < 100000:
                #print('In Here')
                if int(Hrs_time_str[0]) == 0:
                    time_diff_idx = (int(Hrs_time)+240000) - (int(KTLX_time)+240000)
                else:
                    time_diff_idx = int(Hrs_time) - (int(KTLX_time)+240000)
            else:
                #print('No Im Here')
                if int(Hrs_time_str[0]) == 0:
                    #print('Hello')
                    time_diff_idx = (int(Hrs_time)+240000) - int(KTLX_time)
                else: 
                    #print('Hi')
                    #print(KTLX_time)
                    #print(Hrs_time)
                    time_diff_idx = int(Hrs_time) - int(KTLX_time)
                
            
        if time_diff_idx <= 0:
            break_flag=1
            time_idx = m - 1
            #print(KTLX_time)
            #print(l_KTLX[time_idx])
    return time_idx

def Time_Finder_Elevation(l_KTLX,Horusfile,match_wsr):
    if match_wsr is not True:
        Hrs_time_str = Horusfile[9:15]
        Hrs_time = int(Hrs_time_str)
    else:
        Hrs_time_str = Horusfile[13:19]
        Hrs_time = int(Hrs_time_str)
    break_flag = 0
    #print(Hrs_time_str)
    for m in range(0,len(l_KTLX)):
        #print(m)
        if break_flag == 1:
            break
        if match_wsr is not True:
            KTLX_time = l_KTLX[m][9:15]
            #print(KTLX_time)
            #print(Hrs_time_str)
            #print(int(Hrs_time_str[0]))
            if int(KTLX_time) < 100000:
                if int(Hrs_time_str[0]) == 0:
                    time_diff_idx = (int(Hrs_time)+240000) - (int(KTLX_time)+240000)
                else:
                    time_diff_idx = int(Hrs_time) - (int(KTLX_time)+240000)
            else:
                if int(Hrs_time_str[0]) == 0:
                    time_diff_idx = (int(Hrs_time)+240000) - int(KTLX_time)
                else: 
                    time_diff_idx = int(Hrs_time) - int(KTLX_time)
            #print(time_diff_idx)
        else:
            #print('Down Here')
            KTLX_time = l_KTLX[m][9:15]
            #print(KTLX_time)
            if int(KTLX_time) < 100000:
                #print('In Here')
                if int(Hrs_time_str[0]) == 0:
                    time_diff_idx = (int(Hrs_time)+240000) - (int(KTLX_time)+240000)
                else:
                    time_diff_idx = int(Hrs_time) - (int(KTLX_time)+240000)
            else:
                #print('No Im Here')
                if int(Hrs_time_str[0]) == 0:
                    #print('Hello')
                    time_diff_idx = (int(Hrs_time)+240000) - int(KTLX_time)
                else: 
                    #print('Hi')
                    #print(KTLX_time)
                    #print(Hrs_time)
                    time_diff_idx = int(Hrs_time) - int(KTLX_time)
                
            
        if time_diff_idx <= 0:
            break_flag=1
            time_idx = m - 1
            #print(KTLX_time)
            #print(l_KTLX[time_idx])
    return time_idx
def Comparison_Plots(Horusrad,KTLX_rad,KTLXalgs,lon,lat,figpath,dateHrus,Hrus_time,cbar_ax,sweep_num,sweep_Hrus,save_fig):
    os.chdir(figpath)
    
    x,y,z_Hrus = Horusrad.get_gate_x_y_z(sweep_Hrus)
    # This isn't perfect but will work for the moment
    a,b = Horusrad.fields['Reflectivity']['data'].shape
    #print(a)
    a2,b2 = z_Hrus.shape
    #print(a2)
    z_Hrus_2 = np.matlib.repmat(z_Hrus,int(a/a2),1)
    
    ref_dict = get_metadata('Heights')
    ref_dict['data'] = z_Hrus_2/1000
    Horusrad.add_field('Heights', ref_dict)
    
    # This one needs to be done differently bc of something odd
    #z_WSR_2 = []
    for n in range(0,KTLX_rad.nsweeps):
        x,y,z_WSR = KTLX_rad.get_gate_x_y_z(n)
        if n == 0:
            z_WSR_2 = z_WSR/1000
        else:
            z_WSR_2 = np.concatenate((z_WSR_2,z_WSR/1000),axis=0)
    #a,b = KTLX_rad.fields['reflectivity']['data'].shape
    #print(a,b)
    #a2,b2 = z_WSR.shape
    #print(a2,b2)
    #z_WSR_2 = np.matlib.repmat(z_WSR,int(a/a2),1)
    print(z_WSR_2.shape)
    ref_dict = get_metadata('Heights')
    ref_dict['data'] = z_WSR_2
    KTLX_rad.add_field('Heights', ref_dict)
    
    
    # First, lets plot the heights to see if we are in line with each other 
    fig = plt.figure(figsize=(16,12))
    display = pyart.graph.RadarMapDisplay(KTLX_rad)
    display2 = pyart.graph.RadarMapDisplay(Horusrad)
    ax1 = plt.subplot(121,projection=ccrs.PlateCarree())
    display.plot_ppi_map('Heights',sweep=sweep_num,vmin=0,vmax=30,ax=ax1,cmap='pyart_NWSRef',
                         min_lon = lon[0],max_lon=lon[1],
                         min_lat = lat[0], max_lat=lat[1],
                         colorbar_flag=False)
    #plt.gca().xaxis.set_major_locator(plt.NullLocator())
    
    ax2 = plt.subplot(122,projection=ccrs.PlateCarree())
    display2.plot_ppi_map('Heights',sweep=sweep_Hrus,vmin=0,vmax=30,ax=ax2,cmap='pyart_NWSRef',
                         min_lon=lon[0],max_lon=lon[1],
                         min_lat=lat[0],max_lat=lat[1],
                         colorbar_flag=False)
    
    cax = ax1.inset_axes([cbar_ax[0], cbar_ax[1], 0.89, 0.03],transform=ax1.transData)
    m = cm.ScalarMappable(cmap='pyart_NWSRef')
    m.set_clim(0,30)
    fig.colorbar(m,cax=cax,orientation='horizontal')
    if save_fig:
        plt.savefig(dateHrus+'_'+str(Hrus_time)+'_Height_Comparison',bbox_inches = 'tight')
    else:
        plt.show()
    plt.close()
    
    ### Reflectivity
    fig = plt.figure(figsize=(16,12))
    display = pyart.graph.RadarMapDisplay(KTLX_rad)
    display2 = pyart.graph.RadarMapDisplay(Horusrad)
    ax1 = plt.subplot(121,projection=ccrs.PlateCarree())
    display.plot_ppi_map('reflectivity',sweep=sweep_num,ax=ax1,vmin=-10,vmax=70,cmap='pyart_NWSRef',
                         min_lon = lon[0],max_lon=lon[1],
                         min_lat = lat[0], max_lat=lat[1],
                         colorbar_flag=False)
    #plt.gca().xaxis.set_major_locator(plt.NullLocator())
    
    ax2 = plt.subplot(122,projection=ccrs.PlateCarree())
    display2.plot_ppi_map('Reflectivity',sweep=sweep_Hrus,ax=ax2,vmin=-10,vmax=70,cmap='pyart_NWSRef',
                         min_lon=lon[0],max_lon=lon[1],
                         min_lat=lat[0],max_lat=lat[1],
                         colorbar_flag=False)
    
    cax = ax1.inset_axes([cbar_ax[0], cbar_ax[1], 0.89, 0.03],transform=ax1.transData)
    m = cm.ScalarMappable(cmap='pyart_NWSRef')
    m.set_clim(-10,70)
    fig.colorbar(m,cax=cax,orientation='horizontal')
    if save_fig:
        plt.savefig(dateHrus+'_'+str(Hrus_time)+'_Reflectivity_Comparison',bbox_inches = 'tight')
    else:
        plt.show()
    plt.close()
    #del ax1,ax2
    
    # ### ZDR
    # fig = plt.figure(figsize=(16,12))
    # display = pyart.graph.RadarMapDisplay(KTLX_rad)
    # display2 = pyart.graph.RadarMapDisplay(Horusrad)
    # ax1 = plt.subplot(121,projection=ccrs.PlateCarree())
    # display.plot_ppi_map('differential_reflectivity',sweep=sweep_num,ax=ax1,vmin=-1,vmax=7,cmap='pyart_NWSRef',
    #                      min_lon = lon[0],max_lon=lon[1],
    #                      min_lat = lat[0], max_lat=lat[1],
    #                      colorbar_flag=False)
    # #plt.gca().xaxis.set_major_locator(plt.NullLocator())
    
    # ax2 = plt.subplot(122,projection=ccrs.PlateCarree())
    # display2.plot_ppi_map('ZDR',sweep=sweep_Hrus,ax=ax2,vmin=-1,vmax=7,cmap='pyart_NWSRef',
    #                      min_lon=lon[0],max_lon=lon[1],
    #                      min_lat=lat[0],max_lat=lat[1],
    #                      colorbar_flag=False)
    
    # #plt.show(fig)
    # cax = ax1.inset_axes([cbar_ax[0], cbar_ax[1], 0.89, 0.03],transform=ax1.transData)
    # m = cm.ScalarMappable(cmap='pyart_NWSRef')
    # m.set_clim(-1,7)
    # fig.colorbar(m,cax=cax,orientation='horizontal')
    # if save_fig:
    #     plt.savefig(dateHrus+'_'+str(Hrus_time)+'_ZDR_Comparison',bbox_inches = 'tight')
    # else:
    #     plt.show()
    # plt.close()
    
    # ### Correlation Coefficient
    # fig = plt.figure(figsize=(16,12))
    # display = pyart.graph.RadarMapDisplay(KTLX_rad)
    # display2 = pyart.graph.RadarMapDisplay(Horusrad)
    # ax1 = plt.subplot(121,projection=ccrs.PlateCarree())
    # display.plot_ppi_map('cross_correlation_ratio',sweep=sweep_num,ax=ax1,vmin=0.25,vmax=1,cmap=cc_cmap,
    #                      min_lon = lon[0],max_lon=lon[1],
    #                      min_lat = lat[0], max_lat=lat[1],
    #                      colorbar_flag=False)
    # #plt.gca().xaxis.set_major_locator(plt.NullLocator())
    
    # ax2 = plt.subplot(122,projection=ccrs.PlateCarree())
    # display2.plot_ppi_map('RhoHV',sweep=sweep_Hrus,ax=ax2,vmin=0.25,vmax=1,cmap=cc_cmap,
    #                      min_lon=lon[0],max_lon=lon[1],
    #                      min_lat=lat[0],max_lat=lat[1],
    #                      colorbar_flag=False)
    
    # #plt.show(fig)
    # cax = ax1.inset_axes([cbar_ax[0], cbar_ax[1], 0.89, 0.03],transform=ax1.transData)
    # m = cm.ScalarMappable(cmap=cc_cmap)
    # m.set_clim(0.25,1)
    # fig.colorbar(m,cax=cax,orientation='horizontal')
    # if save_fig:
    #     plt.savefig(dateHrus+'_'+str(Hrus_time)+'_RhoHV_Comparison',bbox_inches = 'tight')
    # else:
    #     plt.show()
    # plt.close()
    
    # ### PhiDP
    # fig = plt.figure(figsize=(16,12))
    # display = pyart.graph.RadarMapDisplay(KTLX_rad)
    # display2 = pyart.graph.RadarMapDisplay(Horusrad)
    # ax1 = plt.subplot(121,projection=ccrs.PlateCarree())
    # display.plot_ppi_map('differential_phase',sweep=sweep_num,ax=ax1,vmin=-180,vmax=180,cmap='pyart_RRate11',
    #                      min_lon = lon[0],max_lon=lon[1],
    #                      min_lat = lat[0], max_lat=lat[1],
    #                      colorbar_flag=False)
    # #plt.gca().xaxis.set_major_locator(plt.NullLocator())
    
    # ax2 = plt.subplot(122,projection=ccrs.PlateCarree())
    # display2.plot_ppi_map('PhiDP',sweep=sweep_Hrus,ax=ax2,vmin=-180,vmax=180,cmap='pyart_RRate11',
    #                      min_lon=lon[0],max_lon=lon[1],
    #                      min_lat=lat[0],max_lat=lat[1],
    #                      colorbar_flag=False)
    
    # #plt.show(fig)
    # cax = ax1.inset_axes([cbar_ax[0], cbar_ax[1], 0.89, 0.03],transform=ax1.transData)
    # m = cm.ScalarMappable(cmap='pyart_RRate11')
    # m.set_clim(-180,180)
    # fig.colorbar(m,cax=cax,orientation='horizontal')
    # if save_fig:
    #     plt.savefig(dateHrus+'_'+str(Hrus_time)+'_PhiDP_Comparison',bbox_inches = 'tight')
    # else:
    #     plt.show()
    # plt.close()
    
    ### Velocity
    fig = plt.figure(figsize=(16,12))
    display = pyart.graph.RadarMapDisplay(KTLX_rad)
    display2 = pyart.graph.RadarMapDisplay(Horusrad)
    ax1 = plt.subplot(121,projection=ccrs.PlateCarree())
    display.plot_ppi_map('velocity',sweep=sweep_num, ax = ax1, vmin=-30, vmax=30,cmap = 'pyart_BuDRd18',
                         min_lon = lon[0],max_lon=lon[1],
                         min_lat = lat[0], max_lat=lat[1],
                         colorbar_flag=False)
    ax2 = plt.subplot(122,projection=ccrs.PlateCarree())
    display2.plot_ppi_map('velocity',sweep=sweep_Hrus,ax=ax2,vmin=-30,vmax=30,cmap='pyart_BuDRd18',
                         min_lon=lon[0],max_lon=lon[1],
                         min_lat=lat[0],max_lat=lat[1],
                         colorbar_flag=False)
    cax = ax1.inset_axes([cbar_ax[0], cbar_ax[1], 0.89, 0.03],transform=ax1.transData)
    m = cm.ScalarMappable(cmap='pyart_BuDRd18')
    m.set_clim(-60,60)
    fig.colorbar(m,cax=cax,orientation='horizontal')
    if save_fig:
        plt.savefig(dateHrus+'_'+str(Hrus_time)+'_Velocity_Comparison',bbox_inches = 'tight')
    else:
        plt.show()
    plt.close()
    
    # ### DHCA
    # # We first need to find the right sweep_number for the algorithm data 
    # idx = KTLX_rad.sweep_start_ray_index['data'][sweep_num]
    # ele = KTLX_rad.elevation['data'][idx].astype(float)
    # eles = KTLXalgs.fixed_angle['data'].astype(float)
    # diff = np.absolute(eles - ele)
    # sweep_alg = diff.argmin()
    # #print(ele)
    # #print(eles)
    # #print(sweep_alg)
    # fig = plt.figure(figsize=(16,12))
    # display = pyart.graph.RadarMapDisplay(KTLXalgs)
    # display2 = pyart.graph.RadarMapDisplay(Horusrad)
    # ax1 = plt.subplot(121,projection=ccrs.PlateCarree())
    # display.plot_ppi_map('DHCA',sweep=sweep_alg, ax = ax1,vmin=0.5,vmax=14.5,cmap=hca_cmap, ticklabs=labs,ticks=ticks,
    #                      min_lon = lon[0],max_lon=lon[1],
    #                      min_lat = lat[0], max_lat=lat[1],
    #                      colorbar_flag=False)
    # ax2 = plt.subplot(122,projection=ccrs.PlateCarree())
    # display2.plot_ppi_map('HCA',sweep=sweep_Hrus,ax=ax2,vmin=0.5,vmax=14.5,cmap=hca_cmap,ticklabs=labs,ticks=ticks,
    #                      min_lon=lon[0],max_lon=lon[1],
    #                      min_lat=lat[0],max_lat=lat[1],
    #                      colorbar_flag=False)
    # cax = ax1.inset_axes([cbar_ax[0], cbar_ax[1], 0.89, 0.03],transform=ax1.transData)
    # m = cm.ScalarMappable(cmap=hca_cmap)
    # m.set_clim(0.5,14.5)
    # fig.colorbar(m,cax=cax,orientation='horizontal')
    # if save_fig:
    #     plt.savefig(dateHrus+'_'+str(Hrus_time)+'_HCA_Comparison',bbox_inches = 'tight')
    # else:
    #     plt.show()
    # plt.close()
    
    # ### MESH
    # fig2 = plt.figure(figsize=(16,12))
    # display = pyart.graph.RadarMapDisplay(KTLXmesh)
    # display2 = pyart.graph.RadarMapDisplay(Horusmesh)
    # ax1 = plt.subplot(121,projection=ccrs.PlateCarree())
    # display.plot_ppi_map('MESH',sweep=0,ax=ax1,vmin=0,vmax=80,cmap='pyart_NWSRef',
    #                      min_lon = lon[0],max_lon=lon[1],
    #                      min_lat = lat[0], max_lat=lat[1],
    #                      colorbar_flag=False)
    # #plt.gca().xaxis.set_major_locator(plt.NullLocator())
    
    # ax2 = plt.subplot(122,projection=ccrs.PlateCarree())
    # display2.plot_ppi_map('MESH',sweep=0,ax=ax2,vmin=0,vmax=80,cmap='pyart_NWSRef',
    #                      min_lon=lon[0],max_lon=lon[1],
    #                      min_lat=lat[0],max_lat=lat[1],
    #                      colorbar_flag=False)
    # #del display
    # #plt.show(fig)
    # cax = ax1.inset_axes([cbar_ax[0], cbar_ax[1], 0.89, 0.03],transform=ax1.transData)
    # m2 = cm.ScalarMappable(cmap='pyart_NWSRef')
    # m2.set_clim(0,80)
    # fig2.colorbar(m2,cax=cax,orientation='horizontal')
    # if save_fig:
    #     plt.savefig(dateHrus+'_'+str(Hrus_time)+'_MESH_Comparison',bbox_inches = 'tight')
    # else:
    #     plt.show()
    # plt.close()
    
def Gridding(Horusrad,KTLX_rad,levels_Horus,grid_lims,lon,lat,grid_x_y,grid_z,grid_level,weighting,roi_func,min_radius,grid_algo,cbar_ax,save_fig):
    # Set up the Gridded Radars for both Horus and our comparison radar (in this case KTLX)
    grid = pyart.map.grid_from_radars(
            (KTLX_rad),
            gridding_algo = grid_algo,
            grid_shape=(grid_z, grid_x_y, grid_x_y),
            grid_limits=grid_lims,
            grid_origin=(Horusrad.latitude['data'][0], Horusrad.longitude['data'][0]),
            weighting_function=weighting,
            roi_func = roi_func,
            min_radius = min_radius,
            fields=["reflectivity","velocity","differential_reflectivity","cross_correlation_ratio","differential_phase"],
        )
    #print(grid)
    
    grid2 = pyart.map.grid_from_radars(
            (Horusrad),
            gridding_algo = grid_algo,
            grid_shape=(grid_z, grid_x_y, grid_x_y),
            grid_limits=grid_lims,
            grid_origin=(Horusrad.latitude['data'][0], Horusrad.longitude['data'][0]),
            weighting_function = weighting,
            roi_func = roi_func,
            min_radius = min_radius,
            fields=["Reflectivity","velocity","ZDR","RhoHV","PhiDP"],
        )
    
    # Now Take these Grids and Plot Them
    
    ## Reflectivity
    fig = plt.figure(figsize=(18, 6))
    ax1 = plt.subplot(121, projection=ccrs.PlateCarree())
    display = pyart.graph.GridMapDisplay(grid)
    pcm = display.plot_grid(
        "reflectivity",
        level = grid_level,
        ax=ax1,
        cmap="pyart_NWSRef",
        vmin=-10,
        vmax=70,
        colorbar_flag=False
        )
    plt.xlim(lon)
    plt.ylim(lat)
    
    ax2 = plt.subplot(122,projection=ccrs.PlateCarree())
    display2 = pyart.graph.GridMapDisplay(grid2)
    display2.plot_grid(
        'Reflectivity',
        level = grid_level,
        ax = ax2, 
        cmap='pyart_NWSRef',
        vmin=-10,
        vmax=70,
        colorbar_flag=False
        )
    plt.xlim(lon)
    plt.ylim(lat)
    
    # This will put the colorbar back on
    cax = ax1.inset_axes([cbar_ax[0], cbar_ax[1], 1.1, 0.05],transform=ax1.transData)
    m = cm.ScalarMappable(cmap='pyart_NWSRef')
    m.set_clim(-10,70)
    fig.colorbar(m,cax=cax,orientation='horizontal')
    if save_fig:
        print('Laura Fix This Still')
    else:
        plt.show(fig)
        
    ## ZDR
    fig = plt.figure(figsize=(18, 6))
    ax1 = plt.subplot(121, projection=ccrs.PlateCarree())
    display = pyart.graph.GridMapDisplay(grid)
    pcm = display.plot_grid(
        "differential_reflectivity",
        level = grid_level,
        ax=ax1,
        cmap="pyart_NWSRef",
        vmin=-1,
        vmax=7,
        colorbar_flag=False
        )
    plt.xlim(lon)
    plt.ylim(lat)
    
    ax2 = plt.subplot(122,projection=ccrs.PlateCarree())
    display2 = pyart.graph.GridMapDisplay(grid2)
    display2.plot_grid(
        'ZDR',
        level = grid_level,
        ax = ax2, 
        cmap='pyart_NWSRef',
        vmin=-1,
        vmax=7,
        colorbar_flag=False
        )
    plt.xlim(lon)
    plt.ylim(lat)
    
    # This will put the colorbar back on
    cax = ax1.inset_axes([cbar_ax[0], cbar_ax[1], 1.1, 0.05],transform=ax1.transData)
    m = cm.ScalarMappable(cmap='pyart_NWSRef')
    m.set_clim(-1,7)
    fig.colorbar(m,cax=cax,orientation='horizontal')
    if save_fig:
        print('Laura Fix This Still')
    else:
        plt.show(fig)
    
    ## Correlation Coefficient
    fig = plt.figure(figsize=(18, 6))
    ax1 = plt.subplot(121, projection=ccrs.PlateCarree())
    display = pyart.graph.GridMapDisplay(grid)
    pcm = display.plot_grid(
        "cross_correlation_ratio",
        level = grid_level,
        ax=ax1,
        cmap=cc_cmap,
        vmin=0.25,
        vmax=1,
        colorbar_flag=False
        )
    plt.xlim(lon)
    plt.ylim(lat)
    
    ax2 = plt.subplot(122,projection=ccrs.PlateCarree())
    display2 = pyart.graph.GridMapDisplay(grid2)
    display2.plot_grid(
        'RhoHV',
        level = grid_level,
        ax = ax2, 
        cmap=cc_cmap,
        vmin=0.25,
        vmax=1,
        colorbar_flag=False
        )
    plt.xlim(lon)
    plt.ylim(lat)
    
    # This will put the colorbar back on
    cax = ax1.inset_axes([cbar_ax[0], cbar_ax[1], 1.1, 0.05],transform=ax1.transData)
    m = cm.ScalarMappable(cmap=cc_cmap)
    m.set_clim(0.25,1)
    fig.colorbar(m,cax=cax,orientation='horizontal')
    if save_fig:
        print('Laura Fix This Still')
    else:
        plt.show(fig)
    
    ## Correlation Coefficient
    fig = plt.figure(figsize=(18, 6))
    ax1 = plt.subplot(121, projection=ccrs.PlateCarree())
    display = pyart.graph.GridMapDisplay(grid)
    pcm = display.plot_grid(
        "differential_phase",
        level = grid_level,
        ax=ax1,
        cmap='pyart_RRate11',
        vmin=-180,
        vmax=180,
        colorbar_flag=False
        )
    plt.xlim(lon)
    plt.ylim(lat)
    
    ax2 = plt.subplot(122,projection=ccrs.PlateCarree())
    display2 = pyart.graph.GridMapDisplay(grid2)
    display2.plot_grid(
        'PhiDP',
        level = grid_level,
        ax = ax2, 
        cmap='pyart_RRate11',
        vmin=-180,
        vmax=180,
        colorbar_flag=False
        )
    plt.xlim(lon)
    plt.ylim(lat)
    
    # This will put the colorbar back on
    cax = ax1.inset_axes([cbar_ax[0], cbar_ax[1], 1.1, 0.05],transform=ax1.transData)
    m = cm.ScalarMappable(cmap=cc_cmap)
    m.set_clim(-180,180)
    fig.colorbar(m,cax=cax,orientation='horizontal')
    if save_fig:
        print('Laura Fix This Still')
    else:
        plt.show(fig)
    
    ## Velocity
    fig = plt.figure(figsize=(18, 6))
    ax1 = plt.subplot(121, projection=ccrs.PlateCarree())
    display = pyart.graph.GridMapDisplay(grid)
        #display_grid2 = pyart.graph.GridMapDisplay(grid2)
        #display_grid3 = pyart.graph.GridMapDisplay(grid3)
    pcm = display.plot_grid(
        "velocity",
        level = grid_level,
        ax=ax1,
        cmap="pyart_BuDRd18",
        vmin=-30,
        vmax=30,
        colorbar_flag=False
        )
    plt.xlim(lon)
    plt.ylim(lat)
    
    ax2 = plt.subplot(122,projection=ccrs.PlateCarree())
    display2 = pyart.graph.GridMapDisplay(grid2)
    display2.plot_grid(
        'velocity',
        level = grid_level,
        ax = ax2, 
        cmap='pyart_BuDRd18',
        vmin=-30,
        vmax=30,
        colorbar_flag=False
        )
    #cbar_ax = fig.add_axes([0.55, -0.15, 0.2, 0.9])
    #display2.plot_colorbar(orientation='horizontal',ax=cbar_ax)
    plt.xlim(lon)
    plt.ylim(lat)
    
    # This stuff should put the colorbar back on (when it actually works)
    cax = ax1.inset_axes([cbar_ax[0], cbar_ax[1], 1.1, 0.05],transform=ax1.transData)
    m = cm.ScalarMappable(cmap='pyart_BuDRd18')
    m.set_clim(-30,30)
    fig.colorbar(m,cax=cax,orientation='horizontal')
    
    #fig.subplots_adjust(right=0.8)
    #cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
    #fig.colorbar(pcm, cax=cbar_ax)
    plt.show(fig)
    
    return grid, grid2

def Texture_Calcs(Horusrad,KTLXrad,comparison,winsize):
    # First, calculate the texture fields
    tex_ref_Horus = pyart.util.texture_along_ray(Horusrad,'Reflectivity',wind_size = winsize)
    tex_ref_dict = get_metadata('Reflectivity Texture')
    tex_ref_dict['data'] = tex_ref_Horus
    Horusrad.add_field('Reflectivity Texture', tex_ref_dict)
    if comparison == True:
        tex_ref_KTLX = pyart.util.texture_along_ray(KTLXrad,'reflectivity')
        KTLXrad.add_field('Reflectivity Texture',tex_ref_KTLX)
    
    # Next, plot the texture fields
    fig = plt.figure()
    ax1 = plt.subplot(111)
    display = pyart.graph.RadarMapDisplay(Horusrad)
    display.plot('Reflectivity Texture',vmin = 0, vmax = 5)
    plt.show()
    
def Horus_RHI(Horusrad,az_plot,date,time,zlim,xlim_r,save_fig,figpath,xlims,ylims):
    os.chdir(figpath)
    xsect = pyart.util.cross_section_ppi(Horusrad,[az_plot])
    points = xsect.get_gate_x_y_z(0)
    #print(points[0])
    #print(points[1])
    #print(points[1][1,:])
    
    # First, plot a PPI of the RHI
    display_ppi = pyart.graph.RadarDisplay(Horusrad)
    fig = plt.figure()
    display_ppi.plot_ppi('Reflectivity',2,vmin=-10,vmax=70,cmap='pyart_NWSRef')
    plt.scatter(points[0][1,:]/1000,points[1][1,:]/1000,s=0.5,c='k')
    plt.xlim(xlims)
    plt.ylim(ylims)
    if save_fig is True:
        plt.savefig('Ref_PPI_Horus_Az'+str(az_plot)+'_'+date+'_'+time,bbox_inches = 'tight')
    else:
        plt.show(fig)
    plt.close(fig)
    
    display_ppi = pyart.graph.RadarDisplay(Horusrad)
    fig = plt.figure()
    display_ppi.plot_ppi('corrected_velocity',2,vmin=-50,vmax=50,cmap='RdBu_r')
    plt.scatter(points[0][1,:]/1000,points[1][1,:]/1000,s=0.5,c='k')
    plt.xlim(xlims)
    plt.ylim(ylims)
    # Define the colorbar from the RadarMapDisplay object
    cbar = display_ppi.cbs[0]
    # Modify the colorbar label and size
    cbar.set_label(label="Radial Velocity (m/s)", fontsize=12)
    # Modify the number of colorbar ticks
    # cbar.set_ticks([-20, 0, 20, 40, 60])
    if save_fig is True:
        plt.savefig('Vel_PPI_Horus_Az'+str(az_plot)+'_'+date+'_'+time,bbox_inches = 'tight')
    else:
        plt.show(fig)
    plt.close(fig)
    
    display_ppi = pyart.graph.RadarDisplay(Horusrad)
    fig = plt.figure()
    display_ppi.plot_ppi('RhoHV',2,vmin=0.25,vmax=1,cmap=cc_cmap)
    plt.scatter(points[0][1,:]/1000,points[1][1,:]/1000,s=0.5,c='k')
    plt.xlim(xlims)
    plt.ylim(ylims)
    if save_fig is True:
        plt.savefig('Rho_PPI_Horus_Az'+str(az_plot)+'_'+date+'_'+time,bbox_inches = 'tight')
    else:
        plt.show(fig)
    plt.close(fig)
    
    display_ppi = pyart.graph.RadarDisplay(Horusrad)
    fig = plt.figure()
    display_ppi.plot_ppi('ZDR',2,vmin=-1,vmax=7,cmap='pyart_NWSRef')
    plt.scatter(points[0][1,:]/1000,points[1][1,:]/1000,s=0.5,c='k')
    plt.xlim(xlims)
    plt.ylim(ylims)
    if save_fig is True:
        plt.savefig('ZDR_PPI_Horus_Az'+str(az_plot)+'_'+date+'_'+time,bbox_inches = 'tight')
    else:
        plt.show(fig)
    plt.close(fig)
    
    display = pyart.graph.RadarDisplay(xsect)
    
    fig = plt.figure()
    display.plot('Reflectivity',vmin=-10,vmax=70,cmap='pyart_NWSRef')
    plt.ylim(zlim)
    plt.xlim(xlim_r)
    if save_fig:
        plt.savefig('Ref_RHI_Horus_Az'+str(az_plot)+'_'+date+'_'+time)
    else:
        plt.show()
    plt.close(fig)
        
    fig = plt.figure()
    display.plot('RhoHV',vmin=0.25,vmax=1,cmap=cc_cmap)
    plt.ylim(zlim)
    plt.xlim(xlim_r)
    if save_fig:
        plt.savefig('Rho_RHI_Horus_Az'+str(az_plot)+'_'+date+'_'+time)
    else:
        plt.show()
    plt.close(fig)
        
    fig = plt.figure()
    display.plot('ZDR',vmin=-1,vmax=7,cmap='pyart_NWSRef')
    plt.ylim(zlim)
    plt.xlim(xlim_r)
    if save_fig:
        plt.savefig('ZDR_RHI_Horus_Az'+str(az_plot)+'_'+date+'_'+time)
    else:
        plt.show()
    plt.close(fig)
    
    fig = plt.figure(figsize = (25,10))
    ax1 = plt.subplot(121)
    display2 = pyart.graph.RadarMapDisplay(Horusrad)
    display2.plot_ppi('ZDR',3,ax=ax1,vmin=-1,vmax=7,cmap='pyart_NWSRef',colorbar_flag=False)
    plt.scatter(points[0][1,:]/1000,points[1][1,:]/1000,s=0.5,c='k')
    plt.xlim(xlims)
    plt.ylim(ylims)
    ax1.tick_params(labelsize=14)
    plt.title(time+' UTC',fontsize=20,fontweight='bold')
    ax2 = plt.subplot(122)
    display.plot('ZDR',ax=ax2,vmin=-1,vmax=7,cmap='pyart_NWSRef',colorbar_flag=False)
    plt.ylim(zlim)
    plt.xlim(xlim_r)
    ax2.tick_params(labelsize=14)
    plt.title(str(az_plot)+' Deg. ',fontsize=20,fontweight='bold')
    #cax = ax1.inset_axes([cbar_ax[0], cbar_ax[1], 0.89, 0.03],transform=ax1.transData)
    #m = cm.ScalarMappable(cmap='pyart_NWSRef')
    #m.set_clim(-1,7)
    if save_fig:
        plt.savefig('ZDR_RHI_PPI_Horus_Az'+str(az_plot)+'_'+date+'_'+time)
    else:
        plt.show()
    plt.close(fig)
    
    fig = plt.figure(figsize = (25,10))
    ax1 = plt.subplot(121)
    display2 = pyart.graph.RadarMapDisplay(Horusrad)
    display2.plot_ppi('corrected_velocity',3,ax=ax1,vmin=-50,vmax=50,cmap='RdBu_r',colorbar_flag=False)
    plt.scatter(points[0][1,:]/1000,points[1][1,:]/1000,s=0.5,c='k')
    plt.xlim(xlims)
    plt.ylim(ylims)
    ax1.tick_params(labelsize=14)
    plt.title(time+' UTC',fontsize=20,fontweight='bold')
    ax2 = plt.subplot(122)
    display.plot('corrected_velocity',ax=ax2,vmin=-50,vmax=50,cmap='RdBu_r',colorbar_flag=False)
    plt.ylim(zlim)
    plt.xlim(xlim_r)
    ax2.tick_params(labelsize=14)
    plt.title(str(az_plot)+' Deg. ',fontsize=20,fontweight='bold')
    #cax = ax1.inset_axes([cbar_ax[0], cbar_ax[1], 0.89, 0.03],transform=ax1.transData)
    #m = cm.ScalarMappable(cmap='pyart_NWSRef')
    #m.set_clim(-1,7)
    if save_fig:
        plt.savefig('Vel_RHI_PPI_Horus_Az'+str(az_plot)+'_'+date+'_'+time)
    else:
        plt.show()
    plt.close(fig)
    
    fig = plt.figure(figsize = (25,10))
    ax1 = plt.subplot(121)
    display2 = pyart.graph.RadarMapDisplay(Horusrad)
    display2.plot_ppi('Reflectivity',3,ax=ax1,vmin=-10,vmax=70,cmap='pyart_NWSRef',colorbar_flag=False)
    plt.scatter(points[0][1,:]/1000,points[1][1,:]/1000,s=0.5,c='k')
    plt.xlim(xlims)
    plt.ylim(ylims)
    ax1.tick_params(labelsize=14)
    plt.title(time+' UTC',fontsize=20,fontweight='bold')
    ax2 = plt.subplot(122)
    display.plot('Reflectivity',ax=ax2,vmin=-10,vmax=70,cmap='pyart_NWSRef',colorbar_flag=False)
    plt.ylim(zlim)
    plt.xlim(xlim_r)
    ax2.tick_params(labelsize=14)
    plt.title(str(az_plot)+' Deg. ',fontsize=20,fontweight='bold')
    #cax = ax1.inset_axes([cbar_ax[0], cbar_ax[1], 0.89, 0.03],transform=ax1.transData)
    #m = cm.ScalarMappable(cmap='pyart_NWSRef')
    #m.set_clim(-1,7)
    if save_fig:
        plt.savefig('Ref_RHI_PPI_Horus_Az'+str(az_plot)+'_'+date+'_'+time)
    else:
        plt.show()
    plt.close(fig)
        
    fig = plt.figure(figsize = (25,10))
    ax1 = plt.subplot(121)
    display.plot('Reflectivity',ax=ax1,vmin=-10,vmax=70,cmap='pyart_NWSRef')
    plt.ylim(zlim)
    plt.xlim(xlim_r)
    ax2 = plt.subplot(122)
    display.plot('corrected_velocity',ax=ax2,vmin=-60,vmax=60,cmap='RdBu_r')
    plt.ylim(zlim)
    plt.xlim(xlim_r)
    if save_fig:
        plt.savefig('VelRef_RHI_Horus_Az'+str(az_plot)+'_'+date+'_'+time)
    else:
        plt.show()
    plt.close(fig)
    
    fig = plt.figure()
    display.plot('corrected_velocity',vmin=-60,vmax=60,cmap='RdBu_r')
    print(f"Colorbars available: {display.cbs}")
    if display.cbs:
        cbar = display.cbs[-1]
        cbar.set_label(label="Radial Velocity (m/s)", fontsize=12)
    else:
        print("No colorbar found!")
    # # Define the colorbar from the RadarMapDisplay object
    # cbar = display.cbs[0]
    # # Modify the colorbar label and size
    # cbar.set_label(label="Radial Velocity (m/s)", fontsize=12)
    plt.ylim(zlim)
    plt.xlim(xlim_r)
    if save_fig:
        plt.savefig('Vel_RHI_Horus_Az'+str(az_plot)+'_'+date+'_'+time)
    else:
        plt.show()
    plt.close(fig)
        
    # fig = plt.figure()
    # display.plot('PhiDP',vmin=-180,vmax=180,cmap='pyart_RRate11')
    # plt.ylim(zlim)
    # plt.xlim(xlim_r)
    # if save_fig:
    #     plt.savefig('PhiDP_RHI_Horus_Az'+str(az_plot)+'_'+date+'_'+time)
    # else:
    #     plt.show()
    # plt.close(fig)
        
    # fig = plt.figure()
    # display.plot('KDP',vmin=-1,vmax=10,cmap='pyart_NWSRef')
    # plt.ylim(zlim)
    # plt.xlim(xlim_r)
    # if save_fig:
    #     plt.savefig('KDP_RHI_Horus_Az'+str(az_plot)+'_'+date+'_'+time)
    # else:
    #     plt.show()
    # plt.close(fig)
        
    # fig = plt.figure()
    # display.plot('HCA',vmin=0.5,vmax=14.5,cmap=hca_cmap,ticklabs=labs,ticks=ticks)
    # plt.ylim(zlim)
    # plt.xlim(xlim_r)
    # if save_fig:
    #     plt.savefig('HCA_RHI_Horus_Az'+str(az_plot)+'_'+date+'_'+time)
    # else:
    #     plt.show()
    # plt.close(fig)
    return xsect
        
def RHI_Comparisons(Horusrad,Horusgrid,KTLXgrid,comp_rad,az,date,time,cbar_ax,figfolder,save_fig):
    os.chdir(figfolder)
    # Horusrad here needs to be an xsect object
    xsect = pyart.util.cross_section_ppi(Horusrad,[az])
    # RHI Comparisons will be done using the gridded radar data 
    lats,lons,alts = xsect.get_gate_lat_lon_alt(0)
    lats = lats[0]
    lons = lons[0]
    Horusloc = np.vstack([lats,lons]).T
    
    start_lat = lats[0]
    start_lon = lons[-1]
    end_lat = lats[-1]
    end_lon = lons[0]
    print(start_lat)
    print(end_lat)
    print(start_lon)
    print(end_lon)
    
    points = xsect.get_gate_x_y_z(0)
    #print(points[0])
    #print(points[1])
    #print(points[1][1,:])
    
    # First, plot a PPI of the RHI
    display_ppi = pyart.graph.RadarDisplay(Horusrad)
    fig = plt.figure()
    display_ppi.plot_ppi('Reflectivity',8,vmin=-10,vmax=70,cmap='pyart_NWSRef')
    plt.scatter(points[0][1,:]/1000,points[1][1,:]/1000,s=0.5,c='k')
    #plt.xlim(xlims)
    #plt.ylim(ylims)
    if save_fig:
        plt.savefig('Ref_PPI_Horus_Az'+str(az)+'_'+date+'_'+time,bbox_inches = 'tight')
    else:
        plt.show()
    #plt.close(fig)
    
    fig = plt.figure()
    display_ppi.plot_ppi('ZDR',8,vmin=-1,vmax=7,cmap='pyart_NWSRef')
    plt.scatter(points[0][1,:]/1000,points[1][1,:]/1000,s=0.5,c='k')
    #plt.xlim(xlims)
    #plt.ylim(ylims)
    if save_fig:
        plt.savefig('Zdr_PPI_Horus_Az'+str(az)+'_'+date+'_'+time,bbox_inches = 'tight')
    else:
        plt.show()
    #plt.close(fig)
            
    # Set up radar displays 
    display_Horus = pyart.graph.GridMapDisplay(Horusgrid)
    display_KTLX = pyart.graph.GridMapDisplay(KTLXgrid)
    ########### 
    # Figures #
    ###########
    fig = plt.figure(figsize=(16,12))
    ax1 = plt.subplot(121, projection=ccrs.PlateCarree())
    display_KTLX.plot_cross_section(
            "reflectivity",
            [start_lat,start_lon],
            [end_lat,end_lon],
            ax = ax1,
            x_axis="lat",
            cmap="pyart_NWSRef",
            vmin=-10,
            vmax=70,
            title=comp_rad +' Reflectivity \n'+ time +'Z',
            colorbar_flag = False
        )
    ax1.set_ylim(0,10)
    ax1.set_aspect(1/20)
    ax2 = plt.subplot(122, projection=ccrs.PlateCarree())
    display_Horus.plot_cross_section(
            "Reflectivity",
            [start_lat,start_lon],
            [end_lat,end_lon],
            ax = ax2,
            x_axis="lat",
            cmap="pyart_NWSRef",
            vmin=-10,
            vmax=70,
            title = 'Horus Reflectivity \n'+ time + 'Z',
            colorbar_flag = False
        )
    ax2.set_ylim(0,10)
    ax2.set_aspect(1/20)
    
    cax = ax1.inset_axes([cbar_ax[1]+0.4, -2, 1.1, 0.6],transform=ax1.transData)
    m = cm.ScalarMappable(cmap='pyart_NWSRef')
    m.set_clim(-10,70)
    fig.colorbar(m,cax=cax,orientation='horizontal')
    
    if save_fig:
        plt.savefig('Reflectivity_Comparison_RHI_'+str(az)+'_'+date+'_'+time)
    else:
        plt.show()
        
    
    fig = plt.figure(figsize=(16,12))
    ax1 = plt.subplot(121, projection=ccrs.PlateCarree())
    display_KTLX.plot_cross_section(
            "velocity",
            [start_lat,start_lon],
            [end_lat,end_lon],
            ax = ax1,
            x_axis="lat",
            cmap="pyart_BuDRd18",
            vmin=-30,
            vmax=30,
            title=comp_rad +' Velocity \n'+ time +'Z',
            colorbar_flag = False
        )
    ax1.set_ylim(0,10)
    ax1.set_aspect(1/20)
    ax2 = plt.subplot(122, projection=ccrs.PlateCarree())
    display_Horus.plot_cross_section(
            "corrected_velocity",
            [start_lat,start_lon],
            [end_lat,end_lon],
            ax = ax2,
            x_axis="lat",
            cmap="pyart_BuDRd18",
            vmin=-30,
            vmax=30,
            title = 'Horus Velocity \n'+ time + 'Z',
            colorbar_flag = False
        )
    ax2.set_ylim(0,10)
    ax2.set_aspect(1/20)
    cax = ax1.inset_axes([cbar_ax[1]+0.4, -2, 1.1, 0.6],transform=ax1.transData)
    m = cm.ScalarMappable(cmap='pyart_BuDRd18')
    m.set_clim(-30,30)
    fig.colorbar(m,cax=cax,orientation='horizontal')
    if save_fig:
        plt.savefig('Velocity_Comparison_RHI_'+str(az)+'_'+date+'_'+time)
    else:
        plt.show()
        
    fig = plt.figure(figsize=(16,12))
    ax1 = plt.subplot(121, projection=ccrs.PlateCarree())
    display_KTLX.plot_cross_section(
            "cross_correlation_ratio",
            [start_lat,start_lon],
            [end_lat,end_lon],
            ax = ax1,
            x_axis="lat",
            cmap=cc_cmap,
            vmin=0.25,
            vmax=1,
            title=comp_rad +' RhoHV \n'+ time +'Z',
            colorbar_flag = False
        )
    ax1.set_ylim(0,10)
    ax1.set_aspect(1/20)
    ax2 = plt.subplot(122, projection=ccrs.PlateCarree())
    display_Horus.plot_cross_section(
            "RhoHV",
            [start_lat,start_lon],
            [end_lat,end_lon],
            ax = ax2,
            x_axis="lat",
            cmap=cc_cmap,
            vmin=0.25,
            vmax=1,
            title = 'Horus RhoHV \n'+ time + 'Z',
            colorbar_flag = False
        )
    ax2.set_ylim(0,10)
    ax2.set_aspect(1/20)
    cax = ax1.inset_axes([cbar_ax[1]+0.4, -2, 1.1, 0.6],transform=ax1.transData)
    m = cm.ScalarMappable(cmap=cc_cmap)
    m.set_clim(0.25,1)
    fig.colorbar(m,cax=cax,orientation='horizontal')
    if save_fig:
        plt.savefig('RhoHV_Comparison_RHI_'+str(az)+'_'+date+'_'+time)
    else:
        plt.show()
        
    fig = plt.figure(figsize=(16,12))
    ax1 = plt.subplot(121, projection=ccrs.PlateCarree())
    display_KTLX.plot_cross_section(
            "differential_reflectivity",
            [start_lat,start_lon],
            [end_lat,end_lon],
            ax = ax1,
            x_axis="lat",
            cmap="pyart_NWSRef",
            vmin=-1,
            vmax=7,
            title=comp_rad +' ZDR \n'+ time +'Z',
            colorbar_flag = False
        )
    ax1.set_ylim(0,10)
    ax1.set_aspect(1/20)
    ax2 = plt.subplot(122, projection=ccrs.PlateCarree())
    display_Horus.plot_cross_section(
            "ZDR",
            [start_lat,start_lon],
            [end_lat,end_lon],
            ax = ax2,
            x_axis="lat",
            cmap="pyart_NWSRef",
            vmin=-1,
            vmax=7,
            title = 'Horus ZDR \n'+ time + 'Z',
            colorbar_flag = False
        )
    ax2.set_ylim(0,10)
    ax2.set_aspect(1/20)
    
    cax = ax1.inset_axes([cbar_ax[1]+0.4, -2, 1.1, 0.6],transform=ax1.transData)
    m = cm.ScalarMappable(cmap='pyart_NWSRef')
    m.set_clim(-1,7)
    fig.colorbar(m,cax=cax,orientation='horizontal')
    if save_fig:
        plt.savefig('ZDR_Comparison_RHI_'+str(az)+'_'+date+'_'+time)
    else:
        plt.show()
    
def SidebySidePlot(radar,dateHrus,Hrus_time,ele,xlims,ylims,figpath,save_fig):
    os.chdir(figpath)
    fig = plt.figure(figsize=(30,12))
    ax1 = plt.subplot(121)
    display = pyart.graph.RadarMapDisplay(radar)
    display.plot_ppi('Reflectivity',ax=ax1,vmin=-10,vmax=70,cmap='pyart_NWSRef')
    plt.xlim(xlims)
    plt.ylim(ylims)
    
    ax2 = plt.subplot(122)
    display = pyart.graph.RadarMapDisplay(radar)
    display.plot_ppi('HCA',ax=ax2,vmin=0.5,vmax=14.5,cmap=hca_cmap,ticklabs=labs,ticks=ticks)
    plt.xlim(xlims)
    plt.ylim(ylims)
    if save_fig is True:
        plt.savefig(dateHrus+'_'+str(Hrus_time)+'_'+ele[0:2]+ele[3:]+'_Reflectivity',bbox_inches = 'tight')
    else:
        plt.show(fig)
    plt.close(fig)
    
    fig = plt.figure(figsize=(30,12))
    ax1 = plt.subplot(121)
    display = pyart.graph.RadarMapDisplay(radar)
    display.plot_ppi('Reflectivity',ax=ax1,vmin=-10,vmax=70,cmap='pyart_NWSRef')
    plt.xlim(xlims)
    plt.ylim(ylims)
    
    ax2 = plt.subplot(122)
    display = pyart.graph.RadarMapDisplay(radar)
    display.plot_ppi('ZDR',ax=ax2,vmin=-1,vmax=7,cmap='pyart_NWSRef')
    plt.xlim(xlims)
    plt.ylim(ylims)
    if save_fig is True:
        plt.savefig(dateHrus+'_'+str(Hrus_time)+'_'+ele[0:2]+ele[3:]+'_ReflectivityZDR',bbox_inches = 'tight')
    else:
        plt.show(fig)
    plt.close(fig)
    
    fig = plt.figure(figsize=(30,12))
    ax1 = plt.subplot(121)
    display = pyart.graph.RadarMapDisplay(radar)
    display.plot_ppi('Reflectivity',ax=ax1,vmin=-10,vmax=70,cmap='pyart_NWSRef')
    plt.xlim(xlims)
    plt.ylim(ylims)
    
    ax2 = plt.subplot(122)
    display = pyart.graph.RadarMapDisplay(radar)
    display.plot_ppi('corrected_velocity',ax=ax2,vmin=-30,vmax=30,cmap='pyart_BuDRd18')
    plt.xlim(xlims)
    plt.ylim(ylims)
    if save_fig is True:
        plt.savefig(dateHrus+'_'+str(Hrus_time)+'_'+ele[0:2]+ele[3:]+'_ReflectivityVelocity',bbox_inches = 'tight')
    else:
        plt.show(fig)
    plt.close(fig)
    
    #plt.show(fig)
    
def fan_plots(radar,radar_single,dateHrus,Hrus_time,xlims,ylims,xlim_r,zlim,figpath,save_fig):
    os.chdir(figpath)
    azs = radar_single.azimuth['data']
    
    
    for i in range(0,len(azs)):
        xsect = pyart.util.cross_section_ppi(radar,[azs[i]])
        points = xsect.get_gate_x_y_z(0)
        
        display_ppi = pyart.graph.RadarDisplay(radar)
        fig = plt.figure(figsize=(30,12))
        ax1 = plt.subplot(121)
        display_ppi.plot_ppi('Reflectivity',3,ax=ax1,vmin=-10,vmax=70,cmap='pyart_NWSRef',colorbar_flag=False)
        plt.scatter(points[0][1,:]/1000,points[1][1,:]/1000,s=0.5,c='k')
        plt.xlim(xlims)
        plt.ylim(ylims)
        ax2 = plt.subplot(122)
        display = pyart.graph.RadarDisplay(xsect)
        display.plot('Reflectivity',ax=ax2,vmin=-10,vmax=70,cmap='pyart_NWSRef')
        plt.ylim(zlim)
        plt.xlim(xlim_r)
        if save_fig is True:
            plt.savefig(dateHrus+'_'+str(Hrus_time)+'_'+str(azs[i])[0:3]+'_ReflectivityFan',bbox_inches = 'tight')
        else:
            plt.show(fig)
        plt.close(fig)
        
        display_ppi = pyart.graph.RadarDisplay(radar)
        fig = plt.figure(figsize=(30,12))
        ax1 = plt.subplot(121)
        display_ppi.plot_ppi('ZDR',3,ax=ax1,vmin=-1,vmax=7,cmap='pyart_NWSRef',colorbar_flag=False)
        plt.scatter(points[0][1,:]/1000,points[1][1,:]/1000,s=0.5,c='k')
        plt.xlim(xlims)
        plt.ylim(ylims)
        ax2 = plt.subplot(122)
        display = pyart.graph.RadarDisplay(xsect)
        display.plot('ZDR',ax=ax2,vmin=-1,vmax=7,cmap='pyart_NWSRef')
        plt.ylim(zlim)
        plt.xlim(xlim_r)
        if save_fig is True:
            plt.savefig(dateHrus+'_'+str(Hrus_time)+'_'+str(azs[i])[0:3]+'_ZDRFan',bbox_inches = 'tight')
        else:
            plt.show(fig)
        plt.close(fig)
        
        display_ppi = pyart.graph.RadarDisplay(radar)
        fig = plt.figure(figsize=(30,12))
        ax1 = plt.subplot(121)
        display_ppi.plot_ppi('corrected_velocity',3,ax=ax1,vmin=-60,vmax=60,cmap='pyart_BuDRd18',colorbar_flag=False)
        plt.scatter(points[0][1,:]/1000,points[1][1,:]/1000,s=0.5,c='k')
        plt.xlim(xlims)
        plt.ylim(ylims)
        ax2 = plt.subplot(122)
        display = pyart.graph.RadarDisplay(xsect)
        display.plot('corrected_velocity',ax=ax2,vmin=-30,vmax=30,cmap='pyart_BuDRd18')
        plt.ylim(zlim)
        plt.xlim(xlim_r)
        if save_fig is True:
            plt.savefig(dateHrus+'_'+str(Hrus_time)+'_'+str(azs[i])[0:3]+'_VelocityFan',bbox_inches = 'tight')
        else:
            plt.show(fig)
        plt.close(fig)

def building_blockage(Horusrad,sweep_Hrus,lon,lat,build_lon,build_lat,dateHrus,Hrus_time,save_fig):
    ### Reflectivity
    fig = plt.figure()
    projection = ccrs.LambertConformal(
            central_latitude=Horusrad.latitude["data"][0],
            central_longitude=Horusrad.longitude["data"][0],
            )
    display2 = pyart.graph.RadarMapDisplay(Horusrad)
    ax2 = fig.add_subplot(111,projection=ccrs.PlateCarree())
    display2.plot_ppi_map('Reflectivity',ax=ax2,vmin=-10,vmax=70,cmap='pyart_NWSRef',
                         projection=projection,
                         min_lon=lon[0],max_lon=lon[1],
                         min_lat=lat[0],max_lat=lat[1]
                         )
    ax2.scatter(build_lon,build_lat,marker='*',s=100,color='black',transform=ccrs.Geodetic())
    if save_fig:
        plt.savefig(dateHrus+'_'+str(Hrus_time)+'_Reflectivity_Comparison',bbox_inches = 'tight')
    else:
        plt.show()
    plt.close()
    
#############################
####       CAPPIS       ####
############################

def read_CAPPI(path,alt,idz,dateHorus,rad_name):
    os.chdir(path)
    os.chdir('ReflectivityCAPPI')
    os.chdir(alt)
    l2=np.sort([f for f in os.listdir('.') if re.search('.netcdf', f)])
    file_path=l2[idz]
    #print(file_path)
    
    with netCDF4.Dataset(file_path,mode='r') as data:
            #print(data)
            Ref=np.array(data['ReflectivityCAPPI'])

            I,J=Ref.shape
            #print(I,J)
            #Az=np.array(nc['Azimuth'])

            # Fill the Radar Object with important fields
            rad = pyart.testing.make_empty_ppi_radar(J, I, 1)
            rad.ngates = J
            rad.nrays = I
            time_dict=get_metadata('Time')
            #if epoch == 1: 
            #    time_1 = datetime.fromtimestamp(int(np.array([data.Time])))
            #    time_1 = time_1 - timedelta(hours=5, minutes=0)
            #    time_dict['data'] = np.array([int(time_1.timestamp())])
            #else:
            time_dict['data']=np.array(np.repeat(data.Time,I))
            #print(data.Time)
            #sys.exit()
            time_dict['units']='seconds since 1970-01-01T00:00:00'
            time_dict['calendar']='gregorian'
            rad.time = time_dict
            if dateHorus == '20240428':
                rad.latitude['data'] = np.array([35.1864])
                rad.longitude['data'] = np.array([-97.4459])
            else:
                rad.latitude['data'] = np.array([data.Latitude])
                rad.longitude['data'] = np.array([data.Longitude])
            rad.sweep_number['data'] = np.array([1])
            rad.azimuth['data'] = np.array(data['Azimuth']) 
            #print(np.array(data['Azimuth']))
            #rad.azimuth['data'] = np.array(60.25)
            #azs = rad.azimuth['data'][0]
            #rad.azimuth['data'][0:-1] = rad.azimuth['data'][1:]
            #rad.azimuth['data'][-1]=azs
            rad.fixed_angle['data']=np.array([alt]).astype(float)
            #print(np.ones(I)*np.array([elevations]).astype(float))
            rad.elevation['data'] = np.ones(I) * np.array([alt]).astype(float)
            #rad.sweep_mode['data'] = np.array(data['sweep_mode'])
            rad.instrument_name=rad_name
            rad.metadata['instrument_name']=rad_name
            parameters = {}
            frequency = {'data': int(2705000000), 'units': 'Hertz'}
            parameters['frequency'] = frequency
            rad.instrument_parameters=parameters
            
            numgates=(len(data.dimensions['Gate']))
            gatewidth=np.array(data['GateWidth'])[0]
            start=data.RangeToFirstGate
            ranges=np.ones(numgates)+start
            for i in range(1,len(ranges)):
                ranges[i]=ranges[i-1]+gatewidth
            
            range_dict=get_metadata('Range')
            range_dict['data']=ranges
            range_dict['Units'] = 'Meters'
            rad.range=range_dict
            
            rad.init_gate_altitude()
            rad.init_gate_longitude_latitude()
            
            # Now get other data - ie Radar Fields
            # Reflectivity
            ref_dict = get_metadata('Reflectivity')
            ref_dict['data'] = np.array(data['ReflectivityCAPPI'])
            ref_dict['data'][ref_dict['data']<-10000] = np.nan
            ref_dict['units']='dBZ'
            rad.add_field('Reflectivity', ref_dict)
            
    #RhoHV
    os.chdir(path)
    os.chdir('RhoHVCAPPI')
    os.chdir(alt)
    l2=np.sort([f for f in os.listdir('.') if re.search('.netcdf', f)])
    file_path=l2[idz]
    with netCDF4.Dataset(file_path,mode='r') as data:
        rho_dict = get_metadata('RhoHV')
        try:
            rho_dict['data'] = np.array(data['RhoHVCAPPI'])
            print(np.nanmin(data['RhoHVCAPPI']))
            rho_dict['data'][rho_dict['data']==-99903] = np.nan
        except:
            rho_dict['data'] = np.array(data['RhoHV_EXTRA'])
        rho_dict['units'] = 'unitless'
        rad.add_field('RhoHV',rho_dict)
    
    #ZDR
    os.chdir(path)
    os.chdir('ZdrCAPPI')
    os.chdir(alt)
    l2=np.sort([f for f in os.listdir('.') if re.search('.netcdf', f)])
    file_path=l2[idz]
    with netCDF4.Dataset(file_path,mode='r') as data:
        zdr_dict = get_metadata('Zdr')
        try:
            zdr_dict['data'] = np.array(data['ZdrCAPPI'])
            zdr_dict['data'][zdr_dict['data']==-99903] = np.nan
        except:
            zdr_dict['data'] = np.array(data['Zdr_EXTRA'])  
        zdr_dict['units'] = 'dB'
        rad.add_field('ZDR',zdr_dict)
        
    # Velocity
    os.chdir(path)
    os.chdir('AliasedVelocityCAPPI')
    os.chdir(alt)
    l2=np.sort([f for f in os.listdir('.') if re.search('.netcdf', f)])
    file_path=l2[idz]
    with netCDF4.Dataset(file_path,mode='r') as data:
        zdr_dict = get_metadata('velocity')
        try:
            zdr_dict['data'] = np.array(data['AliasedVelocityCAPPI'])
            zdr_dict['data'][zdr_dict['data']==-99903] = np.nan
        except:
            zdr_dict['data'] = np.array(data['AliasedVelocity_EXTRA'])
        zdr_dict['units'] = 'm/s'
        rad.add_field('velocity',zdr_dict)
    
    # PhiDP
    os.chdir(path)
    os.chdir('PhiDPCAPPI')
    os.chdir(alt)
    l2=np.sort([f for f in os.listdir('.') if re.search('.netcdf', f)])
    file_path=l2[idz]
    with netCDF4.Dataset(file_path,mode='r') as data:
        phi_dict = get_metadata('PhiDP')
        phi_dict['data'] = np.array(data['PhiDPCAPPI'])
        phi_dict['data'][phi_dict['data']==-99903] = np.nan
        phi_dict['units'] = 'degrees'
        rad.add_field('PhiDP',phi_dict)
        
    os.chdir(path)
    os.chdir('DHCACAPPI')
    os.chdir(alt)
    l2=np.sort([f for f in os.listdir('.') if re.search('.netcdf', f)])
    file_path=l2[idz]
    #print(file_path)
    with netCDF4.Dataset(file_path,mode='r') as data:
        zdr_dict = get_metadata('HCA')
        zdr_dict['data'] = np.array(data['DHCACAPPI'])
        zdr_dict['data'][zdr_dict['data']==-99903] = np.nan
        zdr_dict['units'] = 'unitless'
        rad.add_field('HCA',zdr_dict)
            
    return rad
   
def plot_CAPPI(radar,rad_MESH,rad_Echo,xlims,ylims,figpath,dateHrus,Hrus_time,ele,echo_thresh,save_fig,plot_MESH):
    cdict = {'red':   [[0.0,  0.0, 0.0],
                   [0.5,  0.0, 0.0],
                   [0.85,  1.0, 1.0],
                   [1.0,  1.0, 1.0]],
         'green': [[0.0,  0.0, 0.0],
                   [0.25, 0.0, 0.0],
                   [0.5,  0.0, 0.0],
                   [0.85,  1.0, 1.0],
                   [1.0,  0.0, 0.0]],
         'blue':  [[0.0,  1.0, 1.0],
                   [0.5,  1.0, 1.0],
                   [1.0,  0.1, 0.1]]}
    cc_cmap = colors.LinearSegmentedColormap('testCmap', segmentdata=cdict, N=256)
    ticks = np.arange(1,15)
    labs = ['L/MR','HR','R/Ha','BD','AP','BI','UK','NE','DS','WS','CY','','GR','+']
    hca_cmap = colors.ListedColormap(['#AAAAAA','#0000FE','#0057FF','#00A9FF','#00FFFF','#56FFAA','#ACFD56','#FFFF00','#FEAB02','#FF5700','#FF00FF','#B300B3','#DC143C','#FF00FF'])
    hca_cmap = colors.ListedColormap(['limegreen','green','#DC143C','#FFFF00','dimgray','lightgray','#B300B3','white','#00FFFF','#0000FE','mediumturquoise','#000000','#FF00FF','#000000'])

    os.chdir(figpath)
    fig = plt.figure()
    display = pyart.graph.RadarMapDisplay(radar)
    display.plot_ppi_map('Reflectivity',vmin=-10,vmax=70,cmap='pyart_NWSRef',
                         min_lat = ylims[0],max_lat=ylims[1],min_lon=xlims[0],max_lon=xlims[1])

    if save_fig is True:
        plt.savefig(dateHrus+'_'+str(Hrus_time)+'_'+ele[0:2]+ele[3:]+'_ReflectivityCAPPI',bbox_inches = 'tight')
    else:
        plt.show(fig)
    plt.close(fig)
    #plt.show(fig)
    
    fig = plt.figure()
    display = pyart.graph.RadarMapDisplay(radar)
    display.plot_ppi_map('velocity',vmin=-30, vmax=30,cmap = 'pyart_BuDRd18',
                         min_lat = ylims[0],max_lat=ylims[1],min_lon=xlims[0],max_lon=xlims[1])
    if save_fig is True:
        plt.savefig(dateHrus+'_'+str(Hrus_time)+'_'+ele[0:2]+ele[3:]+'_VelocityCAPPI',bbox_inches = 'tight')
    else:
        plt.show(fig)
    plt.close(fig)
    #plt.show(fig)
    
    fig = plt.figure()
    display.plot_ppi_map('ZDR',vmin=-1,vmax=7,cmap='pyart_NWSRef',
                         min_lat = ylims[0],max_lat=ylims[1],min_lon=xlims[0],max_lon=xlims[1])
    if save_fig is True:
        plt.savefig(dateHrus+'_'+str(Hrus_time)+'_'+ele[0:2]+ele[3:]+'_ZDRCAPPI',bbox_inches = 'tight')
    else:
        plt.show(fig)
    plt.close(fig)
    #plt.show(fig)
    
    fig = plt.figure()
    display.plot_ppi_map('RhoHV',vmin=0.5,vmax=1,cmap = cc_cmap,
                         min_lat = ylims[0],max_lat=ylims[1],min_lon=xlims[0],max_lon=xlims[1])
    if save_fig is True:
        plt.savefig(dateHrus+'_'+str(Hrus_time)+'_'+ele[0:2]+ele[3:]+'_RhoHVCAPPI',bbox_inches = 'tight')
    else:
        plt.show(fig)
    plt.close(fig)
    #plt.show(fig)
    
    fig = plt.figure()
    display.plot_ppi_map('PhiDP',vmin=-180,vmax=180,cmap = 'pyart_RRate11',
                         min_lat = ylims[0],max_lat=ylims[1],min_lon=xlims[0],max_lon=xlims[1])
    if save_fig is True:
        plt.savefig(dateHrus+'_'+str(Hrus_time)+'_'+ele[0:2]+ele[3:]+'_PhiDPCAPPI',bbox_inches = 'tight')
    else:
        plt.show(fig)
    plt.close(fig)
    '''
    fig = plt.figure()
    display.plot_ppi('KDP',vmin=-1,vmax=10,cmap = 'pyart_NWSRef')
    plt.xlim(xlims)
    plt.ylim(ylims)
    if save_fig is True:
        plt.savefig(dateHrus+'_'+str(Hrus_time)+'_'+ele[0:2]+ele[3:]+'_KDP',bbox_inches = 'tight')
    else:
        plt.show(fig)
    plt.close(fig)'''
    
    if plot_MESH is True:
        fig = plt.figure()
        #ax = fig.add_subplot(111)
        display2 = pyart.graph.RadarMapDisplay(rad_MESH)
        display.plot_ppi('Reflectivity',vmin=-10,vmax=70,cmap='pyart_NWSRef',alpha=0.1,colorbar_flag=False)
        display2.plot_ppi('MESH',cmap='pyart_NWSRef',vmin=0,vmax=100)
        plt.xlim(xlims)
        plt.ylim(ylims)
        if save_fig is True:
            plt.savefig(dateHrus+'_'+str(Hrus_time)+'_'+ele[0:2]+ele[3:]+'_MESHRef',bbox_inches = 'tight')
        else:
            plt.show(fig)
        plt.close(fig)
        #plt.show(fig)
        
        fig = plt.figure()
        #ax = fig.add_subplot(111)
        display2 = pyart.graph.RadarMapDisplay(rad_MESH)
        #display.plot_ppi('Reflectivity',vmin=-10,vmax=70,cmap='pyart_NWSRef',alpha=0.1,colorbar_flag=False)
        display2.plot_ppi('MESH',cmap='pyart_NWSRef',vmin=0,vmax=100)
        plt.xlim(xlims)
        plt.ylim(ylims)
        if save_fig is True:
            plt.savefig(dateHrus+'_'+str(Hrus_time)+'_MESH',bbox_inches = 'tight')
        else:
            plt.show(fig)
        plt.close(fig)
    
    for j in range(0,len(echo_thresh)):
        fig = plt.figure()
        display3 = pyart.graph.RadarMapDisplay(rad_Echo)
        display3.plot_ppi('EchoTop_'+echo_thresh[j],cmap='pyart_NWSRef',vmin=0,vmax=20)
        plt.xlim(xlims)
        plt.ylim(ylims)
        if save_fig is True:
            plt.savefig(dateHrus+'_'+str(Hrus_time)+'_EchoTop_'+echo_thresh[j],bbox_inches = 'tight')
        else:
            plt.show(fig)
        plt.close(fig)
    # Sets up the HCA colormaps to be used later
    ticks = np.arange(1,15)
    labs = ['L/MR','HR','R/Ha','BD','AP','BI','UK','NE','DS','WS','CY','','GR','+']
    hca_cmap = colors.ListedColormap(['#AAAAAA','#0000FE','#0057FF','#00A9FF','#00FFFF','#56FFAA','#ACFD56','#FFFF00','#FEAB02','#FF5700','#FF00FF','#B300B3','#DC143C','#FF00FF'])
    hca_cmap = colors.ListedColormap(['limegreen','green','#DC143C','#FFFF00','dimgray','lightgray','#B300B3','white','#00FFFF','#0000FE','mediumturquoise','#000000','#FF00FF','#000000'])
    
    fig = plt.figure()
    #ax = fig.add_subplot(111)
    #display2 = pyart.graph.RadarMapDisplay(rad_MESH)
    display.plot_ppi_map('HCA',vmin=0.5,vmax=14.5,cmap=hca_cmap,ticklabs=labs,ticks=ticks,
                         min_lat = ylims[0],max_lat=ylims[1],min_lon=xlims[0],max_lon=xlims[1])
    #display2.plot_ppi('MESH',cmap='pyart_NWSRef',vmin=0,vmax=10)
    if save_fig is True:
        plt.savefig(dateHrus+'_'+str(Hrus_time)+'_'+ele[0:2]+ele[3:]+'_DHCA',bbox_inches = 'tight')
    else:
        plt.show(fig)
    plt.close(fig)
     
    
def Plot_CAPPI_Comparison(radar,radar_ktlx,lon,lat,figpath,dateHrus,Hrus_time,cbar_ax,save_fig):
    os.chdir(figpath)
    display = pyart.graph.RadarMapDisplay(radar)
    display2 = pyart.graph.RadarMapDisplay(radar_ktlx)
    fig = plt.figure(figsize=(16,12))
    ax1 = plt.subplot(121,projection=ccrs.PlateCarree())
    display2.plot_ppi_map('Reflectivity',sweep=0,ax=ax1,vmin=-10,vmax=70,cmap='pyart_NWSRef',
                         min_lon = lon[0],max_lon=lon[1],
                         min_lat = lat[0], max_lat=lat[1],
                         colorbar_flag=False)
    ax2 = plt.subplot(122,projection=ccrs.PlateCarree())
    display.plot_ppi_map('Reflectivity',sweep=0,ax=ax2,vmin=-10,vmax=70,cmap='pyart_NWSRef',
                         min_lon=lon[0],max_lon=lon[1],
                         min_lat=lat[0],max_lat=lat[1],
                         colorbar_flag=False)
    
    cax = ax1.inset_axes([cbar_ax[0], cbar_ax[1], 0.89, 0.03],transform=ax1.transData)
    m = cm.ScalarMappable(cmap='pyart_NWSRef')
    m.set_clim(-10,70)
    fig.colorbar(m,cax=cax,orientation='horizontal')
    if save_fig:
        plt.savefig(dateHrus+'_'+str(Hrus_time)+'_ReflectivityCAPPI_Comparison',bbox_inches = 'tight')
    else:
        plt.show()
    plt.close()
    