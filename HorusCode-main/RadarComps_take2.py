#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 27 13:17:00 2023
2 August: Expanded so that you can work with multiple cases (i.e. Feb cases)
@author: lshedd123
"""

import matplotlib.pyplot as plt
import os,re
import numpy as np
import netCDF4,gzip
import pyart
import math
from pyart.config import get_metadata
from matplotlib import colors
import time
from matplotlib.widgets import LassoSelector
from matplotlib import path
from datetime import datetime,timezone,timedelta

from six.moves import input
from matplotlib.path import Path
import netCDF4 as nc
import sys
import requests
import geopy.distance

import warnings
warnings.filterwarnings('ignore')

######## USER INPUT VARIABLES ############
#These are the things that may need to change depending on which case you are working with
# set epoch as 1 for the moment if using the 11 May dataset
epoch = 0
# which nexrad you are using
radar_name = 'KTLX'
# the date you want to use (these should match) 
date = '20230808_Processed/output'
dateKTLX = '20230808'

x_lims = [-98.5,-97.5]
y_lims = [35,35.4]

cross_ylim = [0,5]
cross_xlim = [0,30]
#aspect = (cross_ylim[1] - cross_ylim[0]) / (cross_xlim[1] - cross_xlim[0]) + 0.3
aspect = 1

#Limits for Gridding in z,y,x (min/max location in meters)
# Currently set to be 181x181 x,y grids but can be changed 
grid_lims = ((500, 14000), (-70000, 30000), (-80000, 30000)) ## good for 11 May
#grid_lims = ((500, 14000), (-30000, 70000), (-10000, 90000)) ## good for 13 Dec 

figfolder = '/Users/lshedd123/Dropbox/Horus Share Folder/Horus Comparisons/' + date

######## END INPUT VARIABLES #############

###### One Important Note: #############
# KTLX files generally have some extra
# on either end of the Horus times so be 
# mindful of this
#######################################
# Set Up Necessary Definitions
def read_KOUN(idz):
    #Reffolder = '/Volumes/Laura/Data/KOUN_20130519/Reflectivity/'+eles
    os.chdir(Reffolder_KTLX)
    l = np.sort([f for f in os.listdir('.') if re.search('.netcdf', f)])
    file=l[idz]
    filename=os.path.basename(file)
    time_str=filename[9:15]
    date_str=filename[0:8]
    time_cat=date_str+time_str
    time_initial = datetime.strptime(time_cat, "%Y%m%d%H%M%S")
    print('Initial Time: '+str(time_initial))
    indices_lvls = []
    lvls=[]
    indices_lvls = np.append(indices_lvls,0)
    lvls = np.append(lvls,levels_KTLX[0])
        
    return time_initial, time_str, indices_lvls, lvls

def read_KOUN_lvls(lvl,idz,start2,end):
     #print(lvl) 
     folder = '/Volumes/MiniHorus/Data/'+radar_name+'_'+dateKTLX+'/Reflectivity/'+str(lvl)
     os.chdir(folder)
     l = np.sort([f for f in os.listdir('.') if re.search('.netcdf', f)])
     file=l[idz]
     #print('Def: '+file)
     #print(os.getcwd())
     #print(file)
 
     #with gzip.open(file) as gz:
     with netCDF4.Dataset(file,mode='r') as data:
            # print(data)
             Ref=np.array(data['Reflectivity'])
             #Ref[Ref==-99900]=np.nan
             I,J=Ref.shape
     
     # Fill the Radar Object with important fields
             rad = pyart.testing.make_empty_ppi_radar(J, I, 1)
             
             #print(J,I)
             time_dict=get_metadata('Time')
             time_dict['data']=np.array([data.Time]*I)
             time_dict['units']='seconds since 1970-01-01T00:00:00Z'
             time_dict['calendar']='gregorian'
             rad.time = time_dict
             rad.latitude['data'] = np.array([data.Latitude])
             rad.longitude['data'] = np.array([data.Longitude])
             rad.fixed_angle['data'] = np.array([data.Elevation])
             rad.sweep_number['data'] = np.array([k])
             rad.sweep_start_ray_index['data'] = np.array([start2])
             rad.sweep_end_ray_index['data'] = np.array([end])
             rad.altitude['data']=np.array([data.Height])
             rad.azimuth['data'] = np.array(data['Azimuth'])
             rad.elevation['data']=np.ones(len(data['Azimuth']))*np.array([data.Elevation])
             rad.elevation['units']='degrees'
             azs = rad.azimuth['data'][0]
             rad.azimuth['data'] = np.array(data['Azimuth'])
             rad.instrument_name='KTLX'
             rad.metadata['instrument_name']='KTLX'
             parameters = {}
             frequency = {'data': int(2705000000), 'units': 'Hertz'}
             parameters['frequency'] = frequency
             rad.instrument_parameters=parameters
             rad.scan_type='sector'
             
             # Lets Have Fun With Ranges
             numgates=(len(data.dimensions['Gate']))
             gatewidth=np.array(data['GateWidth'])[0]
             start=data.RangeToFirstGate
             ranges=np.ones(numgates)+start
             for i in range(1,len(ranges)):
                 ranges[i]=ranges[i-1]+gatewidth
             
             
             range_dict=get_metadata('Range')
             range_dict['data']=ranges
             range_dict['units']='meters'
             rad.range=range_dict
             
             rad.init_gate_altitude()
             rad.init_gate_longitude_latitude()
             rad.init_gate_x_y_z()
             rad.get_elevation(0)
             
             # Now get other data - ie Radar Fields
             ref_dict = get_metadata('Reflectivity')
             ref_dict['data'] = Ref
             #print(Ref.shape)
             ref_dict['units']='dBZ'
             rad.fields = {'Reflectivity': ref_dict}
             
     folder='/Volumes/MiniHorus/Data/'+radar_name+'_'+dateKTLX+'/DHCA/'+str(lvl)
     os.chdir(folder)
     l4=np.sort([f for f in os.listdir('.') if re.search('.netcdf', f)])
     #print(l4)
     file=l4[idz]
     #print(file)
     
     with netCDF4.Dataset(file,mode='r') as data:
         vel_dict=get_metadata('DHCA')
         vel_dict['data']=np.array(data['DHCA'])
         #print(vel_dict['data'].shape)
         rad.add_field('DHCA',vel_dict)
     '''
     folder='/Volumes/MiniHorus/Data/'+radar_name+'_'+dateKTLX+'/RhoHV/'+str(lvl)      
     os.chdir(folder)
     l2=np.sort([f for f in os.listdir('.') if re.search('.netcdf', f)])
     file=l2[idz]
     print(file)
     
     #with gzip.open(file) as gz:
     with netCDF4.Dataset(file,mode='r') as data:
            cc_dict=get_metadata('Correlation_Coefficient')
            cc_dict['data']=np.array(data['RhoHV'])
            cc_dict['units']='unitless'
            try:
                #print(np.array(data['RhoHV']).shape)
                rad.add_field('RhoHV',cc_dict)
            except:
                #print('CC Failed, Now Correcting')
                # Note, we defined I and J for Ref
                A,B = np.shape(np.array(data['RhoHV']))
                CC_nocorr = np.array(data['RhoHV'])
                x_miss = I-A
                y_miss = J-B
                #print(x_miss,y_miss)
                if y_miss != 0:
                    adds = np.ones([I,y_miss])*np.nan
                    CC_corr = np.hstack([CC_nocorr,adds])
                    CC_nocorr = CC_corr
                    #print('Whoops')
                #for m in range((0,x_miss)):
                if x_miss != 0:
                    add = np.ones([x_miss,J])*np.nan
                    CC_corr = np.vstack([CC_nocorr,add])
                cc_dict['data'] = CC_corr
                rad.add_field('RhoHV',cc_dict)
     
     folder='/Volumes/MiniHorus/Data/'+radar_name+'_'+date+'/Zdr/'+str(lvl)
     os.chdir(folder)
     l3=np.sort([f for f in os.listdir('.') if re.search('.netcdf', f)])
     file=l3[idz]
     
     #with gzip.open(file) as gz:
     with netCDF4.Dataset(file,mode='r') as data:
            zdr_dict=get_metadata('Differential_Reflectivity')
            zdr_dict['data']=np.array(data['Zdr'])
            zdr_dict['units']='dB'
            try:
                rad.add_field('ZDR',zdr_dict)
            except:
                #print('ZDR Failed, Now Correcting')
                # Note, we defined I and J for Ref
                A,B = np.shape(np.array(data['Zdr']))
                Zdr_nocorr = np.array(data['Zdr'])
                x_miss = I-A
                y_miss = J-B
                if y_miss != 0:
                    adds = np.ones([I,y_miss])*np.nan
                    Zdr_corr = np.hstack([Zdr_nocorr,adds])
                    Zdr_nocorr = Zdr_corr
                #for m in range((0,x_miss)):
                    if x_miss !=0:
                        add = np.ones([x_miss,J])*np.nan
                        Zdr_corr = np.vstack([Zdr_nocorr,add])
                zdr_dict['data'] = Zdr_corr
                rad.add_field('ZDR',zdr_dict)
    
    # Now Get the Actual Elevations
     #xx,yy,zz = rad.get_gate_x_y_z(0)
     #print(zz[0,:)
     #rad.gate_z = zz
     #rad.elevation['data'] = zz[0,:]
     #print(len(zz[:,0]))
     #rad.elevation['units']='meters'
             
     folder='/Volumes/MiniHorus/Data/'+radar_name+'_'+date+'/AliasedVelocity/'+str(lvl)
     os.chdir(folder)
     l3=np.sort([f for f in os.listdir('.') if re.search('.netcdf', f)])
     file=l3[idz]
     
     #with gzip.open(file) as gz:
     with netCDF4.Dataset(file,mode='r') as data:
            vel_dict=get_metadata('Aliased Velocity')
            vel_dict['data']=np.array(data['AliasedVelocity'])
            try:
                rad.add_field('Aliased Velocity',vel_dict)
            except:
                #print('Velocity Failed, Now Correcting')
                # Note, we defined I and J for Ref
                A,B = np.shape(np.array(data['AliasedVelocity']))
                vel_nocorr = np.array(data['AliasedVelocity'])
                x_miss = I-A
                y_miss = J-B
                if y_miss != 0:
                    #print('Whoops')
                    adds = np.ones([I,y_miss])*np.nan
                    vel_corr = np.hstack([vel_nocorr,adds])
                    vel_nocorr = vel_corr
                #for m in range((0,x_miss)):
                    if x_miss != 0:
                        add = np.ones([x_miss,J])*np.nan
                        vel_corr = np.vstack([vel_nocorr,add])
                vel_dict['data'] = vel_corr
                rad.add_field('Aliased Velocity',vel_dict)
     
    folder='/Volumes/Laura/Data/KOUN_20170326/Velocity/'+str(lvl)
     os.chdir(Velfolder)
     l4=np.sort([f for f in os.listdir('.') if re.search('.netcdf', f)])
     #print(l4)
     file=l4[idz]
     
     with netCDF4.Dataset(file,mode='r') as data:
         vel_dict=get_metadata('Velocity')
         vel_dict['data']=np.array(data['AliasedVelocity'])
         rad.add_field('Velocity',vel_dict)'''
         
     return rad
 
def read_KOUN_lvls_dp(lvl,idz,start2,end):
     #print(lvl) 
     folder = '/Volumes/MiniHorus/Data/'+radar_name+'_'+dateKTLX+'/RhoHV/'+str(lvl)
     os.chdir(folder)
     l = np.sort([f for f in os.listdir('.') if re.search('.netcdf', f)])
     file=l[idz]
     #print('Def: '+file)
     #print(os.getcwd())
    # print(file)
     
     #with gzip.open(file) as gz:
     with netCDF4.Dataset(file,mode='r') as data:
            # print(data)
             Ref=np.array(data['RhoHV'])
             Ref[Ref==-99900]=np.nan
             I,J=Ref.shape
             #print(Ref.shape)
     
     # Fill the Radar Object with important fields
             rad = pyart.testing.make_empty_ppi_radar(J, I, 1)
             time_dict=get_metadata('Time')
             time_dict['data']=np.array([data.Time]*I)
             time_dict['units']='seconds since 1970-01-01T00:00:00Z'
             time_dict['calendar']='gregorian'
             rad.time = time_dict
             rad.latitude['data'] = np.array([data.Latitude])
             rad.longitude['data'] = np.array([data.Longitude])
             rad.fixed_angle['data'] = np.array([data.Elevation])
             rad.sweep_number['data'] = np.array([k])
             rad.sweep_start_ray_index['data'] = np.array([start2])
             rad.sweep_end_ray_index['data'] = np.array([end])
             rad.altitude['data']=np.array([data.Height])
             rad.azimuth['data'] = np.array(data['Azimuth'])
             rad.elevation['data']=np.ones(len(data['Azimuth']))*np.array([data.Elevation])
             rad.elevation['units']='degrees'
             azs = rad.azimuth['data'][0]
             rad.azimuth['data'] = np.array(data['Azimuth'])
             rad.instrument_name='KTLX'
             rad.metadata['instrument_name']='KTLX'
             parameters = {}
             frequency = {'data': int(2705000000), 'units': 'Hertz'}
             parameters['frequency'] = frequency
             rad.instrument_parameters=parameters
             rad.scan_type='sector'
             
             # Lets Have Fun With Ranges
             numgates=(len(data.dimensions['Gate']))
             gatewidth=np.array(data['GateWidth'])[0]
             start=data.RangeToFirstGate
             ranges=np.ones(numgates)+start
             for i in range(1,len(ranges)):
                 ranges[i]=ranges[i-1]+gatewidth
             
             
             range_dict=get_metadata('Range')
             range_dict['data']=ranges
             range_dict['units']='meters'
             rad.range=range_dict
             
             rad.init_gate_altitude()
             rad.init_gate_longitude_latitude()
             rad.init_gate_x_y_z()
             rad.get_elevation(0)
             
             # Now get other data - ie Radar Fields
             ref_dict = get_metadata('RhoHV')
             ref_dict['data'] = Ref
             #print(Ref.shape)
             ref_dict['units']='unitless'
             rad.fields = {'RhoHV': ref_dict}
         
     folder='/Volumes/MiniHorus/Data/'+radar_name+'_'+dateKTLX+'/Zdr/'+str(lvl)
     os.chdir(folder)
     l3=np.sort([f for f in os.listdir('.') if re.search('.netcdf', f)])
     file=l3[idz]
     #print(str(lvl))
     #print(file)
     
     #with gzip.open(file) as gz:
     with netCDF4.Dataset(file,mode='r') as data:
            zdr_dict=get_metadata('Differential_Reflectivity')
            zdr_dict['data']=np.array(data['Zdr'])
            zdr_dict['data'][zdr_dict['data']<-99900]=np.nan
            zdr_dict['units']='dB'
            try:
                rad.add_field('ZDR',zdr_dict)
            except:
                print('ZDR Failed, Now Correcting')
                # Note, we defined I and J for Ref
                A,B = np.shape(np.array(data['Zdr']))
                Zdr_nocorr = np.array(data['Zdr'])
                x_miss = I-A
                y_miss = J-B
                if y_miss != 0:
                    adds = np.ones([I,y_miss])*np.nan
                    Zdr_corr = np.hstack([Zdr_nocorr,adds])
                    Zdr_nocorr = Zdr_corr
                #for m in range((0,x_miss)):
                    if x_miss !=0:
                        add = np.ones([x_miss,J])*np.nan
                        Zdr_corr = np.vstack([Zdr_nocorr,add])
                        Zdr_corr[Zdr_corr<-99900]=np.nan
                zdr_dict['data'] = Zdr_corr
                rad.add_field('ZDR',zdr_dict)

             
     '''folder='/Volumes/MiniHorus/Data/'+radar_name+'_'+date+'/AliasedVelocity/'+str(lvl)
     os.chdir(folder)
     l3=np.sort([f for f in os.listdir('.') if re.search('.netcdf', f)])
     file=l3[idz]
     
     #with gzip.open(file) as gz:
     with netCDF4.Dataset(file,mode='r') as data:
            vel_dict=get_metadata('Aliased Velocity')
            vel_dict['data']=np.array(data['AliasedVelocity'])
            try:
                rad.add_field('Aliased Velocity',vel_dict)
            except:
                #print('Velocity Failed, Now Correcting')
                # Note, we defined I and J for Ref
                A,B = np.shape(np.array(data['AliasedVelocity']))
                vel_nocorr = np.array(data['AliasedVelocity'])
                x_miss = I-A
                y_miss = J-B
                if y_miss != 0:
                    #print('Whoops')
                    adds = np.ones([I,y_miss])*np.nan
                    vel_corr = np.hstack([vel_nocorr,adds])
                    vel_nocorr = vel_corr
                #for m in range((0,x_miss)):
                    if x_miss != 0:
                        add = np.ones([x_miss,J])*np.nan
                        vel_corr = np.vstack([vel_nocorr,add])
                vel_dict['data'] = vel_corr
                rad.add_field('Aliased Velocity',vel_dict)'''
     
     folder='/Volumes/MiniHorus/Data/'+radar_name+'_'+dateKTLX+'/AliasedVelocity/'+str(lvl)
     os.chdir(folder)
     l4=np.sort([f for f in os.listdir('.') if re.search('.netcdf', f)])
     #print(l4)
     file=l4[idz]
     #print(file)
     
     with netCDF4.Dataset(file,mode='r') as data:
         vel_dict=get_metadata('Velocity')
         vel_dict['data']=np.array(data['AliasedVelocity'])
         #print(vel_dict['data'].shape)
         rad.add_field('Velocity',vel_dict)
         
     return rad

def read_Horus(idz):
    file_path=idz
    #print(file_path)
    filename=os.path.basename(file_path)
    #print(filename)
    time_str=filename[9:15]
    time_initial = datetime.strptime(time_str, "%H%M%S")
    #print('Initial Time: '+str(time_initial))
    with netCDF4.Dataset(file_path,mode='r') as data:
            #print(data)
            ele = np.linspace(0,30,31)
            
            #print(data.dimensions)
            #netCDF4.renameDimension(data,'Azimuth','Elevation')
            try:
                Ref=np.array(data['Reflectivity'])
            except: 
                Ref = np.array(data['Reflectivity_EXTRA']) 
            I,J=Ref.shape
            #Az=np.array(nc['Azimuth'])

            # Fill the Radar Object with important fields
            rad = pyart.testing.make_empty_ppi_radar(J, I, 1)
            time_dict=get_metadata('Time')
            if epoch == 1: 
                time_1 = datetime.fromtimestamp(int(np.array([data.Time])))
                time_1 = time_1 - timedelta(hours=5, minutes=0)
                time_dict['data'] = np.array([int(time_1.timestamp())])
            else:
                time_dict['data']=np.array([data.Time])
            #print(data.Time)
            #sys.exit()
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
            rad.fixed_angle['data']=np.array([data.Elevation])
            rad.elevation['data'] = np.array([data.Elevation])
            #rad.sweep_mode['data'] = np.array(data['sweep_mode'])
            rad.instrument_name='Horus'
            rad.metadata['instrument_name']='Horus'
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
            ref_dict = get_metadata('Reflectivity')
            try:
                ref_dict['data'] = np.array(data['Reflectivity'])
            except:
                ref_dict['data'] = np.array(data['Reflectivity_EXTRA'])
            ref_dict['units']='dBZ'
            rad.add_field('Reflectivity', ref_dict)
            
            #RhoHV
            os.chdir(Rhofolder)
            l2=np.sort([f for f in os.listdir('.') if re.search('.netcdf', f)])
            file_path=idz
            with netCDF4.Dataset(file_path,mode='r') as data:
                rho_dict = get_metadata('RhoHV')
                try:
                    rho_dict['data'] = np.array(data['RhoHV'])
                except:
                    rho_dict['data'] = np.array(data['RhoHV_EXTRA'])
                rho_dict['units'] = 'unitless'
                rad.add_field('RhoHV',rho_dict)
            
            #ZDR
            os.chdir(Zdrfolder)
            l2=np.sort([f for f in os.listdir('.') if re.search('.netcdf', f)])
            file_path=idz
            with netCDF4.Dataset(file_path,mode='r') as data:
                zdr_dict = get_metadata('Zdr')
                try:
                    zdr_dict['data'] = np.array(data['Zdr'])
                except:
                    zdr_dict['data'] = np.array(data['Zdr_EXTRA'])  
                zdr_dict['units'] = 'dB'
                rad.add_field('ZDR',zdr_dict)
                
            # Velocity
            os.chdir(Velfolder)
            l2=np.sort([f for f in os.listdir('.') if re.search('.netcdf', f)])
            file_path=idz
            with netCDF4.Dataset(file_path,mode='r') as data:
                zdr_dict = get_metadata('velocity')
                try:
                    zdr_dict['data'] = np.array(data['AliasedVelocity'])
                except:
                    zdr_dict['data'] = np.array(data['AliasedVelocity_EXTRA'])
                zdr_dict['units'] = 'm/s'
                rad.add_field('velocity',zdr_dict)
                
            # HCA
            os.chdir(Hcafolder)
            l2=np.sort([f for f in os.listdir('.') if re.search('.netcdf', f)])
            file_path=l2[t]
            with netCDF4.Dataset(file_path,mode='r') as data:
                zdr_dict = get_metadata('HCA')
                zdr_dict['data'] = np.array(data['DHCA'])
                zdr_dict['units'] = 'unitless'
                rad.add_field('HCA',zdr_dict)
    return rad
# Sets up a rhohv colormap to be used
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

# The elevation angles on Horus are all different so this allows for indexing 
# into the folders based on the specified date
if dateKTLX == '20230511' or dateKTLX == '20231004' or dateKTLX == '20230808':
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
elif dateKTLX == '20230216':
        lvls1 = ['00.00','01.33','02.67','04.00','05.33','06.67','08.00','09.33','10.67',
                 '12.00','13.33','14.67','16.00','17.33','18.67']
        eles = [00.00,01.33,02.67,04.00,05.33,06.67,08.00,09.33,10.67,
                 12.00,13.33,14.67,16.00,17.33,18.67]
elif dateKTLX == '20221213':
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

#Horus Folders
Reffolder = '/Volumes/MiniHorus/Data/Horus_'+date+'/Reflectivity/'+lvls1[0]
Rhofolder = '/Volumes/MiniHorus/Data/Horus_'+date+'/output/RhoHV/'+lvls1[0]
Zdrfolder = '/Volumes/MiniHorus/Data/Horus_'+date+'/output/Zdr/'+lvls1[0]
Velfolder = '/Volumes/MiniHorus/Data/Horus_'+date+'/output/AliasedVelocity/'+lvls1[0]
#Meshfolder = '/Volumes/MiniHorus/Data/Horus_20230511/MESH/at'+dhca_lvl
Hcafolder = '/Volumes/MiniHorus/Data/Horus_'+date+'/output/DHCA360/'+lvls1[0]
Kdpfolder = '/Volumes/MiniHorus/Data/Horus_'+date+'/output/DKDP/'+lvls1[0]
rzdrfolder = '/Volumes/MiniHorus/Data/Horus_'+date+'/output/RZDR/'+lvls1[0]

# Sets up the HCA colormaps to be used later
ticks = np.arange(1,15)
labs = ['L/MR','HR','R/Ha','BD','AP','BI','UK','NE','DS','WS','CY','','GR','+']
hca_cmap = colors.ListedColormap(['#AAAAAA','#0000FE','#0057FF','#00A9FF','#00FFFF','#56FFAA','#ACFD56','#FFFF00','#FEAB02','#FF5700','#FF00FF','#B300B3','#DC143C','#FF00FF'])
hca_cmap = colors.ListedColormap(['limegreen','green','#DC143C','#FFFF00','dimgray','lightgray','#B300B3','white','#00FFFF','#0000FE','mediumturquoise','#000000','#FF00FF','#000000'])

# Here is the stuff we want now
clicking = 'True'

levels_KTLX = ['00.90','01.30','01.80','02.40','03.10','04.00','05.10','06.40','08.00','10.00','12.50','15.60','19.50']
lvls1_KTLX = ['00.90','01.30','01.80','02.40','03.10','04.00','05.10','06.40','08.00','10.00','12.50','15.60','19.50']

# Set Up Some User Specified Parameters
os.chdir(Reffolder)
l_Hrus_index = np.sort([f for f in os.listdir('.') if re.search('.netcdf', f)])
for j in range(0,len(l_Hrus_index)):
    print(str(j)+'. '+l_Hrus_index[j])
timet = int(input('What Time Do You Want To Start With? '))
for n in range(0,len(lvls1_KTLX)):
    print(str(n)+'. '+lvls1_KTLX[n])
eles = int(input('What Elevation Angle Do you Want to Use as a Baseline? Use Number. '))

#while clicking=='True':
    #ask = input('Do You Want To Plot? [y:1, n:0] ')
    #ask=1
    
    #if ask ==1:
for t in range(timet,len(l_Hrus_index)):
        os.chdir(Reffolder)
        l_Hrus = np.sort([f for f in os.listdir('.') if re.search('.netcdf', f)])
        Hrus_file = l_Hrus[t]
        print(Hrus_file)
        Hrusrads = []  
        Hrusrads2 =[]
        
        for j in range(0,len(lvls1)):
            Reffolder = '/Volumes/MiniHorus/Data/Horus_'+date+'/Reflectivity/'+lvls1[j]
            Rhofolder = '/Volumes/MiniHorus/Data/Horus_'+date+'/RhoHV/'+lvls1[j]
            Zdrfolder = '/Volumes/MiniHorus/Data/Horus_'+date+'/Zdr/'+lvls1[j]
            Velfolder = '/Volumes/MiniHorus/Data/Horus_'+date+'/AliasedVelocity/'+lvls1[j]
            #Meshfolder = '/Volumes/MiniHorus/Data/Horus_20230511/MESH/at'+lvls1[i]
            Hcafolder = '/Volumes/MiniHorus/Data/Horus_'+date+'/DHCA/'+lvls1[j]
            Kdpfolder = '/Volumes/MiniHorus/Data/Horus_'+date+'/DKDP/'+lvls1[j]
    
            os.chdir(Reffolder)
    
            Hrusrad = read_Horus(Hrus_file)
    
            if j == 0:
                base_Hrus = Hrusrad
                Hrs_time = int(Hrus_file[9:15])
                Hrs_time_str = Hrus_file[9:15]
    
            Hrusrads.append(Hrusrad)
            Hrusrads2.append(Hrusrad)
    
        while len(Hrusrads) > 1:
            bob = Hrusrads.pop(1)
        #print(bob.nrays)
            Hrusrads[0] = pyart.util.join_radar(Hrusrads[0], bob)
        
        Horusradar=Hrusrads[0]
        # Fix the Azimuths (11 MAY SPECIFIC)
        if dateKTLX == '20230511':
            az_len = len(Horusradar.azimuth['data'])
            if t >= 1510:
                Horusradar.azimuth['data'] = np.ones(az_len) * 285
                base_Hrus.azimuth['data'] = np.array([285])
            else:
                Horusradar.azimuth['data'] = np.ones(az_len) * 270
                base_Hrus.azimuth['data'] = np.array([270])
        
        # Fix the Azimuths for 13 Dec Case
        if dateKTLX == '20221213':
            az_len = len(Horusradar.azimuth['data'])
            Horusradar.azimuth['data'] = np.ones(az_len) * 60
            base_Hrus.azimuth['data'] = np.array([60])
        elif dateKTLX == '20231004':
            az_len = len(Horusradar.azimuth['data'])
            Horusradar.azimuth['data'] = np.ones(az_len) * 60
            base_Hrus.azimuth['data'] = np.array([60])
        #elif dateKTLX == '20230808':
        #    Horusradar.azimuth['data'] = Horusradar.azimuth['data']+30
        #    base_Hrus.azimuth['data'] = base_Hrus.azimuth['data']+30

        Horuslat = base_Hrus.latitude['data']
        Horuslon = base_Hrus.longitude['data']
        
        # Get the Points
        lats,lons,alts = base_Hrus.get_gate_lat_lon_alt(0)
        lats = lats[0]
        lons = lons[0]
        
        #######################
        # Now for the KTLX Data
        #########################
        filepath = '/Volumes/MiniHorus/Data/'+radar_name+'_'+dateKTLX
        Reffolder_KTLX = filepath + '/Reflectivity/00.90'
        os.chdir(Reffolder_KTLX)
        l_KTLX = np.sort([f for f in os.listdir('.') if re.search('.netcdf', f)])
        break_flag = 0
        for m in range(0,len(l_KTLX)):
            if break_flag == 1:
                break
            KTLX_time = l_KTLX[m][9:15]
            if int(KTLX_time) < 100000:
                if int(Hrs_time_str[0:2]) == 0:
                    time_diff_idx = (int(Hrs_time)+240000) - (int(KTLX_time)+240000)
                else:
                    time_diff_idx = int(Hrs_time) - (int(KTLX_time)+240000)
            else:
                if int(Hrs_time_str[0:2]) == 0:
                    time_diff_idx = (int(Hrs_time)+240000) - int(KTLX_time)
                else: 
                    time_diff_idx = int(Hrs_time) - int(KTLX_time)
            #print(time_diff_idx)
            if time_diff_idx <= 0:
                break_flag=1
                time_idx = m 
                print(KTLX_time)
                #print(l_KTLX[time_idx])
                #print(Hrs_time)
        #sys.exit()
            
        time_initial,time_save, indices_lvls,lvls = read_KOUN(time_idx)
        rads = []
        rads_dp= []
        KTLXrads2=[]
        indices_lvls=time_idx
        # Skip 1, Plus 11 for next time step
        for n in range(1,len(levels_KTLX)):
            #os.chdir('/Volumes/MiniHorus/Data')
            folder='/Volumes/MiniHorus/Data/'+radar_name+'_'+dateKTLX+'/Reflectivity/'+str(levels_KTLX[n])
            os.chdir(folder)
            l2 = np.sort([f for f in os.listdir('.') if re.search('.netcdf', f)])
            for j in range(0,len(l2)):
                    idz1=j
                    file=l2[idz1]
                    filename2=os.path.basename(file)
                    time_str1=filename2[9:15]
                    #print(time_str1)
                    date_str1=filename2[0:8]
                    time_cat1=date_str1+time_str1
                    time_aloft=datetime.strptime(time_cat1, "%Y%m%d%H%M%S")
                    time_diff = time_aloft - time_initial
                    if time_diff.total_seconds() < 0:
                            continue
                    elif time_diff.total_seconds() < 450 or time_diff.total_seconds==0:
                            indices_lvls = np.append(indices_lvls,j)
                            #filenames = np.append(filenames,filename2)
                            lvls = np.append(lvls,levels_KTLX[n])
                            #time_save = time_initial
                            #print(time_str1)
                    else:
                            continue
                        
        #start=0
        start_new=0
        start_new_dp=0
        #end_new=0
        start_rad=[]
        end_rad=[]
        #print(time_save)
        
        for k in range(0,len(lvls)):
                #idx3 = int(indices_lvls[k])
                #print(k)
                idx3 = indices_lvls
                Reffolder_KTLX = '/Volumes/MiniHorus/Data/'+radar_name+'_'+dateKTLX+'/Reflectivity/'+str(lvls[k])
                Zdrfolder_KTLX = '/Volumes/MiniHorus/Data/'+radar_name+'_'+dateKTLX+'/Zdr/'+str(lvls[k])
                os.chdir(Reffolder_KTLX)
                l = np.sort([f for f in os.listdir('.') if re.search('.netcdf', f)])
                file=l[idx3]
                file = file[0]

                #with gzip.open(file) as gz:
                with netCDF4.Dataset(file,mode='r') as data:
        
                        Ref=np.array(data['Reflectivity'])
                        I,J=Ref.shape
                
                os.chdir(Zdrfolder)
                l_KTLX1 = np.sort([f for f in os.listdir('.') if re.search('.netcdf', f)])
                file=l_KTLX1[idx3]
                file = file[0]

                #with gzip.open(file) as gz:
                with netCDF4.Dataset(file,mode='r') as data:
        
                        Zdr=np.array(data['Zdr'])
                        I2,J2=Ref.shape

                if k == 0:
                    start2=0
                    start_rad=int(start2)
                    start_rad_dp = int(start2)
                    end=I-1
                    end_dp = I2 - 1
                    end_rad=int(end)
                    end_rad_dp = int(end_rad)
                    rad2 = read_KOUN_lvls(lvls1_KTLX[eles],idx3[0],start2,end)
                    #rad2 = pyart.io.read(file)
                    #sys.exit()
                    rad3 = read_KOUN_lvls_dp(lvls1_KTLX[eles],idx3[0],start2,end)
                else:
                    start2=end+1
                    start2_dp = end_dp+1
                    #print(start2)
                    start_rad=np.append(start_rad,int(start2))
                    start_rad_dp = np.append(start_rad_dp,int(start2_dp))
                    end=end+I
                    end_dp = end_dp + I2
                    end_rad=np.append(end_rad,int(end))
                    end_rad_dp = np.append(end_rad_dp,int(end_dp))

               # print(start,end)
                rad = read_KOUN_lvls(lvls[k],idx3[k],start2,end)
                rad_dp = read_KOUN_lvls_dp(lvls[k],idx3[k],start2,end)
                start_new = len(rad.fields['Reflectivity']['data'][:,0])
                start_new_dp = len(rad_dp.fields['ZDR']['data'][:,0])
                #end_new=(len(rad.fields['reflectivity_horizontal']['data'][:,0])-1)
           # rad = read_KOUN_lvls(levels[n],timet)
                rads.append(rad)
                rads_dp.append(rad_dp)
                KTLXrads2 = np.append(KTLXrads2,rad)
        
        #Concatenate all of the individual Radar objects into one
        while len(rads) > 1:
            #print('Hiding in Here')
            bob = rads.pop(1)
            #print(bob.nrays)
            rads[0] = pyart.util.join_radar(rads[0], bob)
            
        while len(rads_dp) > 1:
            #print('Hiding in Here')
            bob = rads_dp.pop(1)
            #print(bob.nrays)
            rads_dp[0] = pyart.util.join_radar(rads_dp[0], bob)
        
        rad=rads[0]
        rad_dp = rads_dp[0]
        
        
        rad_len = end_rad - start_rad + 1
        rad_len_dp = end_rad_dp - start_rad_dp + 1
        rad.sweep_start_ray_index['data']=start_rad
        rad_dp.sweep_start_ray_index['data'] = start_rad_dp
        rad.sweep_end_ray_index['data']=end_rad
        rad_dp.sweep_end_ray_index['data']=end_rad_dp
        rad.rays_per_sweep['data']=rad_len
        rad_dp.rays_per_sweep['data'] = rad_len_dp
        p,q = rad.fields['Reflectivity']['data'].shape
        p2,q2 = rad_dp.fields['ZDR']['data'].shape
        rad.nrays = p
        rad_dp.nrays = p2
        
        #clicking = 'False'
        
        ####### MPING API REQUEST #########
        # This next section can be commented out if you do not want/care for the mping reports 
        ###################################
        if t == timet:
            reqheaders = {
                'content-type':'application/json',
                'Authorization': 'Token 526911e6a3b40b1ed9d74b1b6bd4b73d5f72ef01'
                }
            
            reqparams = {
                'category':'Hail',
                'year': '2023',
                'month': '08',
                'day': '08'
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
        
        import cartopy.crs as ccrs
        #######################
        ##### NEW FIGURE ######
        #######################
        
        os.chdir(figfolder)
        
        
        # This Figure just plots a PPI map of KTLX data along with the lat/lon from the Horus file
        '''display = pyart.graph.RadarMapDisplay(rad2)
        projection = ccrs.LambertConformal(
            central_latitude=rad2.latitude["data"][0],
            central_longitude=rad2.longitude["data"][0],
            )
        figPPI=plt.figure(figsize=[10,6])
        ax1 = figPPI.add_subplot(111,projection=ccrs.LambertConformal())
        msk=display.plot_ppi_map('DHCA',ax=ax1,cmap=hca_cmap,vmin=0.5, vmax=14.5,
                                 projection=projection,min_lat = y_lims[0],max_lat=y_lims[1],
                                 min_lon=x_lims[0],max_lon=x_lims[1])
        ax1.scatter(Horuslon,Horuslat,color='k',marker='*',s=100,transform=ccrs.Geodetic())
        # If you remove the API request, the next line needs to be commented out b/c code will break
        ax1.plot(lons,lats,color='k',transform=ccrs.Geodetic())
        props = dict(boxstyle='round', facecolor='white', alpha=0.8)
        ax1.text(-97.4,35.22,'Horus',transform = ccrs.Geodetic(),fontsize=14,bbox=props)
        #now plot reports
        ax1.scatter(lon,lat,transform=ccrs.Geodetic(),color='k')
        #plt.xlim(x_lims[0],x_lims[1])
        #plt.ylim(y_lims[0],y_lims[1])
        #plt.savefig(date_str1+time_save+'_PPI_'+str(int(angle))+'_Storm'+storm_num,bbox_inches = 'tight')
        #figPPI.show()
        plt.close()'''
        
        ########################################################################
        # Now, we need to extract the data points from KTLX (bc Horus is an RHI)  
        ########################################################################
        
        KTLXlats,KTLXlons,KTLXalts = rad2.get_gate_lat_lon_alt(0)
        KTLXlats= KTLXlats.flatten()
        KTLXlons = KTLXlons.flatten()
        
        Horusloc = np.vstack([lats,lons]).T
        KTLXloc = np.vstack([KTLXlats,KTLXlons]).T
        
        p = path.Path(Horusloc)
        # Here the radius gets adjusted for a level of tolerance on matching points
        ind = p.contains_points(KTLXloc, radius=0.004)
        I = np.shape(KTLXlats.flatten())
        array = np.zeros((I))
        
        def updateArray(array, indices):
            lin = np.arange(array.size)
            newArray = array.flatten()
            newArray[indices] = 1
            return newArray.reshape(array.shape)
        
        array = updateArray(array, ind)
        arrayy=KTLXloc[ind]
        arrayy2 = arrayy
        arrayy3 = arrayy[arrayy[:,1].argsort()]
        
        arrayy = []
        arrayy = arrayy3
        distances = []
        
        for r in range(0,len(arrayy[:,1])):
            if r == 0:
                dist = 0
                distances = np.append(distances,dist)
            else:
                dist = geopy.distance.geodesic(arrayy[0,:],arrayy[r,:]).km
                distances = np.append(distances,dist)
        
        #######################
        ###### NEW FIGURE #####
        #######################
        # Plot the extracted points on a KTLX PPI to make sure they are where you want them to be
        # SO, the dots here are KTLX LAT/LON NOT HORUS
        # if they aren't or if the box is too big... go back and check the radius being used
        
        display = pyart.graph.RadarMapDisplay(rad2)
        projection = ccrs.LambertConformal(
            central_latitude=rad2.latitude["data"][0],
            central_longitude=rad2.longitude["data"][0],
            )
        figPPI=plt.figure(figsize=[10,6])
        ax1 = figPPI.add_subplot(111,projection=ccrs.LambertConformal())
        msk=display.plot_ppi_map('Reflectivity',ax=ax1,cmap='pyart_NWSRef',vmin=-10, vmax=70,
                                 projection=projection,min_lat = y_lims[0],max_lat=y_lims[1],
                                 min_lon=x_lims[0],max_lon=x_lims[1])
        ax1.scatter(Horuslon,Horuslat,color='k',marker='*',s=100,transform=ccrs.Geodetic())
        ax1.scatter(arrayy[:,1],arrayy[:,0],color='blue',transform=ccrs.Geodetic(),s=10)
        props = dict(boxstyle='round', facecolor='white', alpha=0.8)
        ax1.text(-97.4,35.22,'Horus',transform = ccrs.Geodetic(),fontsize=14,bbox=props)
        plt.savefig(dateKTLX+'_'+str(Hrs_time)+'_Ref_PPI',bbox_inches = 'tight')
        #figPPI.show()
        plt.close()
        
        #########################################
        ##### KTLX DATA EXTRACTION METHOD 1 #####
        #########################################
        # The most effective means for gridding is by using the raw (V06) files
        os.chdir('/Volumes/MiniHorus/Data/KTLX_'+dateKTLX+'_Raw')
        files = np.sort([f for f in os.listdir('.') if re.search('V06', f)])
        break_flag = 0
        for m in range(0,len(files)):
            if break_flag == 1:
                break
            KTLX_time = files[m][13:19]
            if int(Hrs_time_str[0:2]) == 0:
                if int(KTLX_time[0:2]) == 0:
                    time_diff_idx = (int(Hrs_time)+240000) - (int(KTLX_time)+240000)
                else:
                    time_diff_idx = (int(Hrs_time)+240000) - int(KTLX_time)
            else:
                time_diff_idx = int(Hrs_time) - int(KTLX_time)
            #print(time_diff_idx)
            if time_diff_idx <= 0:
                break_flag=1
                time_idx = m - 1
        radartest = pyart.io.read(files[time_idx])
        
        ### First, we try gridding with the V06 Files ###
        grid = pyart.map.grid_from_radars(
            (radartest),
            grid_shape=(len(levels_KTLX), 181, 181),
            grid_limits=grid_lims,
            grid_origin=(rad.latitude['data'][0], rad.longitude['data'][0]),
            fields=["reflectivity","velocity"],
        )
        
        # Now we try gridding with the radar objects that were read in with read_KOUN()
        # This one uses data from read_KOUN()
        '''grid2 = pyart.map.grid_from_radars(
            (rad),
            grid_shape=(len(lvls), 500, 500),
            grid_limits=grid_lims,
            grid_origin=(rad.latitude['data'][0], rad.longitude['data'][0]),
            fields=["DHCA",'Reflectivity'],
        )'''
        
        # This one uses data from read_KOUN_dp()
        '''grid3 = pyart.map.grid_from_radars(
            (rad_dp),
            grid_shape=(len(lvls), 500, 500),
            grid_limits=grid_lims,
            grid_origin=(rad3.latitude['data'][0], rad3.longitude['data'][0]),
            fields=["RhoHV",'ZDR',"Velocity"],
        )'''
        
        ##################
        ##### FIGURE #####
        ##################
        
        os.chdir(figfolder)
        
        # This plots the V06 Gridded KTLX data in a PPI 
        # and the extracted cross section using the Horus start and end
        start=arrayy[0,:]
        end=arrayy[-1,:]
         ### FIGURE 1 
        fig = plt.figure(figsize=(18, 6))
        ax1 = plt.subplot(121, projection=ccrs.PlateCarree())
        display = pyart.graph.GridMapDisplay(grid)
        #display_grid2 = pyart.graph.GridMapDisplay(grid2)
        #display_grid3 = pyart.graph.GridMapDisplay(grid3)
        display.plot_grid(
            "reflectivity",
            ax=ax1,
            cmap="pyart_NWSRef",
            vmin=-10,
            vmax=70,
        )
        
        # Plot our start and end points, as well as a line in between the two
        ax1.scatter(start[1], start[0], color="tab:blue", label="Start")
        ax1.scatter(end[1], end[0], color="black", label="End")
        ax1.plot([start[1], end[1]], [start[0], end[0]], color="k", linestyle=":")
        
        #Now Plot Horus Data Points (not the actual data, just location) on Here Too
        ax1.scatter(Horuslon,Horuslat,color='k',marker='*',s=100,transform=ccrs.Geodetic())
        ax1.plot(lons,lats,color='k',transform=ccrs.Geodetic())
        plt.legend(loc="upper right")
        
        # Add a cross section, using our start and end points, and set our x-axis as latitude (lat)
        ax2 = plt.subplot(122)
        display.plot_cross_section(
            "reflectivity",
            start,
            end,
            x_axis="lon",
            cmap="pyart_NWSRef",
            vmin=-10,
            vmax=70,
        )
        plt.savefig(dateKTLX+'_'+str(Hrs_time)+'_Ref_PPI_RHI_V06',bbox_inches = 'tight')
        plt.close()
        
        #########################
        ###### NEXT FIGURE ######
        #########################
        # Now, do the same thing but with the radar object gridding
        '''
        fig = plt.figure(figsize=(18, 6))
        ax1 = plt.subplot(121, projection=ccrs.PlateCarree())
        display3 = pyart.graph.GridMapDisplay(grid2)
        display3.plot_grid(
            "DHCA",
            ax=ax1,
            cmap=hca_cmap,
            vmin=0.5,
            vmax=14.5,
        )
        # Plot our start and end points, as well as a line in between the two
        ax1.scatter(start[1], start[0], color="tab:blue", label="Start")
        ax1.scatter(end[1], end[0], color="black", label="End")
        ax1.plot([start[1], end[1]], [start[0], end[0]], color="k", linestyle=":")
        
        #Now Plot Horus on Here Too
        ax1.scatter(Horuslon,Horuslat,color='k',marker='*',s=100,transform=ccrs.Geodetic())
        ax1.plot(lons,lats,color='k',transform=ccrs.Geodetic())
        
        plt.legend(loc="upper right")
        
        # Add a cross section, using our start and end points, and set our x-axis as latitude (lat)
        ax2 = plt.subplot(122)
        display3.plot_cross_section(
            "DHCA",
            start,
            end,
            x_axis="lon",
            cmap=hca_cmap,
            vmin=0.5,
            vmax=14.5,
        )
        '''
        #Horusradar.range['data'] = Horusradar.range['data']*-1
        lon_plot = Horusradar.gate_longitude['data']
        lon_plot = lon_plot[0,:]
        xticks = ax2.get_xticks()
        x,y,z= Horusradar.get_gate_x_y_z(0)
        x = x/1000
        x_tick = []
        x_tick_lab = []
        x_tick_ktlx = []
        x_tick_lab_ktlx = []
        
        # Now for some later on Horus plots, we will need the xticks which is what this loop does
        for i in range(0,len(x[0,:])):
            if i%250 == 0:
                #print(i)
                x_tick = np.append(x_tick,x[0,i])
                x_tick_lab=np.append(x_tick_lab,str(round(lon_plot[i],2)))
                
        for i in range(0,len(arrayy[:,0])):
            if i%50 == 0:
                #print(i)
                x_tick_ktlx = np.append(x_tick_ktlx,distances[i])
                x_tick_lab_ktlx=np.append(x_tick_lab_ktlx,str(round(arrayy[i,1],2)))
        
        #plt.savefig(dateKTLX+'_'+str(Hrs_time)+'_DHCA_PPI_RHI',bbox_inches = 'tight')
        #plt.close()'''
        
        ######## NEXT FIGURE #########
        # This one just plots the RHI of Horus Data 
        #############################
        xsect=pyart.util.cross_section_ppi(Horusradar,[Horusradar.azimuth['data'][0]])
        display2=pyart.graph.RadarMapDisplay(xsect)
        #fig=plt.figure()
        #display2.plot_rhi('Reflectivity',vmin=20,vmax=70,cmap='pyart_NWSRef',reverse_xaxis=False)
        #plt.close()
        
        ################################################
        # Joint Figures of the KTLX Gridded Cross Section
        # Along with the Horus Data
        ################################################
        aspect2 = 1
        fig = plt.figure(figsize=(16,12))
        ax1 = plt.subplot(121)
        pcm = display.plot_cross_section(
            "reflectivity",
            start,
            end,
            x_axis="lon",
            cmap="pyart_NWSRef",
            vmin=-10,
            vmax=70,
            title= radartest.time['units'][14:] + '\n KTLX Reflectivity Cross Section',
            colorbar_flag=False
        )
        ax1.set_ylim(cross_ylim)
        ax1.set_aspect(aspect2/100)
        
        ax2 = plt.subplot(122)
        if Horusradar.azimuth['data'][0] > 271:
            display2.plot_rhi('Reflectivity',vmin=-10,vmax=70,cmap='pyart_NWSRef',reverse_xaxis=True,colorbar_flag=False)
        else:
            display2.plot_rhi('Reflectivity',vmin=-10,vmax=70,cmap='pyart_NWSRef',reverse_xaxis=False,colorbar_flag=False)
        ax2.set_ylim(cross_ylim)
        ax2.set_xticks(x_tick)
        ax2.set_xticklabels(x_tick_lab)
        ax2.set_aspect(aspect2)
        
        cax = ax1.inset_axes([-80, -6, 50, 1.25],transform=ax2.transData)
        display2.plot_colorbar(pcm,cax=cax,orient='horizontal')
        plt.savefig(dateKTLX+'_'+str(Hrs_time)+'_Ref_Comp_Gridded',bbox_inches = 'tight')
        plt.close()
        
        fig = plt.figure(figsize=(16,12))
        ax1 = plt.subplot(121)
        pcm = display.plot_cross_section(
            "velocity",
            start,
            end,
            x_axis="lon",
            cmap="pyart_BuDRd18",
            vmin=-30,
            vmax=30,
            title= radartest.time['units'][14:] + '\n KTLX Velocity Cross Section',
            colorbar_flag=False
        )
        ax1.set_ylim(cross_ylim)
        ax1.set_aspect(aspect2/100)
        
        ax2 = plt.subplot(122)
        if Horusradar.azimuth['data'][0] > 271:
            display2.plot_rhi('velocity',vmin=-30,vmax=30,cmap='pyart_BuDRd18',reverse_xaxis=True,colorbar_flag=False)
        else:
            display2.plot_rhi('velocity',vmin=-30,vmax=30,cmap='pyart_BuDRd18',reverse_xaxis=False,colorbar_flag=False)
        ax2.set_ylim(cross_ylim)
        ax2.set_xticks(x_tick)
        ax2.set_xticklabels(x_tick_lab)
        ax2.set_aspect(aspect2)
        
        cax = ax1.inset_axes([-80, -6, 50, 1.25],transform=ax2.transData)
        display2.plot_colorbar(pcm,cax=cax,orient='horizontal')
        plt.savefig(dateKTLX+'_'+str(Hrs_time)+'_Velocity_Comp_Gridded',bbox_inches = 'tight')
        plt.close()
        '''
        fig = plt.figure(figsize=(18,7))
        ax1 = plt.subplot(121)
        pcm = display_grid2.plot_cross_section(
            "DHCA",
            start,
            end,
            x_axis="lon",
            cmap=hca_cmap,
            vmin = 0.5,
            vmax=14.5,
            title= radartest.time['units'][14:] + '\n KTLX DHCA Cross Section',
            colorbar_flag=False
        )
        ax1.set_ylim(cross_ylim)
        #ax1.set_aspect(aspect)
        
        ax2 = plt.subplot(122)
        if Horusradar.azimuth['data'][0] > 271:
            display2.plot_rhi('HCA',vmin=0.5,vmax=14.5,cmap=hca_cmap,reverse_xaxis=True)
        else:
            display2.plot_rhi('HCA',vmin=0.5,vmax=14.5,cmap=hca_cmap,reverse_xaxis=False)
        ax2.set_ylim(cross_ylim)
        ax2.set_xticks(x_tick)
        ax2.set_xticklabels(x_tick_lab)
        #ax2.set_aspect(aspect)
        
        #cax = ax1.inset_axes([17, -6, 50, 1.6],transform=ax1.transData)
        #fig.colorbar(pcm,cax=cax,orientation='horizontal')
        plt.savefig(dateKTLX+'_'+str(Hrs_time)+'_DHCA_Comp_Gridded',bbox_inches = 'tight')
        plt.close()'''
        '''
        fig = plt.figure(figsize=(18,7))
        ax1 = plt.subplot(121)
        pcm = display_grid3.plot_cross_section(
            "ZDR",
            start,
            end,
            x_axis="lon",
            cmap="pyart_HomeyerRainbow",
            vmin=-1,
            vmax=70,
            title= radartest.time['units'][14:] + '\n KTLX ZDR Cross Section',
            colorbar_flag=False
        )
        ax1.set_ylim(cross_ylim)
        #ax1.set_aspect(aspect)
        
        ax2 = plt.subplot(122)
        if Horusradar.azimuth['data'][0] > 271:
            display2.plot_rhi('ZDR',vmin=-1,vmax=7,cmap='pyart_HomeyerRainbow',reverse_xaxis=True)
        else:
            display2.plot_rhi('ZDR',vmin=-1,vmax=7,cmap='pyart_HomeyerRainbow',reverse_xaxis=False)
        ax2.set_ylim(cross_ylim)
        ax2.set_xticks(x_tick)
        ax2.set_xticklabels(x_tick_lab)
        #ax2.set_aspect(aspect)
        
        #cax = ax1.inset_axes([17, -6, 50, 1.6],transform=ax1.transData)
        #fig.colorbar(pcm,cax=cax,orientation='horizontal')
        plt.savefig(dateKTLX+'_'+str(Hrs_time)+'_ZDR_Comp_Gridded',bbox_inches = 'tight')
        plt.close()
        
        fig = plt.figure(figsize=(18,7))
        ax1 = plt.subplot(121)
        pcm = display_grid3.plot_cross_section(
            "RhoHV",
            start,
            end,
            x_axis="lon",
            cmap=cc_cmap,
            vmin=0.5,
            vmax=1,
            title= radartest.time['units'][14:] + '\n KTLX RhoHV Cross Section',
            colorbar_flag=False
        )
        ax1.set_ylim(cross_ylim)
        #ax1.set_aspect(aspect)
        
        ax2 = plt.subplot(122)
        if Horusradar.azimuth['data'][0] > 271:
            display2.plot_rhi('RhoHV',vmin=0.5,vmax=1,cmap=cc_cmap,reverse_xaxis=True)
        else:
            display2.plot_rhi('RhoHV',vmin=0.5,vmax=1,cmap=cc_cmap,reverse_xaxis=False)
        ax2.set_ylim(cross_ylim)
        ax2.set_xticks(x_tick)
        ax2.set_xticklabels(x_tick_lab)
        #ax2.set_aspect(aspect)
        
        #cax = ax1.inset_axes([17, -6, 50, 1.6],transform=ax1.transData)
        #fig.colorbar(pcm,cax=cax,orientation='horizontal')
        plt.savefig(dateKTLX+'_'+str(Hrs_time)+'_RhoHV_Comp_Gridded',bbox_inches = 'tight')
        plt.close()'''
        #sys.exit()
        ########################################################################
        # Lets try another extraction method for the data using a pyart function
        ########################################################################
        for k in range(0,len(arrayy[:,0])):
            ds = pyart.util.columnsect.get_field_location(rad, arrayy[k,0], arrayy[k,1])
            ds_dp = pyart.util.columnsect.get_field_location(rad_dp,arrayy[k,0],arrayy[k,1])
            dhcas = ds['DHCA'].data
            refs = ds['Reflectivity'].data
            zdrs = ds_dp['ZDR'].data
            rhos = ds_dp['RhoHV'].data
            vels = ds_dp['Velocity'].data
            heights = ds['height'].data
            
            if k == 0:
                extracted_hca = dhcas
                extracted_ref = refs
                extracted_zdr = zdrs
                extracted_rho = rhos
                extracted_vel = vels
                height_extract = heights
            else:
                extracted_hca = np.column_stack([extracted_hca,dhcas])
                extracted_ref = np.column_stack([extracted_ref,refs])
                extracted_zdr = np.column_stack([extracted_zdr,zdrs])
                extracted_rho = np.column_stack([extracted_rho,rhos])
                extracted_vel = np.column_stack([extracted_vel,vels])
                height_extract = np.column_stack([height_extract,heights])
                
        # Need to Flip the Arrays because otherwise they plot upside down
        extracted_hca = np.flipud(extracted_hca)
        extracted_ref = np.flipud(extracted_ref)
        extracted_zdr = np.flipud(extracted_zdr)
        extracted_rho = np.flipud(extracted_rho)
        extracted_vel = np.flipud(extracted_vel)
        height_extract = np.flipud(height_extract)/1000
        
        #####################################################
        ###### Plots With The Column Extraction Method ######
        #####################################################
        xs = np.arange(0,len(height_extract[0,:]))
        xs = np.tile(xs,len(height_extract[:,0]))
        xs = np.reshape(xs,np.shape(height_extract))
        ### DHCA
        fig = plt.figure(figsize = (16,10))
        ax1 = plt.subplot(121)
        pcm = ax1.pcolormesh(distances,height_extract,extracted_hca,cmap=hca_cmap,vmin=0.5,vmax=14.5)
        ax1.set_ylim(cross_ylim)
        ax1.set_xlim(cross_xlim)
        ax1.set_xticks(x_tick_ktlx)
        ax1.set_xticklabels(x_tick_lab_ktlx)
        ax1.set_xticks(x_tick_ktlx)
        plt.title(radartest.time['units'][14:] + '\n KTLX DHCA Cross Section')
        ax1.set_aspect(aspect)
        ax1.set_ylabel('')
        ax1.set_xlabel('')
        
        ax2 = plt.subplot(122)
        if Horusradar.azimuth['data'][0] > 271:
            display2.plot_rhi('HCA',vmin=0.5,vmax=14.5,cmap=hca_cmap,reverse_xaxis=True,ticklabs=labs,ticks=ticks,colorbar_flag=False)
        else:
            display2.plot_rhi('HCA',vmin=0.5,vmax=14.5,cmap=hca_cmap,reverse_xaxis=False,ticklabs=labs,ticks=ticks,colorbar_flag=False)
        ax2.set_ylim(cross_ylim)
        ax2.set_xlim(cross_xlim)
        ax2.set_xticks(x_tick)
        ax2.set_xticklabels(x_tick_lab)
        ax2.set_xticks(x_tick)
        ax2.set_aspect(aspect)
        ax2.set_xlabel('')
        ax2.set_ylabel('')
        
        #cax = ax1.inset_axes([10, -6, 50, 1.6],transform=ax1.transData)
        #cbar = fig.colorbar(pcm,cax=cax,orientation='horizontal',ticks=ticks)
        #cbar.ax.set_xticklabels(labs)
        plt.savefig(dateKTLX+'_'+str(Hrs_time)+'_DHCA_Comp',bbox_inches = 'tight')
        plt.close()
        #sys.exit()
        ### Ref
        fig = plt.figure(figsize = (16,10))
        ax1 = plt.subplot(121)
        pcm = ax1.pcolormesh(distances,height_extract,extracted_ref,cmap='pyart_NWSRef',vmin=-10,vmax=70)
        ax1.set_ylim(cross_ylim)
        #ax1.set_xlim(cross_xlim)
        ax1.set_xticks(x_tick_ktlx)
        ax1.set_xticklabels(x_tick_lab_ktlx)
        plt.title(radartest.time['units'][14:] + '\n KTLX Reflectivity Cross Section')
        ax1.set_aspect(aspect)
        ax1.set_ylabel('')
        ax1.set_xlabel('')
        
        ax2 = plt.subplot(122)
        if Horusradar.azimuth['data'][0] > 271:
            display2.plot_rhi('Reflectivity',vmin=-10,vmax=70,cmap='pyart_NWSRef',reverse_xaxis=True,colorbar_flag=False)
        else:
            display2.plot_rhi('Reflectivity',vmin=-10,vmax=70,cmap='pyart_NWSRef',reverse_xaxis=False,colorbar_flag=False)
        ax2.set_ylim(cross_ylim)
        ax2.set_xlim([-40,-13])
        ax2.set_xticks(x_tick)
        ax2.set_xticklabels(x_tick_lab)
        #ax2.set_xticks(x_tick)
        ax2.set_aspect(aspect)
        ax2.set_xlabel('')
        ax2.set_ylabel('')
        
        cax = ax1.inset_axes([10, -6, 50, 1.6],transform=ax1.transData)
        fig.colorbar(pcm,cax=cax,orientation='horizontal')
        plt.savefig(dateKTLX+'_'+str(Hrs_time)+'_Ref_Comp',bbox_inches = 'tight')
        plt.close()
        ### ZDR
        fig = plt.figure(figsize = (16,10))
        ax1 = plt.subplot(121)
        pcm = ax1.pcolormesh(distances,height_extract,extracted_zdr,cmap='pyart_HomeyerRainbow',vmin=-1,vmax=7)
        ax1.set_ylim(cross_ylim)
        ax1.set_xlim(cross_xlim)
        ax1.set_xticks(x_tick_ktlx)
        ax1.set_xticklabels(x_tick_lab_ktlx)
        plt.title(radartest.time['units'][14:] + '\n KTLX ZDR Cross Section')
        ax1.set_aspect(aspect)
        ax1.set_ylabel('')
        ax1.set_xlabel('')
        
        ax2 = plt.subplot(122)
        if Horusradar.azimuth['data'][0] > 271:
            display2.plot_rhi('ZDR',vmin=-1,vmax=7,cmap='pyart_HomeyerRainbow',reverse_xaxis=True,colorbar_flag=False)
        else:
            display2.plot_rhi('ZDR',vmin=-1,vmax=7,cmap='pyart_HomeyerRainbow',reverse_xaxis=False,colorbar_flag=False)
        ax2.set_ylim(cross_ylim)
        ax2.set_xlim(cross_xlim)
        ax2.set_xticks(x_tick)
        ax2.set_xticklabels(x_tick_lab)
        ax2.set_xticks(x_tick)
        ax2.set_aspect(aspect)
        ax2.set_xlabel('')
        ax2.set_ylabel('')
        
        cax = ax1.inset_axes([10, -6, 50, 1.6],transform=ax1.transData)
        fig.colorbar(pcm,cax=cax,orientation='horizontal')
        plt.savefig(dateKTLX+'_'+str(Hrs_time)+'_Zdr_Comp',bbox_inches = 'tight')
        plt.close()
        
        ### RhoHV
        fig = plt.figure(figsize = (16,10))
        ax1 = plt.subplot(121)
        #ax1.pcolormesh(arrayy[:,1],height_extract,extracted_rho,cmap=cc_cmap,vmin=0.5,vmax=1)
        pcm = ax1.pcolormesh(distances,height_extract,extracted_rho,cmap=cc_cmap,vmin=0.5,vmax=1)
        ax1.set_ylim(cross_ylim)
        ax1.set_xlim(cross_xlim)
        ax1.set_xticks(x_tick_ktlx)
        ax1.set_xticklabels(x_tick_lab_ktlx)
        plt.title(radartest.time['units'][14:] + '\n KTLX RhoHV Cross Section')
        ax1.set_aspect(aspect)
        ax1.set_ylabel('')
        ax1.set_xlabel('')
        
        ax2 = plt.subplot(122)
        if Horusradar.azimuth['data'][0] > 271:
            display2.plot_rhi('RhoHV',vmin=0.5,vmax=1,cmap=cc_cmap,reverse_xaxis=True,colorbar_flag=False)
        else:
            display2.plot_rhi('RhoHV',vmin=0.5,vmax=1,cmap=cc_cmap,reverse_xaxis=False,colorbar_flag=False)
        ax2.set_ylim(cross_ylim)
        ax2.set_xlim([-40,-13])
        ax2.set_xticks(x_tick)
        ax2.set_xticklabels(x_tick_lab)
        ax2.set_xticks(x_tick)
        ax2.set_aspect(aspect)
        ax2.set_xlabel('')
        ax2.set_ylabel('')
        
        #cax = ax1.inset_axes([0.9, 0.01, 0.5, 0.06],transform=ax1.transData)
        cax = ax1.inset_axes([10, -6, 50, 1.6],transform=ax1.transData)
        fig.colorbar(pcm,cax=cax,orientation='horizontal')
        
        plt.savefig(dateKTLX+'_'+str(Hrs_time)+'_Rho_Comp',bbox_inches = 'tight')
        plt.close()
        
        ### Velocity
        fig = plt.figure(figsize = (16,10))
        ax1 = plt.subplot(121)
        pcm = ax1.pcolormesh(distances,height_extract,extracted_vel,cmap='pyart_BuDRd18',vmin=-30,vmax=30)
        ax1.set_ylim(cross_ylim)
        #ax1.set_xlim(cross_xlim)
        ax1.set_xticks(x_tick_ktlx)
        ax1.set_xticklabels(x_tick_lab_ktlx)
        plt.title(radartest.time['units'][14:] + '\n KTLX Velocity Cross Section')
        ax1.set_aspect(aspect)
        ax1.set_ylabel('Distance Above Radar (km)')
        
        ax2 = plt.subplot(122)
        if Horusradar.azimuth['data'][0] > 271:
            display2.plot_rhi('velocity',vmin=-30,vmax=30,cmap='pyart_BuDRd18',reverse_xaxis=True,colorbar_flag=False)
        else:
            display2.plot_rhi('velocity',vmin=-30,vmax=30,cmap='pyart_BuDRd18',reverse_xaxis=False,colorbar_flag=False)
        ax2.set_ylim(cross_ylim)
        #ax2.set_xlim(cross_xlim)
        ax2.set_xticks(x_tick)
        ax2.set_xticklabels(x_tick_lab)
        ax2.set_xticks(x_tick)
        ax2.set_aspect(aspect)
        ax2.set_xlabel('')
        ax2.set_ylabel('')
        cax = ax1.inset_axes([10, -6, 50, 1.6],transform=ax1.transData)
        fig.colorbar(pcm,cax=cax,orientation='horizontal')
        
        plt.savefig(dateKTLX+'_'+str(Hrs_time)+'_Vel_Comp',bbox_inches = 'tight')
        plt.close()
        #sys.exit()
        
        refl = Horusradar.fields['Reflectivity']['data']
        #fig = plt.figure()
        #plt.scatter(extracted_ref[0,:],refl[2,:])
        #plt.ylim(10,70)
        #plt.xlim(10,70)
        #plt.show()
        
        refl[refl<-1000] = np.nan
        fig = plt.figure()
        plt.hist(refl.flatten(),bins=20)
        plt.xlim(-10,70)
        plt.title('Horus Reflectivity \n ' + str(Hrs_time))
        plt.savefig(dateKTLX+'_'+str(Hrs_time)+'_HorusRefHist',bbox_inches = 'tight')
        plt.close()
        
        extracted_ref[extracted_ref<-1000] = np.nan
        fig = plt.figure()
        plt.hist(extracted_ref.flatten(),bins=20)
        plt.xlim(-10,70)
        plt.title('KTLX Reflectivity \n ' + str(Hrs_time))
        plt.savefig(dateKTLX+'_'+str(Hrs_time)+'_KTLXRefHist',bbox_inches = 'tight')
        plt.close()
        sys.exit()
        if t == len(l_Hrus):
            clicking = 'False'
#######################################################################
# ALL STUFF BELOW IS NOT NEEDED IN THE PRESENT BUT COULD BE EVENTUALLY 
#######################################################################
# Now we want to pull data along that line we've made from KTLX
# So we need to work with each elevation angle bc the size of the arrays may (are) be different
'''for h in range(0,len(lvls)):
    elevation = levels[h]
    radar = KTLXrads2[h]
    if h == 0:
        KTLXlats2,KTLXlons2,KTLXalts = radar.get_gate_lat_lon_alt(0)
        KTLXlats2= KTLXlats2.flatten()
        KTLXlons2 = KTLXlons2.flatten()
        KTLXloc2 = np.vstack([KTLXlats2,KTLXlons2]).T

        ind2 = p.contains_points(KTLXloc2, radius=0.004)
        I2 = np.shape(KTLXlats2.flatten())
        array2 = np.zeros((I2))

        array2 = updateArray(array2, ind2)
        arrayy2=KTLXloc2[ind2]
        
        Ref = radar.fields['Reflectivity']['data']
        Ref = Ref.flatten()
        Ref_stack = Ref[ind2]
        #sys.exit()

    else:
        #radar.sweep_number = 0
        KTLXlats3 = radar.gate_latitude['data']
        KTLXlons3 = radar.gate_longitude['data']
        KTLXlats3= KTLXlats3.flatten()
        KTLXlons3 = KTLXlons3.flatten()
        KTLXloc3 = np.vstack([KTLXlats3,KTLXlons3]).T

        ind3 = p.contains_points(KTLXloc3, radius=0.004)
        I3 = np.shape(KTLXlats3.flatten())
        array3 = np.zeros((I3))

        array3 = updateArray(array3, ind3)
        arrayy3=KTLXloc3[ind3]
        
        Ref = radar.fields['Reflectivity']['data']
        Ref = Ref.flatten()
        Ref_click2 = Ref[ind3]
        Ref_stack = np.vstack([Ref_stack,Ref_click2])'''
        
### FIGURE 5 - dont need this one
'''
fig = plt.figure(figsize=(18, 6))
ax1 = plt.subplot(121)
display3.plot_cross_section(
    "DHCA",
    start,
    end,
    x_axis="lon",
    cmap=hca_cmap,
    vmin=0.5,
    vmax=14.5,
    title='KTLX HCA Cross Section'
)
ax1.set_ylim(0,14)

ax2 = plt.subplot(122)
display2.plot_rhi('HCA',vmin=0.5,vmax=14.5,cmap=hca_cmap,reverse_xaxis=False)
ax2.set_ylim(0,14)'''
#ax2.set_xticklabels(x_tick_lab)

#ax2.set_xticklabels(lon_plot)

##xticks = ax2.get_xticks()
#x,y,z= Horusradar.get_gate_x_y_z(0)
#figure = plt.figure()
#ax1 = plt.subplot(111)
#display=pyart.graph.GridMapDisplay(grid)
#display.plot_cross_section('Reflectivity',start,end,x_axis='lat',cmap='pyart_NWSRef',vmin=-20,vmax=70)