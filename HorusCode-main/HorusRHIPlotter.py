#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Horus RHI Plotter
Created on Wed Jan 25 12:30:26 2023

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
from datetime import datetime

from six.moves import input
from matplotlib.path import Path
import netCDF4 as nc
import sys
from datetime import datetime,timezone,timedelta

epoch = 0
def read_KOUN(idz):
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
            Ref=np.array(data['Reflectivity'])
            Ref[Ref < -1000] = np.nan
            I,J=Ref.shape
            #Az=np.array(nc['Azimuth'])

            # Fill the Radar Object with important fields
            rad = pyart.testing.make_empty_ppi_radar(J, I, 1)
            #print(J,I)
            time_dict=get_metadata('Time')
            if epoch == 1: 
                time_1 = datetime.fromtimestamp(int(np.array([data.Time])))
                time_1 = time_1 - timedelta(hours=5, minutes=0)
                time_dict['data'] = np.array([int(time_1.timestamp())])
            else:
                time_dict['data']=np.array([data.Time])
            time_dict['units']='seconds since 1970-01-01T00:00:00Z'
            time_dict['calendar']='gregorian'
            rad.time = time_dict
            rad.latitude['data'] = np.array([data.Latitude])
            rad.longitude['data'] = np.array([data.Longitude])
            rad.sweep_number['data'] = np.array([1])
            rad.azimuth['data'] = np.array(data['Azimuth'])
            azs = rad.azimuth['data'][0]
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
            ref_dict['data'] = np.array(data['Reflectivity'])
            ref_dict['units']='dBZ'
            rad.add_field('Reflectivity', ref_dict)
            
            #RhoHV
            os.chdir(Rhofolder)
            l2=np.sort([f for f in os.listdir('.') if re.search('.netcdf', f)])
            file_path=idz
            with netCDF4.Dataset(file_path,mode='r') as data:
                rho_dict = get_metadata('RhoHV')
                Rho = np.array(data['RhoHV'])
                Rho[Rho < -1000] = np.nan
                rho_dict['data'] = Rho
                rho_dict['units'] = 'unitless'
                rad.add_field('RhoHV',rho_dict)
            
            #ZDR
            os.chdir(Zdrfolder)
            l2=np.sort([f for f in os.listdir('.') if re.search('.netcdf', f)])
            file_path=idz
            with netCDF4.Dataset(file_path,mode='r') as data:
                zdr_dict = get_metadata('Zdr')
                Zdr = np.array(data['Zdr'])
                Zdr[Zdr < -1000] = np.nan
                zdr_dict['data'] = Zdr
                zdr_dict['units'] = 'dB'
                rad.add_field('ZDR',zdr_dict)
               
            #HCA
            
            os.chdir(Hcafolder)
            l2=np.sort([f for f in os.listdir('.') if re.search('.netcdf', f)])
            file_path=l2[j]
            #print(file_path)
            with netCDF4.Dataset(file_path,mode='r') as data:
                #print(data)
                    #global hca_dict
                    hca_dict = get_metadata('HCA')
                    hca_dict['data'] = np.array(data['DHCA360'][1,:].reshape((1,J)))
                    #sys.exit()
                    #print(data['DHCA360'][0,:].shape)
                    hca_dict['units'] = 'unitless'
                    rad.add_field('HCA',hca_dict)
           
            #KDP
            
            os.chdir(Kdpfolder)
            l2=np.sort([f for f in os.listdir('.') if re.search('.netcdf', f)])
            file_path=l2[j]
            with netCDF4.Dataset(file_path,mode='r') as data:
                    kdp_dict = get_metadata('DKDP')
                    kdp_dict['data'] = np.array(data['DKDP360'][1,:].reshape((1,J)))
                    kdp_dict['units'] = 'deg/km'
                    rad.add_field('DKDP',kdp_dict)
                    
            # PhiDP
            os.chdir(Phifolder)
            l2=np.sort([f for f in os.listdir('.') if re.search('.netcdf', f)])
            file_path=l2[j]
            with netCDF4.Dataset(file_path,mode='r') as data:
                    kdp_dict = get_metadata('PhiDP')
                    PhiDP = np.array(data['PhiDP'])
                    PhiDP[PhiDP < -1000] = np.nan
                    kdp_dict['data'] = PhiDP
                    kdp_dict['units'] = 'deg'
                    rad.add_field('PhiDP',kdp_dict)
                    
            #Velocity
            os.chdir(Velfolder)
            l2=np.sort([f for f in os.listdir('.') if re.search('.netcdf', f)])
            file_path=l2[j]
            with netCDF4.Dataset(file_path,mode='r') as data:
                    kdp_dict = get_metadata('Velocity')
                    Vel = np.array(data['AliasedVelocity'])
                    Vel[Vel < -1000] = np.nan
                    kdp_dict['data'] = Vel
                    kdp_dict['units'] = 'm/s'
                    rad.add_field('Velocity',kdp_dict)
            
            #RZDR
            '''
            os.chdir(rzdrfolder)
            l2=np.sort([f for f in os.listdir('.') if re.search('.netcdf', f)])
            file_path=l2[j]
            with netCDF4.Dataset(file_path,mode='r') as data:
                #print(data)
                r_dict = get_metadata('RZDR')
                r_dict['data'] = np.array(data['RZDR'])
                r_dict['units'] = 'mm/hr'
                rad.add_field('RZDR',r_dict)'''
                                                 
    return rad

#lvls1 = ['00.00','01.00','02.00','03.00','04.00','05.00','06.00','07.00','08.00','09.00','10.00'
          # ,'11.00','12.00','13.00','14.00','15.00','16.00','17.00','18.00','19.00','20.00','21.00',
         #  '22.00','23.00','24.00','25.00','26.00','27.00','28.00','29.00','30.00']
#lvls1 = ['00.00','01.33','02.67','04.00','05.33','06.67','08.00','09.33','10.67','12.00'
 #              ,'13.33','14.67','16.00','17.33','18.67']
#removed the 0.5 elevation angle for the time being because of some timing issues in needing to plot
# removed the 0.5, 1.0,1.5 to test GC impacts
lvls1 = ['00.50','01.00','01.50','02.00','02.50','03.00','03.50','04.00','04.50','05.00','05.50','06.00',
                   '06.50','07.00','07.50','08.00','08.50','09.00','09.50','10.00','10.50','11.00','11.50','12.00',
                   '12.50','13.00','13.50','14.00', '14.50','15.00','15.50','16.00','16.50','17.00','17.50','18.00',
                   '18.50','19.00','19.50','20.00','20.50','21.00','21.50','22.00','22.50','23.00','23.50','24.00',
                   '24.50','25.00','25.50','26.00','26.50','27.00','27.50','28.00','28.50','29.00','29.50','30.00',
                   '30.50']
#eles = [0.00,1.00,2.00,3.00,4.00,5.00,6.00,7.00,8.00,9.00,10.00
          # ,11.00,12.00,13.00,14.00,15.00,16.00,17.00,18.00,19.00,20.00,21.00,
          # 22.00,23.00,24.00,25.00,26.00,27.00,28.00,29.00,30.00]
#eles = [00.00,01.33,02.67,04.00,05.33,06.67,08.00,09.33,10.67,12.00
#               ,13.33,14.67,16.00,17.33,18.67]
eles = [00.50,01.00,01.50,02.00,02.50,03.00,03.50,04.00,04.50,05.00,05.50,06.00,
                   06.50,07.00,07.50,08.00,08.50,09.00,09.50,10.00,10.50,11.00,11.50,12.00,
                   12.50,13.00,13.50,14.00, 14.50,15.00,15.50,16.00,16.50,17.00,17.50,18.00,
                   18.50,19.00,19.50,20.00,20.50,21.00,21.50,22.00,22.50,23.00,23.50,24.00,
                   24.50,25.00,25.50,26.00,26.50,27.00,27.50,28.00,28.50,29.00,29.50,30.00,
                   30.50]

date = '20230808'
dhca_lvl = '01.00'

Reffolder = '/Volumes/MiniHorus/Data/Horus_20230808_Processed/output/Reflectivity/'+lvls1[0]
Rhofolder = '/Volumes/MiniHorus/Data/Horus_20230808_Processed/output/RhoHV/'+lvls1[0]
Zdrfolder = '/Volumes/MiniHorus/Data/Horus_20230808_Processed/output/Zdr/'+lvls1[0]
Velfolder = '/Volumes/MiniHorus/Data/Horus_20230808_Processed/output/AliasedVelocity/'+lvls1[0]
Phifolder = '/Volumes/MiniHorus/Data/Horus_220230808_Processed/output/PhiDP/'+lvls1[0]
#Meshfolder = '/Volumes/MiniHorus/Data/Horus_20230511/MESH/at'+dhca_lvl
#Hcafolder = '/Volumes/MiniHorus/Data/Horus_20230511_correct/DHCA/'+lvls1[0]
#Kdpfolder = '/Volumes/MiniHorus/Data/Horus_20230511_correct2/DKDP/'+lvls1[0]
#rzdrfolder = '/Volumes/MiniHorus/Data/Horus_20230511_correct2/RZDR/'+lvls1[0]

'''Reffolder = '/Users/lshedd123/Dropbox/Horus Share Folder/Horus wdssii files/Horus_20230511_Downsampled/Reflectivity/'+lvls1[0]
Rhofolder = '/Users/lshedd123/Dropbox/Horus Share Folder/Horus wdssii files/Horus_20230511_Downsampled/RhoHV/'+lvls1[0]
Zdrfolder = '/Users/lshedd123/Dropbox/Horus Share Folder/Horus wdssii files/Horus_20230511_Downsampled/Zdr/'+lvls1[0]
Velfolder = '/Users/lshedd123/Dropbox/Horus Share Folder/Horus wdssii files/Horus_20230511_Downsampled/AliasedVelocity/'+lvls1[0]
#Meshfolder = '/Volumes/MiniHorus/Data/Horus_20230511/MESH/at'+dhca_lvl
Hcafolder = '/Users/lshedd123/Dropbox/Horus Share Folder/Horus wdssii files/Horus_20230511_Downsampled/DHCA/'+lvls1[0]
Phifolder = '/Users/lshedd123/Dropbox/Horus Share Folder/Horus wdssii files/Horus_20230511_Downsampled/PhiDP/'+lvls1[0]
Kdpfolder = '/Users/lshedd123/Dropbox/Horus Share Folder/Horus wdssii files/Horus_20230511_Downsampled/DKDP/'+lvls1[0]'''

#figfolder = '/Volumes/MiniHorus/Horus/Horus Figures'
figfolder = '/Users/lshedd123/Dropbox/Horus Share Folder/Horus Figures/Horus_20230808_Processed'

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
cc_cmap = colors.LinearSegmentedColormap('testCmap', segmentdata=cdict, N=512)
hca_cmap = colors.ListedColormap(['#AAAAAA','#0000FE','#0057FF','#00A9FF','#00FFFF','#56FFAA','#ACFD56','#FFFF00','#FEAB02','#FF5700','#FF00FF','#B300B3','#DC143C','#FF00FF'])
hca_cmap = colors.ListedColormap(['limegreen','green','#DC143C','#FFFF00','dimgray','lightgray','#B300B3','white','#00FFFF','#0000FE','mediumturquoise','#000000','#FF00FF','#000000'])
ticks = np.arange(1,15)
labs = ['L/MR','HR','R/Ha','BD','AP','BI','UK','NE','DS','WS','CY','','GR','+']

os.chdir(Reffolder)
l = np.sort([f for f in os.listdir('.') if re.search('.netcdf', f)])
for j in range(0,len(l)):
    print(str(j)+'. '+l[j])
    
in_idx = input('Choose Your Start Time: ')
plot_type='Ref'
#storm_num = input('What Do You Want to Call The Storm? ')
clicking='True'
azimuth_select = 1

rads = []
for j in range(int(in_idx),len(l)):
#temporary one while the HCA is still running
#for j in range(int(in_idx),43):
    print(j)
# Loop Over All Elevation Angles at a Given Time
    for i in range(0,len(lvls1)):
        #print(i)
        Reffolder = '/Volumes/MiniHorus/Data/Horus_20230808_Processed/output/Reflectivity/'+lvls1[i]
        Rhofolder = '/Volumes/MiniHorus/Data/Horus_20230808_Processed/output/RhoHV/'+lvls1[i]
        Zdrfolder = '/Volumes/MiniHorus/Data/Horus_20230808_Processed/output/Zdr/'+lvls1[i]
        Velfolder = '/Volumes/MiniHorus/Data/Horus_20230808_Processed/output/AliasedVelocity/'+lvls1[i]
        #Meshfolder = '/Volumes/MiniHorus/Data/Horus_20230511/MESH/at'+lvls1[i]
        Hcafolder = '/Volumes/MiniHorus/Data/Horus_20230808_Processed/output/DHCA360/'+lvls1[i]
        Kdpfolder = '/Volumes/MiniHorus/Data/Horus_20230808_Processed/output/DKDP360/'+lvls1[i]
        Phifolder = '/Volumes/MiniHorus/Data/Horus_20230808_Processed/output/PhiDP/'+lvls1[i]
        
        '''Reffolder = '/Users/lshedd123/Dropbox/Horus Share Folder/Horus wdssii files/Horus_20230511_Downsampled/Reflectivity/'+lvls1[i]
        Rhofolder = '/Users/lshedd123/Dropbox/Horus Share Folder/Horus wdssii files/Horus_20230511_Downsampled/RhoHV/'+lvls1[i]
        Zdrfolder = '/Users/lshedd123/Dropbox/Horus Share Folder/Horus wdssii files/Horus_20230511_Downsampled/Zdr/'+lvls1[i]
        Velfolder = '/Users/lshedd123/Dropbox/Horus Share Folder/Horus wdssii files/Horus_20230511_Downsampled/AliasedVelocity/'+lvls1[i]
        #Meshfolder = '/Volumes/MiniHorus/Data/Horus_20230511/MESH/at'+dhca_lvl
        Phifolder = '/Users/lshedd123/Dropbox/Horus Share Folder/Horus wdssii files/Horus_20230511_Downsampled/PhiDP/'+lvls1[i]
        Hcafolder = '/Users/lshedd123/Dropbox/Horus Share Folder/Horus wdssii files/Horus_20230511_Downsampled/DHCA/'+lvls1[i]
        Kdpfolder = '/Users/lshedd123/Dropbox/Horus Share Folder/Horus wdssii files/Horus_20230511_Downsampled/DKDP/'+lvls1[i]'''
        
        os.chdir(Reffolder)
        l = np.sort([f for f in os.listdir('.') if re.search('.netcdf', f)])
        if i == 0:
            file = l[int(j)]
            
            save_time = file[9:15]
            l2 = file
            #idx_2 = np.where(l==l2)[0][0]
        else:
            l = np.sort([f for f in os.listdir('.') if re.search('.netcdf', f)])[0]
            #print(l)
            #sys.exit()
            l2 = np.sort([f for f in os.listdir('.') if re.search(save_time, f)])[0]
            #idx_2 = np.where(l==l2)[0][0]
        
        if len(l2) == 0:
            print('Skipping Iteration')
            continue
            
        rad = read_KOUN(l2)
        
        rads.append(rad)
    
    while len(rads) > 1:
        bob = rads.pop(1)
        #print('bob')
        #print(bob.nrays)
        rads[0] = pyart.util.join_radar(rads[0], bob)
                
    radar=rads[0]
    angle = radar.azimuth['data'][0]
    xlim1 = (min(radar.range['data']) - 5)/1000
    xlim2 = (max(radar.range['data']) + 5)/1000    
    
    times= radar.time['data'][0]
    ts = datetime.fromtimestamp(times,timezone.utc).strftime('%Y-%m-%d %H:%M:%S')
    
    # PLOT
    os.chdir(figfolder)
    ylim1 = 0
    ylim2 = 12
    xlim1 = 0
    xlim2 = 20
    aspects = abs((ylim1-ylim2)/(xlim2-xlim1))
    #print('IM DOWN HERE')

    xsect=pyart.util.cross_section_ppi(radar,[angle])
    display=pyart.graph.RadarDisplay(xsect)
    
    fig=plt.figure()
    display.plot('Reflectivity',vmin=-10,vmax=70,cmap='pyart_NWSRef')
    plt.title(ts+'Z \n Reflectivity')
    #im = plt.gca().collections[0]
    #plt.colorbar(im,fraction=0.046)
    #plt.show()
    #sys.exit()
    #plt.ylim([ylim1,ylim2])
    #plt.xlim([xlim1,xlim2])
    #display.set_aspect_ratio(aspect_ratio = aspects)
    plt.savefig(date+save_time+'_Ref_RHI')
    plt.close()
    #plt.xlim([xlim1,xlim2])
    #plt.ylim([0,20])
    
    fig=plt.figure()
    display.plot('ZDR',vmin=-1,vmax=7,cmap='pyart_HomeyerRainbow')
    plt.title(ts+'Z \n ZDR')
    #plt.ylim([ylim1,ylim2])
    #plt.xlim([xlim1,xlim2])
    #display.set_aspect_ratio(aspect_ratio = aspects)
    plt.savefig(date+save_time+'_Zdr_RHI')
    plt.close()
    
    fig=plt.figure()
    display.plot('RhoHV',vmin=0.5,vmax=1.02,cmap=cc_cmap)
    plt.title(ts+'Z \n RhoHV')
    #plt.ylim([ylim1,ylim2])
    #plt.xlim([xlim1,xlim2])
    #display.set_aspect_ratio(aspect_ratio = aspects)
    plt.savefig(date+save_time+'_RhoHV_RHI')
    plt.close()
    #plt.xlim([xlim1,xlim2])
    #plt.ylim([0,20])
    
    fig=plt.figure()
    display.plot('HCA',cmap=hca_cmap,vmin=0.5,vmax=14.5,ticklabs=labs,ticks=ticks)
    plt.title(ts+'Z \n HCA')
    #plt.ylim([ylim1,ylim2])
    #plt.xlim([xlim1,xlim2])
    #display.set_aspect_ratio(aspect_ratio = aspects)
    plt.savefig(date+save_time+'_HCA_RHI')
    plt.close()
    
    fig=plt.figure()
    display.plot('DKDP',cmap='pyart_RefDiff',vmin=-2,vmax=10)
    plt.title(ts+'Z \n KDP')
    #plt.ylim([ylim1,ylim2])
    #plt.xlim([xlim1,xlim2])
    #display.set_aspect_ratio(aspect_ratio = aspects)
    plt.savefig(date+save_time+'_KDP_RHI')
    plt.close()
    
    fig=plt.figure()
    display.plot('PhiDP',cmap='pyart_Wild25',vmin=-180,vmax=180)
    plt.title(ts+'Z \n PhiDP')
    #plt.ylim([ylim1,ylim2])
    #plt.xlim([xlim1,xlim2])
    #display.set_aspect_ratio(aspect_ratio = aspects)
    plt.savefig(date+save_time+'_PhiDP_RHI')
    plt.close()
    
    fig=plt.figure()
    display.plot('Velocity',cmap='pyart_BuDRd18',vmin=-30,vmax=30)
    plt.title(ts+'Z \n Velocity')
    #plt.ylim([ylim1,ylim2])
    #plt.xlim([xlim1,xlim2])
    #display.set_aspect_ratio(aspect_ratio = aspects)
    plt.savefig(date+save_time+'_Velocity_RHI')
    plt.close()
    #sys.exit()
    radar=[]
    rads=[]
    
    #sys.exit()
    