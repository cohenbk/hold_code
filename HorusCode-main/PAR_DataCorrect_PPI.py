#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 13 12:23:03 2023

@author: lshedd123
"""
import os, re,sys
from netCDF4 import Dataset
import numpy as np
import netCDF4

filepath = '/Volumes/GROGU/Horus_Figures_From_Laura/output_corrected'
os.chdir(filepath)
variables = ['Reflectivity','AliasedVelocity','PhiDP','RhoHV','SpectrumWidth','Zdr']
eles = ['00.62','01.88', '03.12', '04.38', '05.62', '06.88', '08.12', '09.38', '10.62', '11.88',
                '13.12', '14.38', '15.62', '16.88', '18.12', '19.38']

os.chdir(variables[0]); os.chdir(eles[0])
l = np.sort([f for f in os.listdir('.') if re.search('.netcdf', f)])

for k in range(0,len(variables)):
    for j in range(0,len(eles)):
        os.chdir(filepath)
        os.chdir(variables[k])
        os.chdir(eles[j])
        l = np.sort([f for f in os.listdir('.') if re.search('.netcdf', f)])
        for i in range(0,len(l)):
            os.chdir(filepath)
            os.chdir(variables[k])
            os.chdir(eles[j])
            file = l[i]
            filename = '20240428-'+file[9:15]+'.netcdf'
        
            #f = Dataset(file,'r+')
            
            with netCDF4.Dataset(file,mode='r') as data:
                Ref = data.variables[variables[k]]
                if variables[k] == 'Reflectivity':
                    #print(np.nanmin(Ref['data']))
                    Ref = Ref[:] - 5
                    #print(np.nanmin(Ref['data']))
                #elif variables[k] == 5:
                #    Ref['data'] = Ref['data'] + 1
                Az = np.array(data['Azimuth'])
                Beam = np.array(data['BeamWidth'])
                AzSpace = np.array(data['AzimuthalSpacing'])
                            #GateW = np.array(data['GateWidth'])[j]*window
                GateW = np.array(data['GateWidth'])
                            #print(GateW)
                            #sys.exit()
                RadT = np.array(data['RadialTime'])
                Nyquist = np.array(data['NyquistVelocity'])
                            #Ref1 = smooth_field[j]
                Ref1 = Ref[:]
                gate_dim = Ref1.shape[1]
                #Ref2 = Ref[:]
                            #sys.exit()
                attr = data.__dict__
                attr['Elevation'] = float(eles[j])
                #sys.exit()
                
                os.chdir('/Volumes/GROGU/Horus_Figures_From_Laura/output_corrected/'+variables[k]+'/'+eles[j])
                ncfile = netCDF4.Dataset(filename,mode='w',format='NETCDF4_CLASSIC')
                       #print(ncfile)
                        
                ncfile.createDimension('Gate',gate_dim)
                ncfile.createDimension('Azimuth',len(Az))
                        
                        #attr = data.__dict__
                        #print(attr)
                        #sys.exit()
                ncfile.setncatts(attr)
                        #ncfile.setncattr('Elevation', ele[j])
                     
                Azimuth = ncfile.createVariable('Azimuth',np.float32,'Azimuth')
                Azimuth.units = 'Degrees'
                Azimuth[:] = Az
                     
                Beams = ncfile.createVariable('BeamWidth',np.float32,'Azimuth')
                Beams.units = 'Degrees'
                Beams[:] = Beam
                     
                AzSpacing = ncfile.createVariable('AzimuthalSpacing',np.float32,'Azimuth')
                AzSpacing.units = 'Degrees'
                AzSpacing[:]=AzSpace
                     
                Gate = ncfile.createVariable('GateWidth',np.float32,'Azimuth')
                Gate.units = 'Meters'
                Gate[:] = GateW
                     
                RadialTime = ncfile.createVariable('RadialTime',np.float32,'Azimuth')
                RadialTime.units = 'Milliseconds'
                RadialTime[:] = RadT
                     
                NyquistV = ncfile.createVariable('NyquistVelocity',np.float32,'Azimuth')
                NyquistV.units = 'MetersPerSecond'
                NyquistV[:] = Nyquist
                #k = 0     
                Reflect = ncfile.createVariable(variables[k],np.float32,('Azimuth','Gate'))
                if variables[k] == 0:
                            Reflect.units = 'dBZ'
                elif variables[k] == 1:
                            Reflect.units = 'm/s'
                elif variables[k] == 2:
                            Reflect.units = 'degrees'
                elif variables[k] == 3:
                            Reflect.units = 'unitless'
                elif variables[k] == 4:
                            Reflect.units = 'm/s'
                elif variables[k] == 5:
                            Reflect.units = 'dB'
                Reflect[:] = Ref1
                        
                ncfile.close()

#print(f)