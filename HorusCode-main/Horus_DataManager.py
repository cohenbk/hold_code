#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 22 13:51:34 2023

@author: lshedd123
"""
import os,sys,re
import numpy as np

firstmove = 2
var = 'CSPW'

base_path = '/Volumes/MiniHorus/Data/Horus_20230511_correct2'
num_done_path = '/Volumes/MiniHorus/Data/Horus_20230511_correct2/DHCA360'
hold_folder = '/Volumes/MiniHorus/Data/Horus_20230511_hold'
os.chdir(hold_folder)
#variables = os.listdir() # use this if firstmove is not 1
variables = ['AliasedVelocity']
#os.chdir(num_done_path+'/00.50')
#num_done = len(np.sort([f for f in os.listdir('.') if re.search('.netcdf', f)]))
#files_done = np.sort([f for f in os.listdir('.') if re.search('.netcdf', f)])
#files_path = '/Volumes/MiniHorus/Data/Horus_20230511/output/20230511_AlreadyProcessed

folders = ['00.50','01.00','01.50','02.00','02.50','03.00','03.50','04.00','04.50','05.00','05.50','06.00',
                   '06.50','07.00','07.50','08.00','08.50','09.00','09.50','10.00','10.50','11.00','11.50','12.00',
                   '12.50','13.00','13.50','14.00', '14.50','15.00','15.50','16.00','16.50','17.00','17.50','18.00',
                   '18.50','19.00','19.50','20.00','20.50','21.00','21.50','22.00','22.50','23.00','23.50','24.00',
                   '24.50','25.00','25.50','26.00','26.50','27.00','27.50','28.00','28.50','29.00','29.50','30.00',
                   '30.50']

#folders=['00.00','01.33','02.67','04.00','05.33','06.67','08.00','09.33','10.67',
                 #'12.00','13.33','14.67','16.00','17.33','18.67']
os.chdir(base_path)

if firstmove == 1:
    for i in range(0,len(folders)):
        print(folders[i])
        os.chdir(num_done_path+'/'+folders[i])
        num_done = len(np.sort([f for f in os.listdir('.') if re.search('.netcdf', f)]))
        os.chdir(num_done_path+'/00.50')
        files_done = np.sort([f for f in os.listdir('.') if re.search('.netcdf', f)])
        
        #Want to Make sure the Hold Folder is There First
        path = '/Volumes/MiniHorus/Data/Horus_20230511_hold/'+var+'/'+folders[i]
        isExist = os.path.exists(path)
        if not isExist:
            os.makedirs(path)
        
        os.chdir(base_path)
        #for j in range(2,num_done):
        for j in range(0,2):
            #move_file = files_done[j]
            os.chdir(var)
            folder = folders[i]
            os.chdir(folder)
            #sys.exit()
        #os.chdir(folder)
        
            #os.system('mv '+move_file+ ' '+hold_folder+'/'+var+'/'+folder)
            #For output variables use this: 
            os.system('mv ' +base_path +'/' + var +'/'+ folder+'/* ' +hold_folder+'/'+var+'/'+folder)
            #sys.exit()
            os.chdir(base_path)
else:
    #for k in range(4,len(variables)):
    for k in range(0,len(variables)):
        var = variables[k]
        print(var)
        hold_folder2 = '/Volumes/MiniHorus/Data/Horus_20230511_hold/'+var
        for i in range(0,len(folders)):
            os.chdir(hold_folder2)
            #if k == 3:
                #sys.exit()
        #os.chdir(var)
            folder = folders[i]
            #sys.exit()
            #os.chdir(folder)
        
            os.system('mv '+folder+'/* /Volumes/MiniHorus/Data/Horus_20230511_correct2/'+var+'/'+folder)
        
            os.chdir(base_path)
    #sys.exit()