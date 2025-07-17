#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jan 28 14:28:57 2024

@author: brandoncohen
"""
import glob 
from PIL import Image

var = 'Reflectivity'
RHI = 0
az = '221'
lvls1 = ['0075','0225','0375','0525','0675','0825','0975','1125','1275','1425','1575',
                 '1725','1875','2025','2175','2325','2475','2625','2775','2925']
lvls2 = ['00.75','02.25','03.75','05.25','06.75','08.25','09.75','11.25','12.75','14.25','15.75',
                 '17.25','18.75','20.25','21.75','23.25','24.75','26.25','27.75','29.25']

if RHI == 0:
    for i in range(0,len(lvls1)):
        images = []
        figures = sorted(glob.glob("/Volumes/BRANDON/27April_Figs/"+var+"/"+lvls1[i]+"/*.png"))
        #figures = sorted(glob.glob("/Volumes/GROGU/Horus_Figures_From_Laura/Figures/20240428/"+var+"/"+lvls2[i]+"/*.png"))
        #figures = sorted(glob.glob(f"/Volumes/GROGU/Horus_Figures_From_Laura/Figures/20240428/*_{lvls1[i]}_{var}.png"))
        # figures = sorted(glob.glob("/Volumes/GROGU/Horus_Figures_From_Laura/Figures/20240428/VelRef_RHI_Horus_Az232_20240428_*.png"))
        #figures = sorted(glob.glob("/Volumes/GROGU/Horus_Figures_From_Laura/Figures/20240428/Vel_RHI_PPI_Horus_Az232_20240428_*.png"))
        #figures = sorted(glob.glob("/Volumes/GROGU/Horus_Figures_From_Laura/Figures/20240428/Vel_PPI_Horus_Az232_20240428_*.png"))
        #figures = sorted(glob.glob("/Volumes/GROGU/Horus_Figures_From_Laura/Figures/20240428/ZDR_RHI_Horus_Az232_20240428_*.png"))
    
        for fig in figures:
            images.append(Image.open(fig))
        
        #images[0].save("/Volumes/GROGU/Horus_Figures_From_Laura/GIFs/28April_VelRef_RHI_Scan13_Az232.gif", save_all=True, append_images=images[1:], duration=200,loop=2)
        #images[0].save("/Volumes/GROGU/Horus_Figures_From_Laura/GIFs/28April_"+var+"_VelRef_RHI_"+lvls1[i]+"_Scan13.gif", save_all=True, append_images=images[1:], duration=200,loop=2)
        images[0].save("/Volumes/BRANDON/27April_Figs/GIFs/28April_"+var+"_"+lvls1[i]+"_Scan13.gif", save_all=True, append_images=images[1:], duration=200,loop=2)

    
if RHI == 1:
    images = []
    figures = sorted(glob.glob("/Volumes/MiniHorus/Horus Figures/20240428/Scan13/RHI/"+var+"/*.png"))
    
    for fig in figures:
        images.append(Image.open(fig))
    
    images[0].save("//Users/lshedd123/Dropbox/Horus Share Folder/Horus Animations/28April_RHI_"+var+"_"+az+"_Scan13.gif", save_all=True, append_images=images[1:], duration=100,loop=1)
    
    '''
    images = []
    figures = sorted(glob.glob("/Users/brandoncohen/Desktop/Horus_gif_holder/May11_Ref_good/*.png"))

    for fig in figures:
        images.append(Image.open(fig))
    
    images[0].save("/Users/brandoncohen/Dropbox/Horus Share Folder/Horus_gifs/May11_2/Ref_good.gif", save_all=True, append_images=images[1:], duration=100,loop=2)
    '''