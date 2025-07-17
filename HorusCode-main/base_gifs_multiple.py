#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 23 12:38:46 2025

@author: brandoncohen
"""

import glob
import os
from PIL import Image

# Variables and levels
variables = ['Reflectivity', 'Velocity', 'ZDR', 'RhoHV', 'PhiDP']
lvls1 = ['0075', '0225', '0375', '0525', '0675', '0825', '0975', '1125', '1275', '1425',
         '1575', '1725', '1875', '2025', '2175', '2325', '2475', '2625', '2775', '2925']

# Input directory where sorted files are located
input_dir = "/Volumes/BRANDON/27April_Figs_27March_2"

# Output directory for the GIFs
output_dir = os.path.join(input_dir, "GIFs")
os.makedirs(output_dir, exist_ok=True)

# Loop through variables and levels
for var in variables:
    for lvl in lvls1:
        # Collect images for this variable and level
        images = []
        figures = sorted(glob.glob(os.path.join(input_dir, var, lvl, "*.png")))
        
        if not figures:
            print(f"No images found for {var} at level {lvl}.")
            continue
        
        # Open each image and add to the list
        for fig in figures:
            try:
                images.append(Image.open(fig))
            except Exception as e:
                print(f"Error opening {fig}: {e}")
                continue
        
        if not images:
            print(f"No valid images to create GIF for {var} at level {lvl}.")
            continue
        
        # Output GIF path
        gif_name = f"28April_{var}_{lvl}_Scan13.gif"
        gif_path = os.path.join(output_dir, gif_name)
        
        # Save the GIF
        try:
            images[0].save(
                gif_path, save_all=True, append_images=images[1:], duration=200, loop=2
            )
            print(f"Successfully created GIF: {gif_path}")
        except Exception as e:
            print(f"Error creating GIF for {var} at level {lvl}: {e}")