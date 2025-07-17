#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 23 12:23:45 2025

@author: brandoncohen
"""

import os
import shutil

# Input directory containing the files
input_dir = "/Volumes/BRANDON/27April_Figs_27March_2"

# Output directory where sorted folders will be created
output_dir = "/Volumes/BRANDON/27April_Figs_27March_2"

# Levels and variables
lvls1 = ['0075', '0225', '0375', '0525', '0675', '0825', '0975', '1125', '1275', '1425',
         '1575', '1725', '1875', '2025', '2175', '2325', '2475', '2625', '2775', '2925']
variables = ['Velocity', 'Reflectivity', 'ZDR', 'RhoHV', 'PhiDP']

# Create output directories
for var in variables:
    for lvl in lvls1:
        folder_path = os.path.join(output_dir, var, lvl)
        os.makedirs(folder_path, exist_ok=True)

# Process and move files
for file_name in os.listdir(input_dir):
    if not file_name.endswith(".png"):  # Skip non-PNG files
        continue
    
    # Parse the file name (e.g., "20240428_021525_0225_Velocity")
    parts = file_name.split('_')
    if len(parts) < 4:
        print(f"Skipping unexpected file name format: {file_name}")
        continue
    
    level = parts[2]  # Extract the level (e.g., "0225")
    variable = parts[3].split('.')[0]  # Extract the variable (e.g., "Velocity")
    
    # Check if level and variable are valid
    if level not in lvls1 or variable not in variables:
        print(f"Skipping file due to unmatched level/variable: {file_name}")
        continue
    
    # Determine destination folder
    dest_folder = os.path.join(output_dir, variable, level)
    dest_path = os.path.join(dest_folder, file_name)
    
    # Move the file
    shutil.move(os.path.join(input_dir, file_name), dest_path)
    print(f"Moved {file_name} to {dest_folder}")
