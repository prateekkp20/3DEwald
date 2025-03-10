# python3 script.py
import os
import sys
from pathlib2 import Path 
import numpy as np
import pandas as pd
from random import randint
import csv
import subprocess
from subprocess import PIPE

def ChargeFile(file_path, new_charges):
    """
    Replaces the charges in the given file with new charges.
    
    Args:
        file_path (str): Path to the file to be modified.
        new_charges (list): List of new charges to replace the existing ones.
    """
    with open(file_path, 'w') as file:
        for charge in new_charges:
            file.write(f"{charge}\n")

def replace_line_in_file(file_path, start_words, new_line):
    """
    Replaces a line in the given file that starts with a specified phrase.
    
    Args:
        file_path (str): Path to the file to be modified.
        start_words (str): The initial words of the line to be replaced.
        new_line (str): The new line to replace the matched line.
    """
    with open(file_path, 'r') as file:
        lines = file.readlines()

    modified_lines = []
    for line in lines:
        if line.startswith(start_words):
            modified_lines.append(new_line + '\n')      # Replace the matching line
        else:
            modified_lines.append(line)                 # Keep other lines unchanged

    with open(file_path, 'w') as file:
        file.writelines(modified_lines)                 # Write modified content back to the file

os.system(f"clear && make clean && make")
os.chdir("run/")

ChargeFilePath = "charge.in" 
Input_File  = "input.in"
Ewald_File = "ewald.in"

NaCl = [1,-1]
ChargeFile(ChargeFilePath, NaCl)
CSV_File = "Output_Final/Exp2/New_NACL.csv"
POSCAR_File = "POSCAR_Files/Ewald3D/NaCl/POSCAR."
Total = 101

# # Just run the Ewald3D files, with simple FFTW
# with open(CSV_File, 'wb') as file:

#     for i in range(1,Total):
#         print(f"Running File POSCAR."+str(i).zfill(3))
#         replace_line_in_file(Input_File,"Posfile = ","Posfile = "+POSCAR_File+str(i).zfill(3))
#         for j in range(0,10):
#             result = subprocess.run(["./coulomb.x"], stdout=PIPE, stderr=PIPE)
#             file.write(result.stdout)  # Write stdout in binary format
#             file.write(result.stderr)  # Write stderr in binary format

# CaCl2 = [2,-1]
# ChargeFile(ChargeFilePath, CaCl2)
# CSV_File = "POSCAR_Files/Ewald3D/CaCl2/out.csv"
# POSCAR_File = "POSCAR_Files/Ewald3D/CaCl2/POSCAR."
# Total = 101

# # Just run the Ewald3D files, with simple FFTW
# with open(CSV_File, 'wb') as file:

#     for i in range(1,Total):
#         print(f"Running File POSCAR."+str(i).zfill(3))
#         replace_line_in_file(Input_File,"Posfile = ","Posfile = "+POSCAR_File+str(i).zfill(3))
#         for j in range(0,10):
#             result = subprocess.run(["./coulomb.x"], stdout=PIPE, stderr=PIPE)
#             file.write(result.stdout)  # Write stdout in binary format
#             file.write(result.stderr)  # Write stderr in binary format

# Variation of errors and time of PPPM with changing parameters(mesh, order, kvector) with various systems
with open(CSV_File, 'wb') as file:
    for i in range(1,Total,8):

        replace_line_in_file(Input_File,"Posfile = ","Posfile = "+POSCAR_File+str(i).zfill(3))

        for gridx in range(16,128,16):
            for gridy in range(16,128,16):
                for gridz in range(16,128,16):

                    replace_line_in_file(Ewald_File,"gx = ","gx = "+str(gridx))
                    replace_line_in_file(Ewald_File,"gy = ","gy = "+str(gridy))
                    replace_line_in_file(Ewald_File,"gz = ","gz = "+str(gridz))

                    for orderx in range(2,12,2):
                        for ordery in range(2,12,2):
                            for orderz in range(2,12,2):

                                replace_line_in_file(Ewald_File,"nx = ","nx = "+str(orderx))
                                replace_line_in_file(Ewald_File,"ny = ","ny = "+str(ordery))
                                replace_line_in_file(Ewald_File,"nz = ","nz = "+str(orderz))

                                for kvecx in range(4,8,1):
                                    for kvecy in range(4,8,1):
                                        for kvecz in range(4,8,1):

                                            replace_line_in_file(Ewald_File,"kx = ","kx = "+str(kvecx))
                                            replace_line_in_file(Ewald_File,"ky = ","ky = "+str(kvecy))
                                            replace_line_in_file(Ewald_File,"kz = ","kz = "+str(kvecz))

                                            for timer in range(0,5):
                                                result = subprocess.run(["./coulomb.x"], stdout=PIPE, stderr=PIPE)
                                                file.write(result.stdout)  # Write stdout in binary format
                                                file.write(result.stderr)  # Write stderr in binary format

