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

Input_File  = "input.in"
CSV_File = "POSCAR_Files/Ewald3D/NaCl/out.csv"
POSCAR_File = "POSCAR_Files/Ewald3D/NaCl/POSCAR."
Total = 8

os.system(f"clear && make clean && make")
os.chdir("run/")

# # Write the header to the CSV file in text mode
# with open(CSV_File, 'w', newline='') as csvfile:
#     csvwriter = csv.writer(csvfile)
#     csvwriter.writerow(['Ecoul', 'Time (sec)', 'Elong', 'Time(sec)', 'Elong(pppm)', 'Time(sec)', 'Total'])

# Just run the Ewald3D files, with simple FFTW
with open(CSV_File, 'wb') as file:

    for i in range(1,Total):
        replace_line_in_file(Input_File,"Posfile = ","Posfile = "+POSCAR_File+str(i).zfill(3))
        for j in range(0,5):
            result = subprocess.run(["./coulomb.x"], stdout=PIPE, stderr=PIPE)
            file.write(result.stdout)  # Write stdout in binary format
            file.write(result.stderr)  # Write stderr in binary format