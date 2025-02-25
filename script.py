import os
from pathlib2 import Path 
import numpy as np
import pandas as pd
from random import randint

# function to edit a line in a file
def replace_line_in_file(file_path, start_words, new_line):
    with open(file_path, 'r') as file:
        lines = file.readlines()

    modified_lines = []
    for line in lines:
        if line.startswith(start_words):
            modified_lines.append(new_line + '\n')
        else:
            modified_lines.append(line)

    with open(file_path, 'w') as file:
        file.writelines(modified_lines)


