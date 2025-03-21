import os
from pathlib2 import Path 
import numpy as np
import pandas as pd
from random import randint

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
