"""
Created on Tue Apr 9 2024
rdvandamme@gmail.com
"""

# READ ME
# This script is designed to take the tRNA file that contains genomic indices and
# sequence information from the Lowe lab (link below) and translate it into a database
# file that is easy to access in csv format. The resulting file contains rows with every
# mature tRNA sequence with the following format.
# trna_name, chromosome_number, start_position, stop_position, sequence

# Loew tRNA files can be downloaded from the link below
# https://gtrnadb.ucsc.edu/genomes/eukaryota/Hsapi38/Hsapi38-seq.html

# Turn off os.environ['OPENBLAS_NUM_THREADS'] = '1' for running jobs

import sys, os, re
os.environ['OPENBLAS_NUM_THREADS'] = '1'
import pandas as pd


work_dir = "/project/zhipengl_72/rvandamm/rvandamm/pipeline/database_files"
os.chdir(work_dir)

# Define functions for script
def get_trna_coords(file_line):
	pattern = r'chr[^.:]+:\s*([^ ]+)'
	chr_pattern = r'[^:]+'
	coords_pattern = r'(\d+)-(\d+)'
	match = re.search(pattern, file_line)
	chromosome_match = re.search(chr_pattern, match.group(0))
	chromosome_num = chromosome_match.group()
	coords_match = re.search(coords_pattern, match.group(1))
	start_position = int(coords_match.group(1)) - 1
	stop_position = int(coords_match.group(2))
	line_elements = file_line.split()
	trna_name = line_elements[0]
	return trna_name, chromosome_num, start_position, stop_position



with open('trna_seq_list.txt', 'r') as infile:
	trna_seq_list = infile.readlines()

# The tRNA sequence is split bnetween lines 2 and 3. Add the sequence tog-
# ether on line 2 and delete line 3.

for i in range(2, len(trna_seq_list), 3):
	trna_seq_list[i-1] = trna_seq_list[i-1].strip() + trna_seq_list[i]
del trna_seq_list[2::3]

trna_info_list = []
for i in range(len(trna_seq_list)):
	if trna_seq_list[i][0] == '>' :
		trna_name, chromosome_num, start_position, stop_position = get_trna_coords(trna_seq_list[i])
		trna_sequence = trna_seq_list[i+1].rstrip()
		trna_info_list.append((trna_name, chromosome_num, start_position, stop_position, trna_sequence))


trna_df = pd.DataFrame(trna_info_list, columns = ['trna_name', 'chromosome_number', 'start_position', 'stop_position', 'sequence'])

trna_df.to_csv('trna_database.csv', index=False)

# for line in trna_seq_list: print(line)

# Generate tRNA dataframe. Name is column 1, tRNA sequence is column 2,
# RNA distinction is column 3.


