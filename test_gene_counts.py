import sys, os, re, time
os.environ['OPENBLAS_NUM_THREADS'] = '1'
import pandas as pd 

start_time = time.time()

if len(sys.argv) < 4:
	print("\nUsage: NAME takes in a gap1 and sam file after mapping and\n"
		"separation into filetypes. Then the left arm and right arms\n"
		"of reads are determined and using their coordinates and the\n"
		"Location of known tRNA genes from gtrnadb, tRNA derived\n"
		"smRNA reads are determined.")
	sys.exit()

trans_sam = sys.argv[1]
prigap1_sam = sys.argv[2]
min_filter = int(sys.argv[3])

# Define functions to be used for script
def get_match_count(input_CIGAR):
	pattern = r'\d+(?=M)'
	match = re.findall(pattern, input_CIGAR)
	return match


# Import trna_database.csv which contains the chromosome and start and stop
# positions of mature human tRNAs. tRNA database can be found at
# http://gtrnadb.ucsc.edu/genomes/eukaryota/Hsapi38/Hsapi38-seq.html.
database = pd.read_csv('/project/zhipengl_72/rvandamm/rvandamm/pipeline/database_files/trna_database.csv')

# Generate a list of dictionaries of tRNA positions. The dictionary is the chromosome number
# and the values are all positions of the chromosome. If for a dictionary, the start and
# stop positions for a read are in the dictionary, the read will be considered a tsRNA.

checkpoint = time.time()
trna_database_start = checkpoint - start_time
print('start of creating tRNA database', trna_database_start/60, 'minutes')

dict_list = []
for i in range(len(database)):
	# Make a list of tRNA nucleotide positions for dictionary entry
	trna_nuc_positions = []
	for nuc_pos in range(int(database.iloc[i,2]), int(database.iloc[i,3])):
		trna_nuc_positions.append(nuc_pos)
	trna_dict = {database.iloc[i,1] : trna_nuc_positions}
	dict_list.append(trna_dict)

# Rename chr1_KI270713v1 to chr1
for chromosome in dict_list:
	if 'chr1_KI270713v1_random' in chromosome: chromosome['chr1'] = chromosome.pop('chr1_KI270713v1_random')

# Get a list of all chromosome key names from the dictionaries
chromosome_names = set()
for chromosome in dict_list: chromosome_names.update(chromosome.keys())

# Create new list of dictionaries for every unique chromosome key. A master chromosome contains
# the dictionary for every trna gene in the chromosome. i.e. master_dict['chr13_master'] will print
# 3 dictionaries where each key is 'chr13' and each set of values corresponds to each position
# that the trna gene occupies in the index.
master_dict = {}
for chromosome in chromosome_names: master_dict[chromosome + '_master'] = []
for trna in dict_list:
	for key, value in trna.items():
		master_dict[key + '_master'].append(trna)

checkpoint = time.time()
trna_database_stop = checkpoint - start_time
print('tRNA database finished', trna_database_stop/60, 'minutes')

# Import trans sam file into a list of tuples
print('beginning importing data from', trans_sam)

# Find index pairs where one or both arms of the read contain gaps, insertions, or deletions
# i.e. anything but [integer]M.
indices_for_removal = []
try:
	with open(trans_sam, 'r') as trans_infile:
		for index, line in enumerate(trans_infile, 0):
			if index % 2 == 0:
				line_elements = line.strip().split('\t')
				try: _ = int(line_elements[5][:-1])
				except ValueError: indices_for_removal.append(index), indices_for_removal.append(index + 1)
			else:
				line_elements = line.strip().split('\t')
				try: _ = int(line_elements[5][:-1])
				except ValueError: indices_for_removal.append(index), indices_for_removal.append(index - 1)
except FileNotFoundError: print(trans_sam, ' File not Found')

# Remove any potential duplicate indices in the list
indices_for_removal = list(set(indices_for_removal))

# Write a new file that removes index pair lines that were indicated to have gaps, insertions, or deletions.
file_pattern = r'^(.*?)\.sam'
match = re.match(file_pattern, trans_sam)
trans_out_name = match.group(1) + '_tempfile.sam'


with open(trans_sam, 'r') as trans_infile:
	filtered_lines = [line for index, line in enumerate(trans_infile, 0) if index not in indices_for_removal]


# Look for tempfile, if it exists delete it so that a new tempfile can be made.
if os.path.exists(trans_out_name):
	try:
		os.remove(trans_out_name)
		print('Previous trans tempfile was deleted')
	except Exception as e:
		print('An error occured while removing the file:', e)
	
with open(trans_out_name, 'w') as out_file:
	for line in filtered_lines: out_file.write(line)

checkpoint = time.time()
reformat_transfile_stop = checkpoint - start_time
print('trans file reformatted', reformat_transfile_stop/60, 'minutes')

# This code is to check that all non simple Matched arms i.e. 29M and not 23M1N15M have been stricken from
# the original input trans file.

#with open(trans_out_name, 'r') as test_file:
#	for line in test_file:
#		line_elements = line.strip().split('\t')
#		try: _ = int(line_elements[5][:-1])
#		except ValueError:
#			print('This file contains arm reads with gaps, insertions, or deletions. Exiting now')
#			sys.exit()
	
#print('The file is cleaned of all gaps, insertions, or deletions in each arms values')

# Store the relevant information from the trans_sam file as a list of tuples. Each read is reduced
# to an 8mer tuple with the following syntax. 
# [(left arm name, left arm chromosome, left arm start position, left arm stop position,
#   right arm name, right arm chromosome, right arm start position, right arm stop position)]
trans_list = []
with open(trans_out_name, 'r') as trans_infile2:
	for index, line in enumerate(trans_infile2, 0):
		line_elements = line.strip().split('\t')
		if index % 2 == 0:
			left_arm = (line_elements[0], line_elements[2], int(line_elements[3]) - 1, int(line_elements[5][:-1]) + int(line_elements[3]) - 1)
		else:
			right_arm = (line_elements[0], line_elements[2], int(line_elements[3]) - 1, int(line_elements[5][:-1]) + int(line_elements[3]) - 1)
			full_read = left_arm + right_arm
			trans_list.append(full_read)

checkpoint = time.time()
trans_file_data = checkpoint - start_time
print('trans file data obtained', trans_file_data/60, 'minutes')



# Open the prigap1_sam file and extract the information as a list of tuples. Each read is reduced
# to an 6mer tuple with the following syntax.
# [(read name, chromosome number, left arm start position, left arm stop position, 
#  right arm start position, right arm stop position)] 
prigap1_list = []
with open(prigap1_sam, 'r') as prigap_file:
	for line in prigap_file:
		line_elements = line.strip().split('\t')
		N_gaps = re.findall(r'(\D|^)(\d+)N', line_elements[5])
		gap_length = int([match[1] for match in N_gaps][0])
		if (gap_length < min_filter): continue
		# Calculate the right arm start and stop position using the left arm
		# start position, the left arm length, and the gap length
		match_lengths = get_match_count(line_elements[5])
		left_arm_length, right_arm_length = int(match_lengths[0]), int(match_lengths[1])
		gap1_read = (line_elements[0], line_elements[2], int(line_elements[3]) - 1,
		int(line_elements[3]) - 1 + left_arm_length, int(line_elements[3]) - 1 + left_arm_length + gap_length,
		int(line_elements[3]) - 1 + left_arm_length + gap_length + right_arm_length)
		prigap1_list.append(gap1_read)

checkpoint = time.time()
prigap1_file_data = checkpoint - start_time
print('prigap1 file data obtained', prigap1_file_data/60, 'minutes')

# Check the trans and prigap1 reads lists and determine if any of the arms in the chimceric reads belong to
# a tsRNA. This is done by checking if the start and stop position of a read is within the genomic
# coordinates of tRNAs indicated by the tRNA database from the Lowe lab found at 
# http://gtrnadb.ucsc.edu/genomes/eukaryota/Hsapi38/Hsapi38-seq.html.

# Check trans file list for tsrna hits. If a read is found to be tRNA derived, append the read name
# to the trans_hits list. Upon finding a hit, exits the loop to prevent similar tRNA derived sequences
# from being counted multiple times.  
trans_hits = []
trans_flag = False

for index, read in enumerate(trans_list, 0):
	if (read[1] + '_master') in master_dict:			
		for inner_dict in master_dict[read[1] + '_master']:
			for trna_coords in inner_dict.values():
				if read[2] in trna_coords and read[3] in trna_coords:
					trans_hits.append(read[0])
					trans_flag = True
					break
			if trans_flag:
				trans_flag = False
				break

	if (read[5] + '_master') in master_dict:
		for inner_dict in master_dict[read[5] + '_master']:
			for trna_coords in inner_dict.values():
				if read[6] in trna_coords and read[7] in trna_coords:
					trans_hits.append(read[0])
					trans_flag = True
					break
			if trans_flag:
				trans_flag = False
				break


checkpoint = time.time()
trans_hits_finished = checkpoint - start_time
print('trans hits obtained', trans_hits_finished/60, 'minutes')

# Check prigap1 file list for tsrna hits. If a read is found to be tRNA derived, append the read name
# to the prigap1 list. Upon finding a hit, exits the loop to prevent similar tRNA derived sequences
# from being counted multiple times.
prigap1_hits = []
prigap1_flag = False

for index, read in enumerate(prigap1_list, 0):
	if (read[1] + '_master') in master_dict:
		for inner_dict in master_dict[read[1] + '_master']:
			for trna_coords in inner_dict.values():
				if read[2] in trna_coords and read[3] in trna_coords:
					prigap1_hits.append(read[0])
					prigap1_flag = True
					break
				if read[4] in trna_coords and read[5] in trna_coords:
					prigap1_hits.append(read[0])
					prigap1_flag = True
					break
			if prigap1_flag:
				prigap1_flag = False
				break


checkpoint = time.time()
prigap1_hits_finished = checkpoint - start_time
print('prigap1 hits obtained', prigap1_hits_finished/60, 'minutes')

# Testing counts for H2AC19 to compare to literature
# Define coordinates for H2AC19 (shifted by 1 at the start and end incase
# there is a shift in genome alignment
H2AC19_coords = ('chr1', 149851061 - 1, 149851624 + 1)
RPL35A_coords = ('chr3', 197950190 - 1, 197956610 + 1)
TIMM10_coords = ('chr11', 57528464 - 1, 57530803 + 1)
 
H2AC19_range = list(range(H2AC19_coords[1], H2AC19_coords[2]))
RPL35A_range = list(range(RPL35A_coords[1], RPL35A_coords[2]))
TIMM10_range = list(range(TIMM10_coords[1], TIMM10_coords[2]))

H2AC19_prigap1_hits = []
RPL35A_prigap1_hits = []
TIMM10_prigap1_hits = []

H2AC19_trans_hits = []
RPL35A_trans_hits = []
TIMM10_trans_hits = []

H2AC19_trf_count = 0
RPL35A_trf_count = 0
TIMM10_trf_count = 0


# Check hits for prigap1 reads
for read in prigap1_list:
	if read[1] == H2AC19_coords[0]:
		if read[2] in H2AC19_range and read[3] in H2AC19_range:
			H2AC19_prigap1_hits.append(read[0])
		if read[4] in H2AC19_range and read[5] in H2AC19_range:
			H2AC19_prigap1_hits.append(read[0])
	if read[1] == RPL35A_coords[0]:
		if read[2] in RPL35A_range and read[3] in RPL35A_range:
			RPL35A_prigap1_hits.append(read[0])
		if read[4] in RPL35A_range and read[5] in RPL35A_range:
			RPL35A_prigap1_hits.append(read[0])
	if read[1] == read[2] in TIMM10_coords[0]:
		if read[2] in TIMM10_range and read[3] in TIMM10_range:
			TIMM10_prigap1_hits.append(read[0])
		if read[4] in TIMM10_range and read[5] in TIMM10_range:
			TIMM10_prigap1_hits.append(read[0])


H2AC19_prigap1_hits = list(set(H2AC19_prigap1_hits))
RPL35A_prigap1_hits = list(set(RPL35A_prigap1_hits))
TIMM10_prigap1_hits = list(set(TIMM10_prigap1_hits))

checkpoint = time.time()
test_genes_found_prigap1 = checkpoint - start_time
print('test gene lists generated from prigap1 files', test_genes_found_prigap1/60, 'minutes')

for read in trans_list:
	if read[1] == H2AC19_coords[0]:
		if read[2] in H2AC19_range:
			if read[3] in H2AC19_range:
				H2AC19_trans_hits.append(read[0])
	if read[5] == H2AC19_coords[0]:
		if read[6] in H2AC19_range:
			if read[7] in H2AC19_range:
				H2AC19_trans_hits.append(read[4])
	if read[1] == RPL35A_coords[0]:
		if read[2] in RPL35A_range:
			if read[3] in RPL35A_range:
				RPL35A_trans_hits.append(read[0])
	if read[5] in RPL35A_coords[0]:
		if read[6] in RPL35A_range:
			if read[7] in RPL35A_range:
				RPL35A_trans_hits.append(read[4])
	if read[1] == TIMM10_coords[0]:
		if read[2] in TIMM10_range:
			if read[3] in TIMM10_range:
				TIMM10_trans_hits.append(read[0])
	if read[5] == TIMM10_coords[0]:
		if read[6] in TIMM10_range:
			if read[7] in TIMM10_range:
				TIMM10_trans_hits.append(read[4])

H2AC19_trans_hits = list(set(H2AC19_trans_hits))
RPL35A_trans_hits = list(set(RPL35A_trans_hits))
TIMM10_trans_hits = list(set(TIMM10_trans_hits))


checkpoint = time.time()
test_genes_found_trans = checkpoint - start_time
print('test gene lists generated from trans files', test_genes_found_trans/60, 'minutes')

H2AC19_total_hits = H2AC19_prigap1_hits + H2AC19_trans_hits
RPL35A_total_hits = RPL35A_prigap1_hits + RPL35A_trans_hits
TIMM10_total_hits = TIMM10_prigap1_hits + TIMM10_trans_hits

total_trf_hits = prigap1_hits + trans_hits

# Iterate over the list of read names for the gene of interest found in each dataset. If the read name is present
# in the tsRNA lists then that read represents an indicated crosslink between a tsRNA and the gene of interest.
# Count how many of these occur for each gene of interest.

# total_trf_hits = H2AC19_total_hits[:10] + total_trf_hits

for read_name in H2AC19_total_hits:
	if read_name in total_trf_hits:
		H2AC19_trf_count += 1

for read_name in RPL35A_total_hits:
	if read_name in total_trf_hits:
		RPL35A_trf_count += 1

for read_name in TIMM10_total_hits:
	if read_name in total_trf_hits:
		TIMM10_trf_count += 1


checkpoint = time.time()
test_gene_tsrna_time = checkpoint - start_time
print('test gene - tsRNA interaction count finished', test_gene_tsrna_time/60, 'minutes')

# Print outputs
print('raw number and percentage of tsRNA hits from trans reads', len(trans_hits), str(round((len(trans_hits)/len(trans_list)) * 100, 2)) + '%')
print('total reads from prigap1 file are', len(prigap1_list))
print('raw number and percentage of tsRNA hits from prigap1 reads', len(prigap1_hits), str(round((len(prigap1_hits)/len(prigap1_list)) * 100, 2)) + '%')

# print statements for verification of genes observed in Kumar et. al. paper
print('Total number of raw H2AC19 reads in data', len(H2AC19_total_hits))
print('Total number of raw RPL35A reads in data', len(RPL35A_total_hits))
print('Total number of raw TIMM10 reads in data', len(TIMM10_total_hits))
print('Total number of H2AC19 - tsRNA interactions in data', H2AC19_trf_count)
print('Total number of RPL35A - tsRNA interactions in data', RPL35A_trf_count)
print('Total number of TIMM10 - tsRNA interactions in data', TIMM10_trf_count)









