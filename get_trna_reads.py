
# enter script to see use case


import sys, os, re, time
os.environ['OPENBLAS_NUM_THREADS'] = '1'
import pandas as pd

start_time = time.time()

if len(sys.argv) < 1:
	print("\nUsage: This script takes in two arguments, a sam file nad a FastQ file.\n"
		"The program creates a unique list from the SAM file and pulls out all\n"
		"reads of the same named areads from the FastQ file to a new file.\n"
		"The purpose of this is to map these reads to a new genome to identify.\n"
		"potential chimeric reads. Input order is fastq sam outfile_name min_match_length min_total_read_length")
	sys.exit()


# fastq_filepath = sys.argv[1]
sam_filepath = sys.argv[1]
out_filepath = sys.argv[2]
min_match_length = int(sys.argv[3])
min_softclip_length = int(sys.argv[4])

# Manual filepaths outdated

fastq_filepath = "/project/zhipengl_72/rvandamm/rvandamm/clash_project/raw_data/fastq_raw_files/trimming/adapt_remov_no_trim/SRR959751_fastq_remov_adapt.fastq"
# sam_filepath = "/project/zhipengl_72/rvandamm/rvandamm/clash_project/raw_data/fastq_raw_files/trimming/mapping_4/SRR959751_trna_mapping/trna_mapping_data/SRR959751_trna_map_dir_Aligned.sortedByCoord.out.sam"
# out_filepath = "/project/zhipengl_72/rvandamm/rvandamm/clash_project/raw_data/fastq_raw_files/trimming/trna_mapped_fastq"


# Function takes in a CIGAR string and returns the largest match from the list of matches.
def match_from_CIGAR(CIGAR):
	matches_string = re.findall(r'(\d+)M', CIGAR)
	softclip_string = re.findall(r'(\d+)S', CIGAR)
	matches_int = [int(i) for i in matches_string]
	softclip_int = [int(i) for i in softclip_string]
	largest_match = max(matches_int)
	if not softclip_int: softclip_int.append(0)
#	print(matches_int)
#	print(softclip_int)
	largest_softclip = max(softclip_int)
	return largest_match, largest_softclip

# Block of code for troubleshooting and testing
'''
example_string = '5S7M1I2M300N18M3S'
print(match_from_CIGAR(example_string))
'''


good_reads = []
with open(sam_filepath, 'r') as in_sam_file:
	for line in in_sam_file:
		if line.startswith('@'): continue
		read = line.split()
		max_match, max_softclip = match_from_CIGAR(read[5])
		if max_match >= min_match_length and max_softclip >= min_softclip_length:
			good_reads.append(read[0])

		
unique_good_reads = set(good_reads)
# print(len(good_reads))
# print(len(unique_good_reads))

checkpoint = time.time()
sam_finished = checkpoint - start_time
print('read names extracted and filtered from SAM file', round(sam_finished/60, 2), 'minutes') 


fastq_names = []
with open(fastq_filepath, 'r') as in_fastq_file:
	for line in in_fastq_file:
		if line.startswith('@SRR'):
			fastq_names.append(line.split()[0][1:])

checkpoint = time.time()
fastq_name_time = checkpoint - start_time
print('names extracted from fastq file', fastq_name_time/60, 'minutes')		

# print(len(fastq_names))
# print(fastq_names[:10])

fastq_sam_intersection = set(fastq_names).intersection(unique_good_reads)


with open(fastq_filepath, 'r') as in_fastq_file, open(out_filepath, 'w') as out_file:
	lines = in_fastq_file.readlines()
	for i in range(len(lines)):
		if lines[i].startswith('@SRR'):
			if lines[i].split()[0][1:] in fastq_sam_intersection:
				out_file.write(lines[i])
				out_file.write(lines[i + 1])
				out_file.write(lines[i + 2])
				out_file.write(lines[i + 3])


checkpoint = time.time()
fastq_writeout = checkpoint - start_time
print('new fastq file has been written', fastq_writeout/60, 'minutes') 


# Test to see if script is working
with open(out_filepath, 'r') as infile:
	for line in infile:
		if line.startswith('@SRR'):
			if line.split()[0][1:] not in fastq_sam_intersection:
				print('TEST FAILED. Read present in file that is not from first round of mapping.')
				print(line.split()[0])
				print(line.split()[0][1:]) 
				sys.exit()

print('Test has passed. All reads in new fastq file present in the first round of mapping.')








