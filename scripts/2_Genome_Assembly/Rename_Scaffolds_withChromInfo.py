import sys
from Bio import SeqIO
import os
import argparse

def file_exists(file_path):
    """Check if a file exists."""
    if not os.path.isfile(file_path):
        raise argparse.ArgumentTypeError(f"File '{file_path}' does not exist.")
    return file_path

def ensure_path_for_file(file_path):
    """Ensure that the directory for a given file path exists."""
    directory = os.path.dirname(file_path)
    if directory and not os.path.exists(directory):
        os.makedirs(directory)
        print(f"Directory '{directory}' created.")
    else:
        print(f"Directory '{directory}' already exists or not needed.")

def get_scaffold_information(fasta_file):
	Seq_Length = []
	Dictionary_Seq = {}
	for record in SeqIO.parse(fasta_file, "fasta"):
		Seq_Length.append([len(record.seq), record.id])
		Dictionary_Seq[record.id] = record.seq
	return(Seq_Length,Dictionary_Seq)

def get_chrom_information(chr_table):
	Dictionary_chr = {}
	with open(chr_table, "r") as infile:
		for line in infile:
			line = line.rstrip()
			if line.startswith("Scaffold"):
				continue
			line_list = line.rsplit("\t")
			chrom_name = line_list[1][3:]
			if len(chrom_name) == 1:
				chrom_name = "0" + chrom_name
			Dictionary_chr[line_list[0]] = chrom_name
	return(Dictionary_chr)



parser = argparse.ArgumentParser(description="Rename scaffolds based on the potential chromosomes they belong to and the length")
parser.add_argument("-f", "--fasta_input", type=file_exists, required=True, help="Reference file in fasta format")
parser.add_argument("-o", "--out_prefix", type=str, default="Out", help="Prefix for the output files. Default: Out")
parser.add_argument("-c", "--chrom", type=file_exists, required=True, help="File name of the table containing Scaffold name \
	and chromosome names.")

#Verify the argument
args = parser.parse_args()

# create directories for output if not existing
ensure_path_for_file(args.out_prefix)


# Get info from fasta
scafflen, dictSeq = get_scaffold_information(args.fasta_input)
# Get info for chrom
chromInfo = get_chrom_information(args.chrom)

# Sort scaffolds based on length

scafflen.sort(reverse=True) # Sort from largest to smallest scaff

Dict_used_chrom = {}
increm_not_chrom = 25

Equivalence_list = []

with open(args.out_prefix+".fasta", "w") as outfile:
	for elem in scafflen:
		print(elem)
		record = elem[1]

		# check if in chromInfo
		if record in chromInfo.keys():
			chrom = chromInfo[record]
			# if yes, check if chrom already used
			if chrom not in Dict_used_chrom.keys():
				Dict_used_chrom[chrom] = 1
				new_name = "scf"+chrom+".1"
			else:
				Dict_used_chrom[chrom] += 1
				new_name = "scf"+chrom+"."+str(Dict_used_chrom[chrom])
		else:
			new_name = "scf"+str(increm_not_chrom)+".1"
			increm_not_chrom += 1

		Equivalence_list.append([record, new_name])
		# And print in the file
		print(">"+new_name, file=outfile)
		print(dictSeq[record], file=outfile)

# And then print the equivalence file

with open(args.out_prefix+".Equivalence.txt", "w") as outfile:
	for elem in Equivalence_list:
		print(elem[0]+"\t"+elem[1], file=outfile)







