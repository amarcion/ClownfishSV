import sys
import os
from Bio import SeqIO

wd = sys.argv[1]
annotation_prefix = sys.argv[2] #braker.ProtOnly
output_prefix = sys.argv[3] #braker.ProtOnly.Filtered
KeepTranscript_file = sys.argv[4] #ProtOnly.GenesToKeep.txt


def rename_chrom(name_chrom):
	name_chrom_list = name_chrom.split("_")
	ID = name_chrom_list[0]
	return(ID)

def filter_and_renameChr_gft(input_file_name, output_file_name, genestokeep):
	found_genes = []

	with open(output_file_name, "w") as output_file:
		with open(input_file_name, "r") as input_file:
			
			for line in input_file:
				line = line.rstrip()
				line_list = line.split("\t")
				#NewChromName = rename_chrom(line_list[0])
				NewChromName = line_list[0]
				line_list[0] = NewChromName

				if len(line_list[-1].split(" ")) == 1:
					if "." in line_list[-1]:
						index = line_list[-1].find(".")
						gene_id = line_list[-1][:index]
					else:
						gene_id = line_list[-1]
				else:
					gene_id = line_list[-1].split(" ")[-1].replace('"', "").replace(";", "")

				if gene_id in genestokeep:
					print("\t".join(line_list), file=output_file)
					found_genes.append(gene_id)

	#found_genes = set(found_genes)
	#found_genes = list(found_genes)
	#print(len(found_genes), len(genes))


def filter_fasta_seq(input_file_name, output_file_name, ListToKeep, entity="Transc"):
	
	with open(output_file_name, "w") as output_file:
		for record in SeqIO.parse(input_file_name, "fasta"):
			if entity == "Transc":
				toKeep = record.id
			else:
				index = record.id.find(".")
				toKeep = record.id[:index]

			if toKeep in ListToKeep:
				print(">"+record.id, file=output_file)
				print(record.seq, file=output_file)

	print("Finished ", output_file_name)


###########################

os.chdir(wd)

###
#Genes and transcripts to keep
###

genes = []
transcripts = []
input_file = open(KeepTranscript_file, "r")
for line in input_file:
	line = line.rstrip()
	index = line.find(".")
	gene = line[:index]
	genes.append(gene)
	transcripts.append(line)
input_file.close()

#########
# Filter the gtf to remove filtered genes
#########

filter_and_renameChr_gft(annotation_prefix+".gtf", output_prefix + ".gft",  genes)

#########
# Filter the amino_acid file to remove filtered genes
#########

filter_fasta_seq(annotation_prefix+".LongestIso.aa", output_prefix+".LongestIso.aa", transcripts)
filter_fasta_seq(annotation_prefix+".aa", output_prefix+".aa", genes, "Genes")

#########
# Filter the coding sequences files to remove filtered genes
#########

filter_fasta_seq(annotation_prefix+".LongestIso.codingseq", output_prefix+".LongestIso.codingseq", transcripts)
filter_fasta_seq(annotation_prefix+".codingseq", output_prefix+".codingseq", genes, "Genes")





