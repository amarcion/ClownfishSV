from Bio import SeqIO
import sys
import os


input_file_name = sys.argv[1]
output_directory = sys.argv[2]
Species_ID = sys.argv[3]

# Get information on chromosomes
path_to_out_dir= output_directory+"/"

Dictionary_Chrom = {}
List_no_scaffold_needed = []

# Get info on scaffolding

with open(Species_ID+".Scaffolds_to_Chromosomes.ForScaffolding.txt", "r") as input_file:
	for line in input_file:
		line_list = line.rstrip().split("\t")
		chrom = line_list[0]

		if len(line_list) == 1:
			continue

		scaffs = line_list[1].split(",")
		if len(scaffs) == 1:
			List_no_scaffold_needed.append(scaffs[0])
		else:
			isExist = os.path.exists(path_to_out_dir+chrom)
			if not isExist:
				os.makedirs(path_to_out_dir+chrom)
			for scaff in scaffs:
				Dictionary_Chrom[scaff] = chrom


# Split and save the scaffolds in corresponding chrom folder and save chrom not to be scaffolded

with open(output_directory+"/"+Species_ID+".ScaffoldsNotToScaffold.fa", "w") as out_notScaff: 

	tot_scaff = 0
	total_no_Scaff = 0

	for record in SeqIO.parse(input_file_name, "fasta"):
		if record.id in List_no_scaffold_needed:
			print(">"+record.id, file=out_notScaff)
			print(record.seq, file=out_notScaff)
			total_no_Scaff += 1
		else:
			chromosome = Dictionary_Chrom[record.id]
			prefix_to_output=path_to_out_dir+chromosome+"/"+Species_ID+"."
			out_file_name = prefix_to_output+record.id+".fa"
			out_file = open(out_file_name, "w")
			print(">"+record.id, file=out_file)
			print(record.seq, file=out_file)
			out_file.close()
			tot_scaff += 1


print(f"FINISHED with file {input_file_name}")
print(f"Nb of scaffold to scaffold: {tot_scaff}")
print(f"Nb of scaffold not to scaffold: {total_no_Scaff}")
print(f"Total number of scaffolds: {tot_scaff+total_no_Scaff}")






