import sys

input_table = sys.argv[1]				# File of cooridnates for whole genome alignments
input_name_length_scaff = sys.argv[2]	# File with scaffold lengths
theshold_coverage = float(sys.argv[3])	# Threshold coverage
theshold_identity = float(sys.argv[4])	# Threshold identity
output_file_name =  sys.argv[5]			# Output file name


# Get Scaffold Length
Scaff_length = {}

with open(input_name_length_scaff, "r") as input_file:
	for line in input_file:
		line = line.rstrip()
		line_list = line.split("\t")
		Scaff_length[line_list[0]] = line_list[1]


# Get info alignments
Align_info = {}
with open(input_table, "r") as input_file:
	for line in input_file:
		line = line.rstrip()
		if line.startswith("/") or line.startswith("[") or line.startswith("NUC") or line =="":
			continue

		line_list = line.split("\t")
		identity = line_list[6]
		Al_length = line_list[4]
		chrom = line_list[-2]
		scaff = line_list[-1]
		if scaff not in Align_info.keys():
			Align_info[scaff] = {}

		if chrom not in Align_info[scaff].keys():
			Align_info[scaff][chrom] = [float(Al_length), float(identity)]
		else:
			mean_id = (Align_info[scaff][chrom][1]+float(identity))/2
			newLegth = Align_info[scaff][chrom][0]+float(Al_length)
			Align_info[scaff][chrom] = [newLegth, mean_id]

# Check the percentage of alignments

with open(output_file_name, "w") as output_file:

	print("Scaffold\tchromosome\tpercentAlignment\tMeanIdentity", file=output_file)
	Potential_misassemblies = {}

	for scaffold in Align_info.keys():
		scaff_length = Scaff_length[scaffold]

		list_align = []
		for chrom in Align_info[scaffold]:

			Perc_length = round((Align_info[scaffold][chrom][0]/float(scaff_length))*100, 2)
			mean_id = Align_info[scaffold][chrom][1]

			if Perc_length > theshold_coverage and mean_id > theshold_identity:
				list_align.append([chrom, str(Perc_length), str(mean_id)])
				what_to_print = scaffold+"\t"+chrom+"\t"+str(Perc_length)+"\t"+str(mean_id)
				print(what_to_print, file=output_file)

		if len(list_align) > 1:
			Potential_misassemblies[scaffold] = list_align



# Print potential miassemblies
for scaffold in Potential_misassemblies:
	temp_chr = []
	temp_Align = []	
	temp_Ident = []	
	for align in Potential_misassemblies[scaffold]:
		temp_chr.append(align[0])
		temp_Align.append(align[1])
		temp_Ident.append(align[2])
	print(scaffold, ",".join(temp_chr), ",".join(temp_Align), ",".join(temp_Align), sep="\t")


