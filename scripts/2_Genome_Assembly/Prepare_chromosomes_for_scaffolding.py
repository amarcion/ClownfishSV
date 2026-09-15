import sys

tot_scaff_file_name = sys.argv[1]
input_name = sys.argv[2]
input_file_scaff_name = sys.argv[3]
output_name = sys.argv[4]


## Get all scaffolds
with open(tot_scaff_file_name, "r") as input_file:
	tot_scaff = [scaff.strip() for scaff in input_file]


### link scaff and chrom together

with open(input_name, "r") as input_file:

	Dictionary_scaff_chrom = {}
	Dictionary_chrom_scaff= {}

	for line in input_file:
		line_list = line.rstrip().split("\t")
		if line_list[0] == "Chr":
			continue

		chrom = line_list[0] 
		scaffolds = line_list[2].split(",")

		Dictionary_chrom_scaff[chrom] = scaffolds
		
		for scaff in scaffolds:
			if scaff not in Dictionary_scaff_chrom.keys():
				Dictionary_scaff_chrom[scaff] = [line_list[0]]
			else:
				Dictionary_scaff_chrom[scaff].append(line_list[0])



# Add information on small scaffoldaligning well
with open(input_file_scaff_name, "r") as input_file_scaff:
	for line in input_file_scaff:
		line_list = line.rstrip().split()
		if line_list[0] == "Scaff":
			continue
		if float(line_list[-1]) >= 25.0:
			scaff = line_list[0]
			chrom = line_list[1]
			if scaff in Dictionary_scaff_chrom.keys():
				if chrom not in Dictionary_scaff_chrom[scaff]:
					Dictionary_scaff_chrom[scaff].append(chrom)
			else:
				Dictionary_scaff_chrom[scaff] = [chrom]


# Update the Dictionary_chrom_scaff
for scaff in Dictionary_scaff_chrom.keys():
	chroms = Dictionary_scaff_chrom[scaff]
	for chrom in chroms:
		if scaff not in Dictionary_chrom_scaff[chrom]:
			Dictionary_chrom_scaff[chrom].append(scaff)


#Get the chromosomes to scaffolds together
New_keys = {}
for chrom in Dictionary_scaff_chrom.values():
	if len(chrom) > 1:
		new_key = "_".join(chrom)
		for single_chr in chrom:
			if single_chr not in Dictionary_scaff_chrom.keys():
				New_keys[single_chr] = new_key
			else:
				New_keys[single_chr].append(new_key)


#Create the file with information on scaffolding
Dictionary_final_scaff_info = {}
for chrom in Dictionary_chrom_scaff.keys():
	variable_for_folder = ""
	if chrom not in New_keys.keys():
		variable_for_folder=chrom
	else:
		variable_for_folder=New_keys[chrom]

	if variable_for_folder not in Dictionary_final_scaff_info.keys():
		Dictionary_final_scaff_info[variable_for_folder] = Dictionary_chrom_scaff[chrom]
	else:
		Dictionary_final_scaff_info[variable_for_folder] += Dictionary_chrom_scaff[chrom]

# Get remaining scaffolds
All_considered_scaff = []
for ch in Dictionary_final_scaff_info.keys():
	scaff = list(set(Dictionary_final_scaff_info[ch]))
	All_considered_scaff += scaff

Other = []
for scaff in tot_scaff:
	if scaff not in All_considered_scaff:
		Other.append(scaff)


# Save information
with open(output_name, "w") as output_file:
	for chrom in Dictionary_final_scaff_info.keys():
		scaff = set(Dictionary_final_scaff_info[chrom])
		print(chrom, ",".join(scaff), sep="\t", file=output_file)
		print(chrom, ",".join(scaff), sep="\t")
	print("Other", ",".join(Other), sep="\t", file=output_file)
	print("Other", ",".join(Other), sep="\t")




