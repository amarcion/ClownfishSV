import sys

input_file_name = sys.argv[1]
file_with_scaff_and_Chrom = sys.argv[2]
out_name_prefix = sys.argv[3]


# Get chromosomes and scaffolds from file
Scaff_and_chroms = {}
with open(file_with_scaff_and_Chrom, "r") as infile:
	for line in infile:
		line = line.rstrip()
		line_list = line.split("\t")
		if line_list[0] not in Scaff_and_chroms.keys():
			Scaff_and_chroms[line_list[0]] = [line_list[1].split(",")]
		else:
			Scaff_and_chroms[line_list[0]].append(line_list[1].split(","))


all_scaff = Scaff_and_chroms.keys()

# Go throug all alignment positions
with open(input_file_name, "r") as input_file:
	temp_informations = []
	for line in input_file:
		line = line.rstrip()
		if line.startswith("/") or line.startswith("[") or line.startswith("NUC") or line =="":
			continue

		line_list = line.split("\t")
		scaff = line_list[-1]
		chrom = line_list[-2]
		if scaff in all_scaff:
			for all_chroms in Scaff_and_chroms[scaff]:
				if chrom in all_chroms:
					scaff_length = line_list[8]
					start = min(int(line_list[2]), int(line_list[3]))
					stop = max(int(line_list[2]), int(line_list[3]))
					temp_informations.append([scaff, start, stop, chrom])
					
# Save the info of the chromosomes
with open(out_name_prefix+".txt", "w") as out:
	temp2 = sorted(temp_informations)
	for element in temp2:
		print(element[0], element[1], element[2], element[3], sep="\t", file=out)


with open(out_name_prefix+".Potential.txt", "w") as out:

	# Summary to identify break points (for first scaffold)
	for scaff in all_scaff: 
		all_chroms = Scaff_and_chroms[scaff]
		# each element of temmp2: scaffold, start, stop, chromosome
		# We look for the first and second chromosome for the considered scaffold
		
		# Get first occurrence of scaffold
		for idx, line in enumerate(temp2): 
			if line[0] == scaff: 
				i = idx
				break
		considered_chrom = temp2[idx][3]

		# We then get the alternative chrom:
		for chroms in all_chroms:
			if considered_chrom in chroms:
					if chroms[0] != considered_chrom:
						alternative_chrom = chroms[0]
					else:
						alternative_chrom = chroms[1]

		print(considered_chrom, alternative_chrom)

		# Get split positions for each chromosomes
		# This consist in adding or removing units at each alignment, depending on the chromosome
		# that is matching 

		list_numeric = [] 
		count = 0
		for element in temp2:
			if element[0] != scaff:
				continue

			if element[3] == considered_chrom:
				count += 1
				list_numeric.append([count, element[0], element[1], element[2], element[3]])
				#print(count, element[0], element[1], element[2], element[3])
			else:
				count -=1 
				list_numeric.append([count, element[0], element[1], element[2], element[3]])
				#print(count, element[0], element[1], element[2], element[3]) 

		# Get maximum
		maximum_split = max(list_numeric)	

		newi = 0
		while newi < len(list_numeric):
			if list_numeric[newi] == maximum_split:

				if newi+1 < len(list_numeric):
					print(list_numeric[newi][1],  "Split 1;", list_numeric[newi][3],  "Split 2:", list_numeric[newi+1][2])
					print(list_numeric[newi][1],  list_numeric[newi][3], list_numeric[newi+1][2], sep="\t", file=out)
				else:
					print(list_numeric[newi][1],  "Split 1;", list_numeric[newi][3],  "Split 2: NA")
					print(list_numeric[newi][1],  list_numeric[newi][3], "NA", sep="\t", file=out)

				break

			newi += 1	



