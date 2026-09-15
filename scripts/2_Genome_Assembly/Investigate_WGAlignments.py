import sys


list_scaffolds_files = sys.argv[1]
align_file = sys.argv[2]
output_prefix = sys.argv[3]
threshold_percentAlign = int(sys.argv[4])

# Get all scaffolds
with open(list_scaffolds_files, "r") as input_scaff:
	list_of_scaffolds = [line.rstrip() for line in input_scaff]

# Get all match with chromosomes

Dictionary_align_chrom = {}	
Dictionary_align_scaff = {}

Dictionary_length_chrom = {}
Dictionary_length_scaff = {}

with open(align_file, "r") as input_align:

	for line in input_align:
		line_list = line.rstrip().split("\t")

		scaff = line_list[-1]
		chrom = line_list[-2]
		length_chrom = line_list[7]
		length_scaff = line_list[8]
		align_chrom = line_list[4]
		align_scaff = line_list[5]


		# Store alignment chromosomes and all alignment length
		if chrom not in Dictionary_align_chrom.keys():
			Dictionary_align_chrom[chrom] = {}
		if scaff not in Dictionary_align_chrom[chrom]:
			Dictionary_align_chrom[chrom][scaff] = [int(align_chrom)]
		else:
			Dictionary_align_chrom[chrom][scaff].append(int(align_chrom))

		# Store alignment scaffolds and all alignment length
		if scaff not in Dictionary_align_scaff.keys():
			Dictionary_align_scaff[scaff] = {}
		if chrom not in Dictionary_align_scaff[scaff]:
			Dictionary_align_scaff[scaff][chrom] = [int(align_scaff)]
		else:
			Dictionary_align_scaff[scaff][chrom].append(int(align_scaff))


		# Store scaffolds and chromosome length
		if scaff not in Dictionary_length_scaff.keys():
			Dictionary_length_scaff[scaff] = int(length_scaff)
		if chrom not in Dictionary_length_chrom.keys():
			Dictionary_length_chrom[chrom] = int(length_chrom)


# Check if some scaffolds do not align
Not_aligned = []
print(f"Total number of scaffolds: {len(list_of_scaffolds)}")
print(f"Number of scaffolds aligning: {len(Dictionary_align_scaff.keys())}")
for scaffold in list_of_scaffolds:
	if scaffold not in Dictionary_align_scaff.keys():
		Not_aligned.append(scaffold)
print(f"{len(Not_aligned)} scaffold not aligned")
print(Not_aligned)

print()

# Save all alignments of scaffolds to chromosomes
with open(output_prefix+".Scaffolds_to_Chromosomes.All_Alignments.txt", "w") as output_file:
	print("Chr\tScaff\tNbHits\tTotal_Align\tPropAlign", file=output_file)
	Dictionary_Chrom_LargeAlign = {}
	for chrom in Dictionary_align_chrom.keys():
		Dictionary_Chrom_LargeAlign[chrom] = []
		for scaff in Dictionary_align_chrom[chrom]:
			nb_hits = len(Dictionary_align_chrom[chrom][scaff])
			tot_align = sum(Dictionary_align_chrom[chrom][scaff])
			prop_align = tot_align/Dictionary_length_chrom[chrom]*100
			print(chrom, scaff, nb_hits, tot_align, prop_align, sep="\t", file=output_file)
			if prop_align > threshold_percentAlign:
				Dictionary_Chrom_LargeAlign[chrom].append(scaff)

# Save only best (> threshold_percentAlign) alignments of scaffolds to chromosomes
with open(output_prefix+".Scaffolds_to_Chromosomes.Best_Alignments."+str(threshold_percentAlign)+"Pct.txt", "w") as output_file:
	print("Chr\tNb_Best_Hits\tBest_scaffold_Hits", file=output_file)
	tot_scaff = []
	for chrom in Dictionary_Chrom_LargeAlign.keys():
		tot_scaff += Dictionary_Chrom_LargeAlign[chrom]
		print(chrom, len(Dictionary_Chrom_LargeAlign[chrom]), ",".join(Dictionary_Chrom_LargeAlign[chrom]), sep="\t", file=output_file)

print(f"Total number of scaffolds with good align to chrom: {len(tot_scaff)}")
print(f"Total number of scaffolds with good align to chrom: {len(set(tot_scaff))}")
print("If the two number above differ, some scaffolds have best hits in different chromosomes ")

if len(tot_scaff) != len(set(tot_scaff)):
	new_list = []
	double = []
	for scaff in tot_scaff:
		if scaff not in new_list:
			new_list.append(scaff)
		else:
			double.append(scaff)
	print("The scaffolds matching more chromosoems are: ", ",".join(double))
print()

# Save all alignments of chromosomes to scaffolds

with open(output_prefix+".Chromosomes_to_Scaffolds.All_Alignments.txt", "w") as output_file:
	print("Scaff\tChrom\tNbHits\tTotal_Align\tPropAlign", file=output_file)
	Dictionary_Scaffolds_LargeAlign = {}
	for scaff in Dictionary_align_scaff.keys():
		Dictionary_Scaffolds_LargeAlign[scaff] = []
		for chrom in Dictionary_align_scaff[scaff]:
			nb_hits = len(Dictionary_align_scaff[scaff][chrom])
			tot_align = sum(Dictionary_align_scaff[scaff][chrom])
			prop_align = tot_align/Dictionary_length_scaff[scaff]*100
			print(scaff, chrom, nb_hits, tot_align, prop_align, sep="\t", file=output_file)
			if prop_align > threshold_percentAlign:
				Dictionary_Scaffolds_LargeAlign[scaff].append(chrom)


# Check only best (> threshold_percentAlign) alignments of chromosomes to scaffolds
with open(output_prefix+".Chromosomes_to_Scaffolds.Best_Alignments."+str(threshold_percentAlign)+"Pct.txt", "w") as output_file:
	print("Scaff\tNb_Best_Hits\tBest_chr_Hits", file=output_file)
	for scaff in Dictionary_Scaffolds_LargeAlign.keys():
		print(scaff, len(Dictionary_Scaffolds_LargeAlign[scaff]), ",".join(Dictionary_Scaffolds_LargeAlign[scaff]), sep="\t", file=output_file)

# Check scaffolds not in the "Best Hits"
with open(output_prefix+".Scaffolds_NotBestHits."+str(threshold_percentAlign)+"Pct.txt", "w") as output_file:
	Missing_scaffolds = []
	for scaff in list_of_scaffolds:
		if scaff not in tot_scaff:
			print(scaff, file=output_file)
			Missing_scaffolds.append(scaff)

print(f"Number of scaffolds without good matches at {threshold_percentAlign}%: {len(Missing_scaffolds)}")
print(f"The scaffolds are {','.join(Missing_scaffolds)}")
print()




