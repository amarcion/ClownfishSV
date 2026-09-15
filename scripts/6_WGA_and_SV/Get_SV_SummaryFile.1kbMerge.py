import sys

vcf_input_file_name = sys.argv[1]
chrom = sys.argv[2]
SVType = sys.argv[3]
overlap_thresh = sys.argv[4]
output_file_prefix =  sys.argv[5]


def get_SV_length_from_INFO(Position, INFO_column_vcf):
	index1= INFO_column_vcf.find("END=")+4
	index2=  INFO_column_vcf.find(";EndB=")
	stop=int(INFO_column_vcf[index1:index2])
	SV_length = abs(stop-int(Position))
	
	return SV_length
#END

def read_vcf(vcf_input_file_name):

	Dictionary_of_SVs = {}

	with open(vcf_input_file_name, "r") as infile:
		for line in infile:
			line = line.rstrip()
			if line.startswith("##"):
				continue
			elif line.startswith("#"):
				species_list = line.split("\t")[9:]	# Species are from 9 columns
				# Additional columns: 0: CHROM, 1: POS, 2: ID, 4: SV_Type, 7: INFO (where length is)
				continue

			line_list = line.split("\t")
			SV_len = get_SV_length_from_INFO(line_list[1], line_list[7])
			genotypes = line_list[9:]
			species_with_SV = [val2 for val1, val2 in zip(genotypes, species_list) if val1 == '1']
			Dictionary_of_SVs[int(line_list[1])] = [line_list[2], line_list[4], SV_len, species_with_SV]

	return Dictionary_of_SVs, species_list
#END

def merge_overlap_SVs(dictonary_SV, overlap_threshold):
	Dictionary_of_SVs_overlap = {}
	List_of_merged_SVs = []

	sorted_positions = sorted(dictonary_SV.keys())
	counter = 0
	for pos in sorted_positions:
		value_min = pos-overlap_threshold
		value_max = pos+overlap_threshold

		flag = True
		needed_pos = 0
		for new_pos in Dictionary_of_SVs_overlap.keys():
			if new_pos >= value_min and new_pos <= value_max:
				needed_pos = new_pos
				flag = False

		if flag:
			counter += 1
			Dictionary_of_SVs_overlap[pos] = dictonary_SV[pos]
			List_of_merged_SVs.append([pos]+dictonary_SV[pos][:-1]+["Kept_"+str(counter)]+["Group_"+str(counter)]+[dictonary_SV[pos][-1]])

		else:
			Dictionary_of_SVs_overlap[needed_pos][-1] =  Dictionary_of_SVs_overlap[needed_pos][-1] + dictonary_SV[pos][-1]
			List_of_merged_SVs.append([pos]+dictonary_SV[pos][:-1]+["Removed_"+str(counter)]+["Group_"+str(counter)]+[dictonary_SV[pos][-1]])

	return(Dictionary_of_SVs_overlap, List_of_merged_SVs)
#END

def remove_single_species_SVs(dictionary_SV):
    dict_more_species = {}
    dict_single_species = {}

    for SV, data in dictionary_SV.items():
        species = data[-1]
        
        # Remove 'PRCRefNo24' if present, without modifying original list
        filtered_species = [sp for sp in species if sp != 'PRCRefNo24']

        if len(set(filtered_species)) > 1:
            dict_more_species[SV] = data
        else:
            dict_single_species[SV] = data

    return dict_more_species, dict_single_species
#End

def format_line(list_to_write, chrom, pos, SVCat=False):
	if not SVCat:
		list_to_write_str = [str(el) for el in list_to_write[:3]]
		what_to_print = chrom + "\t" + str(pos) + "\t" + "\t".join(list_to_write_str)
		species = ",".join(list_to_write[-1])
		what_to_print = what_to_print + "\t" + species

	else:
		list_to_write_str = [str(el) for el in list_to_write[:3]]
		what_to_print = chrom + "\t" + str(pos) + "\t" + "\t".join(list_to_write_str) + "\t" + list_to_write[-1]
		species = ",".join(list_to_write[-2])
		what_to_print = what_to_print + "\t" + species
	return what_to_print
#End

def count_PRC_SVs(dictionary_SVs):
	Dict_count = {"PRC": 0, "PRCRefNo24": 0, "Both": 0}
	for SV in dictionary_SVs.keys():
		if "PRC" in dictionary_SVs[SV][-1] and "PRCRefNo24" in dictionary_SVs[SV][-1]:
			Dict_count["Both"] += 1
		elif "PRC" in dictionary_SVs[SV][-1]:
			Dict_count["PRC"] += 1
		elif "PRCRefNo24" in dictionary_SVs[SV][-1]:
			Dict_count["PRCRefNo24"] += 1

	return Dict_count
#END

def get_SV_cat(dictionary_SVs, dict_singlespecies, dict_multispecies):
	Dict_cat = {}
	for SV in dictionary_SVs.keys():
		if SV in dict_singlespecies.keys():
			Dict_cat[SV] = "SingleSpecies"
		elif SV in dict_multispecies.keys():
			Dict_cat[SV] = "MultiSpecies"
		else:
			print(SV, "Not found")
			Dict_cat[SV] = ""
	return Dict_cat
#END

def count_type_SV_per_species(species, dictionary_SVs, Dict_cat):
	subdic_cat = {"SingleSpecies": 0, "MultiSpecies": 0}
	count_total = 0 
	for SV in dictionary_SVs.keys():
		if species in dictionary_SVs[SV][-1]:
			count_total += 1
			subdic_cat[Dict_cat[SV]] += 1
	return subdic_cat, count_total
#END

print("Analysing chrom:", chrom, ", SVtype:", SVType, "Overlap threshold:", overlap_thresh )

# Get SVs
dict_of_SVs, species_list  = read_vcf(vcf_input_file_name)
# species_list : species considered in this chromsomes
# dict_of_SVs: Dictioanry of SV retrieved from the VCF. it has as key the position and as value the list [SV_ID, SV_Type, SV_len, species_with_SV]

# Remove overlapped SV
dict_of_SVs_overlapped, liste_removedOverlap = merge_overlap_SVs(dict_of_SVs, int(overlap_thresh))
# We then get the dict_of_SVs_overlapped, which containes the SV after merging overlapping sv (overlapping to a given threshold)
# The liste_removedOverlap report the list of SVs that were merged in dict_of_SVs_overlapped and thus not present anymore


# Get SV that are multispecies and the one that are present only in single species (PRC and PRCRef considered as signle species)
dict_of_SVs_overlapped_multspecies, dict_of_SVs_overlapped_singlespecies = remove_single_species_SVs(dict_of_SVs_overlapped)

# get dictionary of categories for each variant
dic_cat = get_SV_cat(dict_of_SVs_overlapped, dict_of_SVs_overlapped_singlespecies, dict_of_SVs_overlapped_multspecies)

# Count Nb consistent PRC SV
Dic_PRC = count_PRC_SVs(dict_of_SVs_overlapped)


######
# SAVE THE RESULTS
######


# Save the results 
# 1. All SV summary with category
fname = output_file_prefix + "."+chrom+"."+SVType+"."+overlap_thresh+"trsh.AllSVs.txt"
with open(fname, "w") as outfile:
	print("Chrom\tPosition\tID\tSVType\tSVLength\tSVCat\tSpecies", file=outfile)
	for SV in dict_of_SVs_overlapped.keys():
		cat_liste = [x for x in dict_of_SVs_overlapped[SV]]
		cat_liste.append(dic_cat[SV])
		line = format_line(cat_liste, chrom, SV, SVCat=True)
		print(line, file=outfile)

# 2. Multi-species SVs in file: output_file_prefix + "."+chr+"."+SVType+".MultiSpecies.txt"
fname = output_file_prefix + "."+chrom+"."+SVType+"."+overlap_thresh+"trsh.MultiSpecies.txt"
with open(fname, "w") as outfile:
	print("Chrom\tPosition\tID\tSVType\tSVLength\tSpecies", file=outfile)
	for SV in dict_of_SVs_overlapped_multspecies.keys():
		line = format_line(dict_of_SVs_overlapped_multspecies[SV],chrom, SV)
		print(line, file=outfile)

# 3. Single species SVs in file: output_file_prefix + "."+chr+"."+SVType+".SingleSpecies.txt"
fname = output_file_prefix + "."+chrom+"."+SVType+"."+overlap_thresh+"trsh.SingleSpecies.txt"
with open(fname, "w") as outfile:
	print("Chrom\tPosition\tID\tSVType\tSVLength\tSpecies", file=outfile)
	for SV in dict_of_SVs_overlapped_singlespecies.keys():
		line = format_line(dict_of_SVs_overlapped_singlespecies[SV],chrom, SV)
		print(line, file=outfile)

# 4. Summary file 
fname = output_file_prefix + "."+chrom+"."+SVType+"."+overlap_thresh+"trsh.OverallSummary.txt"
with open(fname, "w") as outfile:
	print("Number of analyzed species:", len(species_list), file=outfile)
	print("Analyzed species:", ",".join(species_list), file=outfile)
	print("Original number of SV:", len(dict_of_SVs.keys()), file=outfile)
	print("Number of SV without overlap:", len(dict_of_SVs_overlapped.keys()), file=outfile)
	print("Number of SV single species:", len(dict_of_SVs_overlapped_singlespecies.keys()), file=outfile)
	print("Number of SV multi species:", len(dict_of_SVs_overlapped_multspecies.keys()), file=outfile)
	print("Number of SV in PRC and PRCRefNo24", str(Dic_PRC["Both"]), file=outfile)
	print("Number of SV in PRC only", str(Dic_PRC["PRC"]), file=outfile)
	print("Number of SV in PRCRefNo24 only", str(Dic_PRC["PRCRefNo24"]), file=outfile)

# 5. Removed Overlap
fname = output_file_prefix + "."+chrom+"."+SVType+"."+overlap_thresh+"trsh.OverlappedSVs.txt"
with open(fname, "w") as outfile:
	print("Chrom\tPosition\tID\tSVType\tSVLength\tKeptRemoved\tGroupKeptRemoved\tSpecies", file=outfile)
	for el in liste_removedOverlap:
		list_to_write=el[1:]
		pos=el[0]
		list_to_write_str = [str(el) for el in list_to_write[:-1]]
		what_to_print = chrom + "\t" + str(pos) + "\t" + "\t".join(list_to_write_str)
		species = ",".join(list_to_write[-1])
		what_to_print = what_to_print + "\t" + species
		print(what_to_print, file=outfile)


# 6. Number of SV per species
fname = output_file_prefix + "."+chrom+"."+SVType+"."+overlap_thresh+"trsh.SpeciesSummary.txt"
with open(fname, "w") as outfile:
	print("Species\tChrom\tSVType\tTotalCount\tSingleSpecies\tMultiSpecies", file=outfile)
	for sp in species_list:
		temp_cat, tmp_tot = count_type_SV_per_species(sp, dict_of_SVs_overlapped, dic_cat)
		what_to_print = sp +"\t"+chrom +"\t"+ SVType +"\t"+ str(tmp_tot) +"\t"+ str(temp_cat["SingleSpecies"]) \
		+"\t"+ str(temp_cat["MultiSpecies"])
		print(what_to_print, file=outfile)