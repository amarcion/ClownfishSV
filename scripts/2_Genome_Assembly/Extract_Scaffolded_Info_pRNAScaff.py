import sys
import os
from Bio import SeqIO

SpeciesFile = sys.argv[1]
ChromFile = sys.argv[2]
path_to_pRNA_outFolder = sys.argv[3]
path_to_out = sys.argv[4]
step = sys.argv[5]

def get_species(species_file):
	input_file = open(species_file, "r")
	species = [line.rstrip() for line in input_file]
	input_file.close()
	return(species)

def get_expected_chrom(expectedChrom_file):
	Dictionary_results = {}
	with open(expectedChrom_file, "r") as input_file:
		for line in input_file:
			line = line.rstrip()
			if line.startswith("Species"):
				continue

			line_list = line.split("\t")
			if line_list[0] not in Dictionary_results.keys():
				Dictionary_results[line_list[0]] = {}

			Dictionary_results[line_list[0]][line_list[3]] = line_list[1]

	return(Dictionary_results)

def get_nb_scaffolds(fasta_file):
	nb_scaff = 0
	for record in SeqIO.parse(fasta_file, "fasta"):
		nb_scaff += 1
	return(nb_scaff)


def get_scaffolded_sequences(bothPath_file):
	Scaffold_path = []
	Scaffolds = []

	with open(bothPath_file, "r") as input_file:
		for line in input_file:
			line = line.rstrip()
			line_rp = line.replace("/r", "")
			line_rp = line_rp.replace("->N(100)->", " ")
			line_split = line_rp.split(" ")
			#print(line_split)
			#print(Line_split[::-1])
			Scaffold_path.append(tuple(line_split))
			Scaffolds += line_split

	Scaffolds = list(set(Scaffolds))
	return(Scaffold_path, Scaffolds)

def scaffPath_to_chromPath(scaffPath, species, link2Chrom):
	chromPath = []
	for el in scaffPath:
		if el in link2Chrom[species].keys():
			chrom = link2Chrom[species][el]
		else:
			chrom = "NA"

		chromPath.append(chrom)
	return(tuple(chromPath))

###

print("Species", "InitialScaffolds", "pe5000_Nbfinal", "pe7000_Nbfinal", "pe9000_Nbfinal", "pe11000_Nbfinal", "pe15000_Nbfinal", "peMerged")


Species = get_species(SpeciesFile)
ChromResults = get_expected_chrom(ChromFile)

if step=="step1":
	pe=["5000", "7000", "9000","11000","15000"]
else:
	pe=["5000", "7000", "9000","11000","15000", "Merged"]


for sp_considered in Species:

	Summary_path_results = {}
	Summary_path_results_Chrom = {}
	Link_scaffPath_chromPath = {}

	Scaffolded_sequences = {}
	Scaffolds_and_chrom = {}

	for pe_considered in pe:

		path_to_file = path_to_pRNA_outFolder+"/"+sp_considered+".pe"+pe_considered+".out/"

		Scaffolded, Scaffolds = get_scaffolded_sequences(path_to_file+"both.path")

		# Summary of the scaffold path
		for scaffpath in Scaffolded:
			chromPath = scaffPath_to_chromPath(scaffpath, sp_considered, ChromResults)
			
			if scaffpath in Summary_path_results.keys():
				Summary_path_results[scaffpath].append(pe_considered)
				Summary_path_results_Chrom[chromPath].append(pe_considered)

			elif scaffpath[::-1] in Summary_path_results.keys():
				Summary_path_results[scaffpath[::-1]].append(pe_considered)
				Summary_path_results_Chrom[chromPath[::-1]].append(pe_considered)

			else:
				Summary_path_results[scaffpath] = [pe_considered]
				Summary_path_results_Chrom[chromPath] = [pe_considered]
				Link_scaffPath_chromPath[scaffpath] = chromPath
		
		# Summary of the Sequences Scaffoled
		for scaff in Scaffolds:
			if scaff in Scaffolded_sequences:
				Scaffolded_sequences[scaff].append(pe_considered)

			else:
				Scaffolded_sequences[scaff] = [pe_considered]

				if scaff in ChromResults[sp_considered].keys():
					chrom = ChromResults[sp_considered][scaff]
				else:
					chrom = "NA"
				Scaffolds_and_chrom[scaff] = chrom
				

	if step=="step1":
		out_file_name = path_to_out +"/"+ sp_considered + ".SummaryPath.txt"
	else:
		out_file_name = path_to_out +"/"+ sp_considered + ".SummaryPath.S2.txt"
	with open(out_file_name, "w") as out_file:
		if step=="step1":
			what_to_print = "ScaffPath\tpe5000\tpe7000\tpe9000\tpe11000\tpe15000"
		else:
			what_to_print = "ScaffPath\tpe5000\tpe7000\tpe9000\tpe11000\tpe15000\tpeMerged"

		print(what_to_print, file=out_file)
		for scaffpath in Summary_path_results.keys():
			what_to_print = str(scaffpath)
			for el in pe:
				if el in Summary_path_results[scaffpath]:
					what_to_print += "\tYes"
				else:
					what_to_print += "\tNo"
			print(what_to_print, file=out_file)


	if step=="step1":
		out_file_name = path_to_out+"/"+sp_considered + ".SummaryPathWithChrom.txt"
	else:
		out_file_name = path_to_out+"/"+sp_considered + ".SummaryPathWithChrom.S2.txt"
	with open(out_file_name, "w") as out_file:
		if step=="step1":
			what_to_print = "ScaffPath\tChromPath\tpe5000\tpe7000\tpe9000\tpe11000\tpe15000"
		else:
			what_to_print = "ScaffPath\tChromPath\tpe5000\tpe7000\tpe9000\tpe11000\tpe15000\tpeMerged"
		
		print(what_to_print, file=out_file)
		for scaffpath in Summary_path_results.keys():
			chrom_path=Link_scaffPath_chromPath[scaffpath]
			what_to_print = str(scaffpath)+"\t"+str(chrom_path)
			for el in pe:
				if el in Summary_path_results[scaffpath]:
					what_to_print += "\tYes"
				else:
					what_to_print += "\tNo"
			print(what_to_print, file=out_file)

	if step=="step1":
		out_file_name = path_to_out+"/"+sp_considered + ".ScaffoldedChrom.txt"
	else:
		out_file_name = path_to_out+"/"+sp_considered + ".ScaffoldedChrom.S2.txt"
	with open(out_file_name, "w") as out_file:
		if step=="step1":
			what_to_print = "Scaffold\tPotentialChrom\tpe5000\tpe7000\tpe9000\tpe11000\tpe15000"
		else:
			what_to_print = "Scaffold\tPotentialChrom\tpe5000\tpe7000\tpe9000\tpe11000\tpe15000\tpeMerged"

		print(what_to_print, file=out_file)

		for scaff in Scaffolded_sequences.keys():
			chrom = Scaffolds_and_chrom[scaff]
			what_to_print = str(scaff)+"\t"+str(chrom)
			for el in pe:
				if el in Scaffolded_sequences[scaff]:
					what_to_print += "\tYes"
				else:
					what_to_print += "\tNo"
			print(what_to_print, file=out_file)

	# Get scaffolded sequences per species

	#initialNb = get_nb_scaffolds(path_to_pRNA_outFolder+"/"+sp_considered+".pRNA_scaffold.r1.fasta")
	pe5000_scaffolded = get_nb_scaffolds(path_to_pRNA_outFolder+"/"+sp_considered+".pe5000.out/scaffold.fasta")
	pe5000_unscaffolded = get_nb_scaffolds(path_to_pRNA_outFolder+"/"+sp_considered+".pe5000.out/unscaffold.fasta")
	pe5000_final = get_nb_scaffolds(path_to_pRNA_outFolder+"/"+sp_considered+".pe5000.out/P_RNA_scaffold.fasta")

	pe7000_scaffolded = get_nb_scaffolds(path_to_pRNA_outFolder+"/"+sp_considered+".pe7000.out/scaffold.fasta")
	pe7000_unscaffolded = get_nb_scaffolds(path_to_pRNA_outFolder+"/"+sp_considered+".pe7000.out/unscaffold.fasta")
	pe7000_final = get_nb_scaffolds(path_to_pRNA_outFolder+"/"+sp_considered+".pe7000.out/P_RNA_scaffold.fasta")

	pe9000_scaffolded = get_nb_scaffolds(path_to_pRNA_outFolder+"/"+sp_considered+".pe9000.out/scaffold.fasta")
	pe9000_unscaffolded = get_nb_scaffolds(path_to_pRNA_outFolder+"/"+sp_considered+".pe9000.out/unscaffold.fasta")
	pe9000_final = get_nb_scaffolds(path_to_pRNA_outFolder+"/"+sp_considered+".pe9000.out/P_RNA_scaffold.fasta")

	pe11000_scaffolded = get_nb_scaffolds(path_to_pRNA_outFolder+"/"+sp_considered+".pe11000.out/scaffold.fasta")
	pe11000_unscaffolded = get_nb_scaffolds(path_to_pRNA_outFolder+"/"+sp_considered+".pe11000.out/unscaffold.fasta")
	pe11000_final = get_nb_scaffolds(path_to_pRNA_outFolder+"/"+sp_considered+".pe11000.out/P_RNA_scaffold.fasta")

	pe15000_scaffolded = get_nb_scaffolds(path_to_pRNA_outFolder+"/"+sp_considered+".pe15000.out/scaffold.fasta")
	pe15000_unscaffolded = get_nb_scaffolds(path_to_pRNA_outFolder+"/"+sp_considered+".pe15000.out/unscaffold.fasta")
	pe15000_final = get_nb_scaffolds(path_to_pRNA_outFolder+"/"+sp_considered+".pe15000.out/P_RNA_scaffold.fasta")

	if step!="step1":
		peMerged_scaffolded = get_nb_scaffolds(path_to_pRNA_outFolder+"/"+sp_considered+".peMerged.out/scaffold.fasta")
		peMerged_unscaffolded = get_nb_scaffolds(path_to_pRNA_outFolder+"/"+sp_considered+".peMerged.out/unscaffold.fasta")
		peMerged_final = get_nb_scaffolds(path_to_pRNA_outFolder+"/"+sp_considered+".peMerged.out/P_RNA_scaffold.fasta")


	if step=="step1":
		out_file_name = path_to_out+"/"+sp_considered + ".ScaffoldingSummary.txt"
	else:
		out_file_name = path_to_out+"/"+sp_considered + ".ScaffoldingSummary.S2.txt"
	with open(out_file_name, "w") as out_file:
		
		print("File\tNb:ScaffoldedSeq\tNb_UnScaffoldedSeq\tNbFinalScaffolds", file=out_file)
		#print("OriginalFile", "0", "0", initialNb, sep="\t", file=out_file)
		print("pe5000", pe5000_scaffolded, pe5000_unscaffolded, pe5000_final, sep="\t", file=out_file)
		print("pe7000", pe7000_scaffolded, pe7000_unscaffolded, pe7000_final, sep="\t", file=out_file)
		print("pe9000", pe9000_scaffolded, pe9000_unscaffolded, pe9000_final, sep="\t", file=out_file)
		print("pe11000", pe11000_scaffolded, pe11000_unscaffolded, pe11000_final, sep="\t", file=out_file)
		print("pe15000", pe15000_scaffolded, pe15000_unscaffolded, pe15000_final, sep="\t", file=out_file)
		if step!="1":
			print("peMerged", peMerged_scaffolded, peMerged_unscaffolded, peMerged_final, sep="\t", file=out_file)
		
	if step=="step1":
		print(sp_considered, pe5000_final, pe7000_final, pe9000_final, pe11000_final, pe15000_final)
		#print(sp_considered, initialNb, pe5000_final, pe7000_final, pe9000_final, pe11000_final, pe15000_final)
	else:
		print(sp_considered, pe5000_final, pe7000_final, pe9000_final, pe11000_final, pe15000_final, peMerged_final)
		#print(sp_considered, initialNb, pe5000_final, pe7000_final, pe9000_final, pe11000_final, pe15000_final, peMerged_final)




