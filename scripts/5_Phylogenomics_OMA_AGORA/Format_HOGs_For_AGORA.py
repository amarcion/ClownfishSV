import os
import pyham
import logging

os.chdir("/Users/amarcion/Documents/PacBio_Assemblies_and_StructuralVariants/7_OMA_and_AncestralGenome")


# Get species tree
tree_str = pyham.utils.get_newick_string("0c_Trees/ManualSpeciesTree.WithAncestral.nwk", type="nwk")
orthoxml_path =  "1c_OMA_Output/HierarchicalGroups.orthoxml"
# pyham.Ham is the main object that containes all information and functionalities.
ham_analysis = pyham.Ham(tree_str, orthoxml_path, use_internal_name=True)


# Get genomes and genes for each species
species = ["FRE", "EPH", "CLA", "AKA", "SAN", "PRD", "SEB", "POL", "OMA", "ALL", "LAT", "CRP", 
           "MCC", "AKY", "LAZ", "PRC", "AMPPE", "OCE", "AMPOC", "BIA", "ACH", "DTR", "ORENI", "OREAU"]

Dict_genes_extantSpecies = {}
for sp in species:
    Dict_genes_extantSpecies[sp] = ham_analysis.get_extant_genome_by_name(sp).genes
print(Dict_genes_extantSpecies["AKA"][1:10])

# Transform the dictionary to get each gene as key and the species as value for later access:
Dict_genes = {}
for sp in Dict_genes_extantSpecies.keys():
    for gene in Dict_genes_extantSpecies[sp]:
        if gene not in Dict_genes.keys():
            Dict_genes[gene] = sp
        else:
            print(gene)

# See the internal names
# print("Ancestral genomes name using artificial ham names:")
# for ag in ham_analysis.taxonomy.internal_nodes:
#     print("\t- {}".format(ag.name))

# # We then get the list of HOGs for each Ancestral state
for ag in ham_analysis.taxonomy.internal_nodes:
	print(ag.name)
	ancestral_genome = ham_analysis.get_ancestral_genome_by_name(ag.name)
	ancestral_genes = ancestral_genome.genes
	print("number genes at node:", ag.name, len(ancestral_genes))

	Temp_Dict = {}
	dict_ancestral_Clust = ancestral_genome.get_ancestral_clustering()
	for HOG in dict_ancestral_Clust.keys():
		HOGnameComplete = str(HOG)
		index1=HOGnameComplete.find(":")
		index2=HOGnameComplete.find("_")
		HOGname = "HOG"+HOGnameComplete[index1+1:index2].strip("0")
		Temp_Dict[HOGname] = []

		for gene in dict_ancestral_Clust[HOG]:
			species = Dict_genes[gene]
			gene_id = str(gene)[5:-1]
			geneid_species = ham_analysis.get_gene_by_id(gene_id).get_dict_xref()["protId"]
			if species == "AMPOC" or species == "AMPPE" or species == "OREAU" or species == "ORENI" :
				newID = geneid_species.split()[0]
			elif species == "ACH":
				newID = species + "." + geneid_species.split()[0]
			else:
				newID = species + "." + geneid_species

			if ".t" in newID:
				index = newID.rfind(".t")
				finalID = newID[:index]
			else:
				finalID = newID
				
			Temp_Dict[HOGname].append(finalID)

		Temp_Dict[HOGname].sort()


	with open("3_HOGs_for_AGORA/orthologyGroups."+str(ag.name)+".list", "w") as outfile_NoHOG:
		with open("3_HOGs_for_AGORA/orthologyGroups.HOGnames."+str(ag.name)+".list", "w") as outfile_HOG:
			for HOG in Temp_Dict.keys():
				what_to_print_HOG = HOG + "\t" +  " ".join(Temp_Dict[HOG])
				what_to_print_NoHOG = " ".join(Temp_Dict[HOG])
				print(what_to_print_HOG, file=outfile_HOG)
				print(what_to_print_NoHOG, file=outfile_NoHOG)       










