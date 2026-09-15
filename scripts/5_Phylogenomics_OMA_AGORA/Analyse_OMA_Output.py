import os
import pyham
import logging


def build_HOG_dictionary(genome, Dict_genes, ham_analysis):
    """
    Build a dictionary of HOGs for a given genome.
    """
    hog_dict = {}
    ancestral_clusters = genome.get_ancestral_clustering()
    
    for HOG, genes in ancestral_clusters.items():
        HOGnameComplete = str(HOG)
        index1 = HOGnameComplete.find(":")
        index2 = HOGnameComplete.find("_")
        HOGname = "HOG" + HOGnameComplete[index1+1:index2].strip("0")
        hog_dict[HOGname] = []
        
        for gene in genes:
            species = Dict_genes[gene]
            gene_id = str(gene)[5:-1]
            geneid_species = ham_analysis.get_gene_by_id(gene_id).get_dict_xref()["protId"]
            
            if species in ["AMPOC", "AMPPE", "OREAU", "ORENI"]:
                newID = geneid_species.split()[0]
            elif species == "ACH":
                newID = species + "." + geneid_species.split()[0]
            else:
                newID = species + "." + geneid_species
                
            hog_dict[HOGname].append(newID)
        
        hog_dict[HOGname].sort()
    
    return hog_dict

def write_HOG_files(hog_dict, prefix, outdir="."):
    """
    Write HOG dictionary to two files inside the specified directory:
    - <outdir>/HOGgenes_<prefix>.txt (with HOG names)
    - <outdir>/HOGgenes_<prefix>.NoHOGname.txt (without HOG names)

    Parameters:
        hog_dict (dict): Dictionary of HOG -> list of genes
        prefix (str): Prefix for the file names
        outdir (str): Directory to save the files (default: current directory)
    """
    os.makedirs(outdir, exist_ok=True)

    file_with_HOG = os.path.join(outdir, f"HOGgenes_{prefix}.txt")
    file_no_HOG = os.path.join(outdir, f"HOGgenes_{prefix}.NoHOGname.txt")

    with open(file_with_HOG, "w") as outfile_HOG, open(file_no_HOG, "w") as outfile_NoHOG:
        for hog_name, genes in hog_dict.items():
            line_with_HOG = f"{hog_name}\t{' '.join(genes)}"
            line_no_HOG = " ".join(genes)
            print(line_with_HOG, file=outfile_HOG)
            print(line_no_HOG, file=outfile_NoHOG)

out_file_prefix = "2_OMA_HOGs_Summaries/HOGs_AncestralNodes"
os.chdir("/Users/amarcion/Documents/PacBio_Assemblies_and_StructuralVariants/7_OMA_and_AncestralGenome")

# Get species tree
tree_str = pyham.utils.get_newick_string("0c_Trees/ManualSpeciesTree.nwk", type="nwk")
orthoxml_path =  "1c_OMA_Output/HierarchicalGroups.orthoxml"
# pyham.Ham is the main object that containes all information and functionalities.
ham_analysis = pyham.Ham(tree_str, orthoxml_path, use_internal_name=False)

# See the internal names
#print("Ancestral genomes name using artificial ham names:")
#for ag in ham_analysis.taxonomy.internal_nodes:
#    print("\t- {}".format(ag.name))
#print(pyham.utils.previsualize_taxonomy(tree_str))

# Last common ancestor of all considered species:
# FRE/EPH/CLA/AKA/SAN/PRD/SEB/POL/OMA/ALL/LAT/CRP/MCC/AKY/LAZ/PRC/AMPPE/OCE/AMPOC/BIA/ACH/DTR/ORENI/OREAU
outgroups_genome = ham_analysis.get_ancestral_genome_by_name("FRE/EPH/CLA/AKA/SAN/PRD/SEB/POL/OMA/ALL/LAT/CRP/MCC/AKY/LAZ/PRC/AMPPE/OCE/AMPOC/BIA/ACH/DTR/ORENI/OREAU")

# Last common ancestor of damselfish:
# FRE/EPH/CLA/AKA/SAN/PRD/SEB/POL/OMA/ALL/LAT/CRP/MCC/AKY/LAZ/PRC/AMPPE/OCE/AMPOC/BIA/ACH/DTR
Damsel_genome = ham_analysis.get_ancestral_genome_by_name("FRE/EPH/CLA/AKA/SAN/PRD/SEB/POL/OMA/ALL/LAT/CRP/MCC/AKY/LAZ/PRC/AMPPE/OCE/AMPOC/BIA/ACH/DTR")

# Last common ancestor of clownfish:
# FRE/EPH/CLA/AKA/SAN/PRD/SEB/POL/OMA/ALL/LAT/CRP/MCC/AKY/LAZ/PRC/AMPPE/OCE/AMPOC/BIA
Clown_genome = ham_analysis.get_ancestral_genome_by_name("FRE/EPH/CLA/AKA/SAN/PRD/SEB/POL/OMA/ALL/LAT/CRP/MCC/AKY/LAZ/PRC/AMPPE/OCE/AMPOC/BIA")

# Get ancestra genes for the three level
outgroups_ancestral_genes = outgroups_genome.genes
Damsel_ancestral_genes = Damsel_genome.genes
Clown_ancestral_genes = Clown_genome.genes

print("Number ancestral genes considerin outgroups:", len(outgroups_ancestral_genes))
print("Number ancestral genes in Damselfish: ", len(Damsel_ancestral_genes))
print("Number ancestral genes in clownfish: ", len(Clown_ancestral_genes))

# Create the tree profile
treeprofile = ham_analysis.create_tree_profile(outfile="2_OMA_HOGs_Summaries/treeprofile_Outgroups.html")

# Get genomes and genes for each species
species = ["FRE", "EPH", "CLA", "AKA", "SAN", "PRD", "SEB", "POL", "OMA", "ALL", "LAT", "CRP", 
           "MCC", "AKY", "LAZ", "PRC", "AMPPE", "OCE", "AMPOC", "BIA", "ACH", "DTR", "ORENI", "OREAU"]
Dict_genes_extantSpecies = {}
for sp in species:
    Dict_genes_extantSpecies[sp] = ham_analysis.get_extant_genome_by_name(sp).genes

# Transform the dictionary to get each gene as key and the species as value for later access:
Dict_genes = {}
for sp in Dict_genes_extantSpecies.keys():
    for gene in Dict_genes_extantSpecies[sp]:
        if gene not in Dict_genes.keys():
            Dict_genes[gene] = sp
        else:
            print(gene)


# Get the ancestral cluster for each type of ancestral genomes (outgroups, damselfish, clowns). 
Dictionary_HOGs_Outgroups = build_HOG_dictionary(outgroups_genome, Dict_genes, ham_analysis)
Dictionary_HOGs_Damsel    = build_HOG_dictionary(Damsel_genome, Dict_genes, ham_analysis)
Dictionary_HOGs_Clown     = build_HOG_dictionary(Clown_genome, Dict_genes, ham_analysis)

# Write files for each type
write_HOG_files(Dictionary_HOGs_Outgroups, "Outgroups", out_file_prefix)
write_HOG_files(Dictionary_HOGs_Damsel, "Damselfish", out_file_prefix)
write_HOG_files(Dictionary_HOGs_Clown, "Clownfish", out_file_prefix)




