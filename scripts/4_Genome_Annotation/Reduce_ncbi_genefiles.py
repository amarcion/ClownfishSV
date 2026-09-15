import sys

blastFile = sys.argv[1]
path_to_db_files = sys.argv[2]
prefix_to_out_db_files = sys.argv[3]

# Get needed prot
List_references = []

with open(blastFile, "r") as input_file:
	input_file = open(blastFile, "r")
	for line in input_file:
		line = line.rstrip()
		line_list = line.split("\t")
		List_references.append(line_list[1])

List_references = set(List_references)
print("Nb reference sequences: ", len(List_references))

# keep only needed info in gene2accession
temp_gene2acc = []
temp_geneID = []

with open(prefix_to_out_db_files+".gene2accession", "w") as output_file:
	with open(path_to_db_files+"gene2accession", "r") as input_file:
		for line in input_file:
			line = line.rstrip()
			if line.startswith("#"):
				print(line, file=output_file)
				continue
			line_list = line.split("\t")
			if line_list[5] in List_references:
				temp_gene2acc.append(line_list[5])
				temp_geneID.append(line_list[1])
				print(line, file=output_file)

print("Nb total gene2accessions retrived: ", len(temp_gene2acc))
print("Nb total gene id retrieved: ", len(temp_geneID))
temp_gene2acc = set(temp_gene2acc)
temp_geneID = set(temp_geneID)
print("Nb unique gene2accessions retrived: ", len(temp_gene2acc))
print("Nb unique gene id retrieved: ", len(temp_geneID))

# keep only needed info in gene2accession
tmp_GO_acc = []

with open(prefix_to_out_db_files+".gene2go", "w") as output_file:
	with open(path_to_db_files+"gene2go", "r") as input_file:
		for line in input_file:
			line = line.rstrip()
			if line.startswith("#"):
				print(line, file=output_file)
				continue
			line_list = line.split("\t")
			if line_list[1] in temp_geneID:
				tmp_GO_acc.append(line_list[1])
				print(line, file=output_file)

print("Nb of GO term retrieved: ", len(tmp_GO_acc))
tmp_GO_acc = set(tmp_GO_acc)
print("Nb of retrieved genes with GO: ", len(tmp_GO_acc))

# keep only needed info in gene_info
temp_gene_info = []

with open(prefix_to_out_db_files+".gene_info", "w") as output_file:
		with open(path_to_db_files+"gene_info", "r") as input_file:
		for line in input_file:
			line = line.rstrip()
			if line.startswith("#"):
				print(line, file=output_file)
				continue
			line_list = line.split("\t")
			if line_list[1] in temp_geneID:
				temp_gene_info.append(line_list[1])
				print(line, file=output_file)

print("Nb of genes with gene_info: ", len(temp_gene_info))
temp_gene_info = set(temp_gene_info)
print("Nb unique genes with gene_info: ", len(temp_gene_info))















