import sys

stat_file = sys.argv[1]
GeneKept = sys.argv[2]
out_stat_file = sys.argv[3]

# Genes to keep
input_file = open(GeneKept, "r")
GenesKeep = [gene.rstrip() for gene in input_file.readlines()]
input_file.close()

#All transcripts
transKeep = []
for gene in GenesKeep:
	index = gene.find(".")
	trans = gene[:index]
	transKeep.append(trans)


# Filter stat file
input_file = open(stat_file, "r")
out_file = open(out_stat_file, "w")

for line in input_file:
	line = line.rstrip()
	if line.startswith("Transcript"):
		newline = line + "\tKeptFinalTranscripts\tKeptFinalGene"
		print(newline, file=out_file)
		continue
	
	line_list = line.split("\t")
	if line_list[1] in transKeep:
		var1 = "Yes"
	else:
		var1 = "No"

	if line_list[0] in GenesKeep:
		var2 = "Yes"
	else:
		var2 = "No"

	newline = line + "\t" + var1 + "\t" + var2
	print(newline, file=out_file)

input_file.close()
out_file.close()