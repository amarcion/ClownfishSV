import sys
from Bio import SeqIO

all_transcripts_file = sys.argv[1]
longest_transcript_file = sys.argv[2]
blastFile = sys.argv[3]
db_prefix = sys.argv[4]
output_file_name = sys.argv[5]

# Get longest transcripts (genes)
Genes = []
for record in SeqIO.parse(longest_transcript_file, "fasta"):
	Genes.append(record.id)

# Get all transcripts
Transcripts = []
Lengthtransc = {}
for record in SeqIO.parse(all_transcripts_file, "fasta"):
	Transcripts.append(record.id)
	Lengthtransc[record.id] = len(record.seq)*3

# Get Blast results
Dict_transcript_to_blast = {}
with open(blastFile, "r") as input_file:
	for line in input_file:
		line = line.rstrip()
		line_list = line.split("\t")
		Dict_transcript_to_blast[line_list[0]] = line_list[1]

# Get all information for mapping genes and annotation
# gene2Accession
Dict_Annotation_to_GeneID = {}
with  open(db_prefix+".gene2accession", "r") as input_file:
	for line in input_file:
		line = line.rstrip()
		line_list = line.split("\t")
		Dict_Annotation_to_GeneID[line_list[5]] = line_list[1]
# gene2go
Dict_GeneID_to_GO = {}
with open(db_prefix+".gene2go", "r") as input_file:
	for line in input_file:
		line = line.rstrip()
		line_list = line.split("\t")
		if line_list[1] not in Dict_GeneID_to_GO.keys():
			Dict_GeneID_to_GO[line_list[1]] = [line_list[2]]
		else:
			Dict_GeneID_to_GO[line_list[1]].append(line_list[2])
# gene_info
Dict_GeneID_to_Info = {}
with open(db_prefix+".gene_info", "r") as input_file:
	for line in input_file:
		line = line.rstrip()
		line_list = line.split("\t")
		Dict_GeneID_to_Info[line_list[1]] = [line_list[2], line_list[8]]

# Save annotation
with open(output_file_name, "w") as output_file:
	what_to_print = "TranscriptID\tIsLongestIsoform\tlength_bp\tNR_ID_Annotation\tGeneID\tGeneName\tGeneDescription\tGeneOntologies"
	print(what_to_print, file=output_file)

	for tran in Transcripts:
		# Longest transcript
		if tran in Genes:
			IsLongest="Yes"
		else:
			IsLongest="No"

		length_trans = str(Lengthtransc[tran])

		# NR_annotation
		if tran in Dict_transcript_to_blast.keys():
			NR_ID = Dict_transcript_to_blast[tran]
			if NR_ID in Dict_Annotation_to_GeneID.keys():
				GeneID = Dict_Annotation_to_GeneID[NR_ID]
				if GeneID in Dict_GeneID_to_Info.keys():
					GeneName = Dict_GeneID_to_Info[GeneID][0]
					GeneDescription = Dict_GeneID_to_Info[GeneID][1]
				else:
					GeneName = "NA"
					GeneDescription = "NA"

				if GeneID in Dict_GeneID_to_GO.keys():
					GO = ",".join(Dict_GeneID_to_GO[GeneID])
				else:
					GO = "NA"

			else:
				GeneID = "NA"	
				GeneName = "NA"
				GeneDescription = "NA"
				GO = "NA"
		else:
			NR_ID = "NA"
			GeneID = "NA"
			GeneName = "NA"
			GeneDescription = "NA"
			GO = "NA"

		#print annotation
		what_to_print = tran + "\t" + IsLongest +"\t" + length_trans + "\t" + NR_ID + "\t" + GeneID + "\t" + GeneName + "\t" + GeneDescription + "\t" + GO
		print(what_to_print, file=output_file)


