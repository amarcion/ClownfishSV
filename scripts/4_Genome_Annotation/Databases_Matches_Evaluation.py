from Bio import SeqIO
import sys

ProtOnly_LongestIso_fileName = sys.argv[1]
ProtOnly_SwissProt_fileName = sys.argv[2]
ProtOnly_InterProScan_fileName = sys.argv[3]
output_path = sys.argv[4]

def get_transcritps(input_aa_name):
	Transc = []
	for record in SeqIO.parse(input_aa_name, "fasta"):
		Transc.append(record.id)
	return(Transc)

def get_hits(input_hitsfiles_name):
	hitsTransc = []
	with open(input_hitsfiles_name, "r") as input_file:
		for line in input_file:
			line = line.rstrip()
			if "AntiFam" in line:
				continue
			line_list = line.split("\t")
			hitsTransc.append(line_list[0])

	hitsTransc = set(hitsTransc)
	hitsTransc = list(hitsTransc)
	return(hitsTransc)

def get_transcripts_annotations(transcripts, swissport, interproscan):
	Annotation_Dict = { "SwissOnly": 0, "InterProScan":0, "NoMatch":0}
	transcripts = set(transcripts)
	swissport = set(swissport)
	interproscan = set(interproscan)

	SwissAndInterProScan = set(transcripts) & set(swissport) & set(interproscan)
	Annotation_Dict["Swiss+InterProScan"] = len(SwissAndInterProScan)
	SwissOnly = set.difference((transcripts & swissport), SwissAndInterProScan)
	Annotation_Dict["SwissOnly"] = len(SwissOnly)
	InterProScanOnly = set.difference((transcripts & interproscan), SwissAndInterProScan)
	Annotation_Dict["InterProScan"] = len(InterProScanOnly)
	NoMatches = set.difference(transcripts, swissport,interproscan)
	Annotation_Dict["NoMatch"] = len(NoMatches)
	transcripts_to_keep = set.difference(transcripts, NoMatches)
	transcripts_to_keep = list(transcripts_to_keep)
	return(Annotation_Dict, transcripts_to_keep)
			

# Get information and stats for ProtOnly transcripts 
Transc_ProtOnly = get_transcritps(ProtOnly_LongestIso_fileName)
SwissProt_ProtOnly = get_hits(ProtOnly_SwissProt_fileName)
InterProScan_ProtOnly = get_hits(ProtOnly_InterProScan_fileName)

print("ProtOnly: number of trancripts: ", len(Transc_ProtOnly))
print("ProtOnly: number of trancripts with swissprot hits: ", len(SwissProt_ProtOnly))
print("ProtOnly: number of trancripts with interproscan hits: ", len(InterProScan_ProtOnly))

ProtOnly_res, FilteredTrans_ProtOnly = get_transcripts_annotations(Transc_ProtOnly, SwissProt_ProtOnly, InterProScan_ProtOnly)
print("ProtOnly: databases matches results")
for ann in ProtOnly_res.keys():
	print(ann, ProtOnly_res[ann])
print("Total transcripts/genes to keep:", len(FilteredTrans_ProtOnly))

with open(output_path+".GenesToKeep.txt", "w") as output_file:
	for transc in FilteredTrans_ProtOnly:
		print(transc, file=output_file)




