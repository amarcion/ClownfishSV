from Bio import SeqIO
import os
import sys

###

def countLower(Sequence):
	count = 0
	for nuc in Sequence:
		if nuc == "a" or nuc == "t" or nuc == "c" or nuc == "g":
			count += 1
	return(count)

###

wd = sys.argv[1]
input_file_prefix = sys.argv[2]
output_file_prefix = sys.argv[3]


os.chdir(wd)

Dictionary_genes_transcripts = {}
# gene: [t1, t2, etc]

Dictionary_transcripts_info = {}
# transcript : [length, nb_softmasked, seq]

Largest_transcr = {}
# gene : [record.id, length]

transcript_no_softmask = 0
min_length = 10000000
max_length  = 0

for record in SeqIO.parse(input_file_prefix+".codingseq", "fasta"):
	gene_name_index = record.id.find(".")
	gene_name = record.id[:gene_name_index]
	transc_nb = record.id[gene_name_index+1:]
	
	if gene_name not in Dictionary_genes_transcripts.keys():
		Dictionary_genes_transcripts[gene_name] = []
	Dictionary_genes_transcripts[gene_name].append(transc_nb)

	tr_length = len(record.seq)
	if tr_length <= min_length:
		min_length = tr_length
	if tr_length >= max_length:
		max_length = tr_length

	tr_lowcase = countLower(record.seq)
	if tr_lowcase == 0:
		transcript_no_softmask += 1
	Dictionary_transcripts_info[record.id] = [tr_length, tr_lowcase, record.seq]

	if gene_name not in Largest_transcr.keys():
		Largest_transcr[gene_name] = [gene_name+"."+transc_nb, tr_length]
	else:
		if tr_length > Largest_transcr[gene_name][1]:
			Largest_transcr[gene_name] = [gene_name+"."+transc_nb, tr_length]


# Save coding Sequences of longest transcripts
with open(output_file_prefix + ".LongestIso.codingseq", "w") as out_file_seq:
	longest_trans = []
	nb_longestIso_no_softmasked = 0
	min_longestIso =1000000
	max_longestIso = 0

	for gene in Largest_transcr.keys():
		transc = Largest_transcr[gene][0]
		longest_trans.append(transc)
		print(">"+transc, file=out_file_seq)
		print(str(Dictionary_transcripts_info[transc][-1]), file=out_file_seq)

		if Dictionary_transcripts_info[transc][1] == 0:
			nb_longestIso_no_softmasked += 1
		if Largest_transcr[gene][1] <= min_longestIso:
			min_longestIso = Largest_transcr[gene][1]
		if Largest_transcr[gene][1] >= max_longestIso:
			max_longestIso = Largest_transcr[gene][1]

# Save aa sequences of longest transcripts
with open(output_file_prefix + ".LongestIso.aa", "w") as out_file_aa:
	for record in SeqIO.parse(input_file_prefix+".aa", "fasta"):
		if record.id in longest_trans:
			print(">"+record.id, file=out_file_aa)
			print(record.seq, file=out_file_aa)


# Save summary table
with open(output_file_prefix + ".summary", "w") as out_file_summary:
	what_to_print = "Transcript\tGene\tLength\tNbSoftMaskedCar\tProportionSoftMasked\tLongestIsoform"
	print(what_to_print, file=out_file_summary)
	for transcript in Dictionary_transcripts_info.keys():
		indx = transcript.find(".")
		gene = transcript[:indx]
		length = Dictionary_transcripts_info[transcript][0]
		SoftMaskedCar = Dictionary_transcripts_info[transcript][1]
		Prop_soft = round(float(SoftMaskedCar)/float(length), 2)
		if transcript == Largest_transcr[gene][0]:
			longIso = "Yes"
		else:
			longIso = "No"

		what_to_print = transcript + "\t" + gene + "\t" + str(length) + "\t" + str(SoftMaskedCar) + "\t" + str(Prop_soft) + "\t" + longIso
		print(what_to_print, file=out_file_summary)


print("Number of transcripts: ", len(Dictionary_transcripts_info.keys()))
print("Number transcripts without softmasking: ", transcript_no_softmask)
print("Minimum transcript length: ", min_length)
print("Maximum transcript length: ", max_length)
print("Number of genes: ", len(Dictionary_genes_transcripts.keys()))
print("Number LongestIsoform without softmasking: ", nb_longestIso_no_softmasked)
print("Minimum length of Longest isoform: ", min_longestIso)
print("Maximum length of Longest isoform: ", max_longestIso)



