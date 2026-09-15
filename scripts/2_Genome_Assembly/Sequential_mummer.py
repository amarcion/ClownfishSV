import sys
import os
from Bio import SeqIO
from statistics import mean

input_file_sequences = sys.argv[1]
temp_directory = sys.argv[2]
out_prefix = sys.argv[3]
min_coverge_align = float(sys.argv[4])
min_Meanidentity = float(sys.argv[5])


def get_dictionary_sequences(input_file_fasta):
	Dict_seq = {}
	Length_sequences = []

	for record in SeqIO.parse(input_file_fasta, "fasta"):
		Dict_seq[record.id] = record.seq
		Length_sequences.append((len(record.seq), record.id))
	return(Dict_seq, Length_sequences)

def read_delta_file(input_file_delta):
	Query_scaffolds = {}
	Query_scaffolds_identity_length = {}
	input_file = open(input_file_delta, "r")
	for line in input_file:
		line = line.rstrip()
		line_list = line.split("\t")
		scaff = line_list[10]

		# Get alignment statistics
		al_start = min(int(line_list[2]), int(line_list[3]))
		al_stop = max(int(line_list[2]), int(line_list[3]))
		al_length = int(line_list[5])
		scaff_length = int(line_list[8])
		identity = float(line_list[6])

		if scaff not in Query_scaffolds_identity_length.keys():
			Query_scaffolds_identity_length[scaff] = ([identity], scaff_length)
		else:
			Query_scaffolds_identity_length[scaff][0].append(identity)

		if scaff not in Query_scaffolds.keys():
			Query_scaffolds[scaff] = ["0"]*scaff_length
			Query_scaffolds[scaff][al_start:al_stop+1] = ["1"]*al_length
		else:
			Query_scaffolds[scaff][al_start:al_stop + 1] = ["1"] * al_length

	input_file.close()
	return(Query_scaffolds, Query_scaffolds_identity_length)


def get_statistics_Alignmets(Query_scaffolds, Query_scaffolds_identity_length):
	if len(Query_scaffolds.keys()) != 1:
		print("Problem: more than one scaffold!")
	for scaffold in Query_scaffolds.keys():
		Nb_Aligned_bases = Query_scaffolds[scaffold].count("1")
		Perc_bases = Nb_Aligned_bases/Query_scaffolds_identity_length[scaffold][1]
		ident = mean(Query_scaffolds_identity_length[scaffold][0])

	return(Perc_bases, ident)

def is_non_zero_file(fpath):  
    return os.path.isfile(fpath) and os.path.getsize(fpath) > 0


Scaffols_dict, Length_scaffolds = get_dictionary_sequences(input_file_sequences)
Length_scaffolds_Sort = sorted(Length_scaffolds, reverse = True)

# we keep the longest scaffolds in any case and create a file with it
if not os.path.isfile(temp_directory+"reference.tmp.fa"):
	ref_file = open(temp_directory+"reference.tmp.fa", "w")
	long_scaff = Length_scaffolds_Sort[0][1]
	print("Long_scaff: ", long_scaff)
	print(">"+long_scaff, file=ref_file)
	print(Scaffols_dict[long_scaff], file=ref_file)
	ref_file.close()

if not os.path.isfile(temp_directory+"AllAlignments.Coord.tsv"):
	delta_File = open(temp_directory+"AllAlignments.Coord.tsv", "w")
	delta_File.close()

if not os.path.isfile(temp_directory+"AllAlignmentStatistics.tsv"):
	AlignStat_File = open(temp_directory+"AllAlignmentStatistics.tsv", "w")
	print("Scaffold", "Length", "Proportion_Bases_Aligned", "Mean_Identity",  sep="\t", file=AlignStat_File)
	AlignStat_File.close()

if not os.path.isfile(temp_directory+"RemovedScaffolds.fa"):
	RemovedScaff_File = open(temp_directory+"RemovedScaffolds.fa", "w")
	RemovedScaff_File.close()


for element in Length_scaffolds_Sort[1:]:

	# We align each file in increasing length:
	considered_scaff = element[1]
	print(considered_scaff)

	# Create query file
	query_file = open(temp_directory+"query.tmp.fa", "w")
	print(">"+considered_scaff, file=query_file)
	print(Scaffols_dict[considered_scaff], file=query_file)
	query_file.close()

	# Do the alignments
	nucmer_command = "nucmer -b 500 -g 200 -t 12 -p "+temp_directory+"numer.tmp " + temp_directory+"reference.tmp.fa " + temp_directory+"query.tmp.fa"
	show_coords_command = "show-coords -roTlH "+ temp_directory+"numer.tmp.delta > " + temp_directory+"numer.tmp.Coord.tsv"
	print(nucmer_command)
	os.system(nucmer_command)
	print(show_coords_command)
	os.system(show_coords_command)

	# Check if coordinates files is present
	if is_non_zero_file(temp_directory+"numer.tmp.Coord.tsv"):
		AlignScaff_Dic, ScaffIdentity_Dic = read_delta_file(temp_directory+"numer.tmp.Coord.tsv")
		BasesCovered, IdentityAlign = get_statistics_Alignmets(AlignScaff_Dic, ScaffIdentity_Dic)

		# We append the information on the scaffold
		AlignStat_File = open(temp_directory+"AllAlignmentStatistics.tsv", "a")
		what_to_print=considered_scaff+"\t"+str(len(Scaffols_dict[considered_scaff])) + "\t" + str(BasesCovered) + "\t" + str(IdentityAlign) 
		print(what_to_print, file=AlignStat_File)
		AlignStat_File.close()

		if BasesCovered <= min_coverge_align or IdentityAlign <= min_Meanidentity:
			#we add the file to the reference as we keep it as not already represented
			print(considered_scaff + " is kept")

			# We write the file
			ref_file = open(temp_directory+"/reference.tmp.fa", "a")
			print(">"+considered_scaff, file=ref_file)
			print(Scaffols_dict[considered_scaff], file=ref_file)
			ref_file.close()

			# We concatenate the delta files
			command_cat = "cat " + temp_directory+"AllAlignments.Coord.tsv " + temp_directory+"numer.tmp.Coord.tsv > " + temp_directory+"AllAlignments.Coord_2.tsv"
			command_mv = "mv " + temp_directory+"AllAlignments.Coord_2.tsv " + temp_directory+"AllAlignments.Coord.tsv"
			print(command_cat)
			os.system(command_cat)
			print(command_mv)
			os.system(command_mv)

			# we remove the files we do not need for prepare for the next round
			command_rm = "rm " +  temp_directory+"numer.tmp.Coord.tsv " + temp_directory+"numer.tmp.delta " + temp_directory+"query.tmp.fa"
			print(command_rm)
			os.system(command_rm)

		else:
			print(considered_scaff + " is removed")

			# We print on the removed scaffolds
			RemovedScaff_File = open(temp_directory+"RemovedScaffolds.fa", "a")
			print(">"+considered_scaff, file=RemovedScaff_File)
			print(Scaffols_dict[considered_scaff], file=RemovedScaff_File)

			RemovedScaff_File.close()


	# If not coord file exists, it means that the scaffold did not align and thus we append it to the reference		
	else:
		print(considered_scaff + " did not align and it is kept")

		# We write the file
		ref_file = open(temp_directory+"/reference.tmp.fa", "a")
		print(">"+considered_scaff, file=ref_file)
		print(Scaffols_dict[considered_scaff], file=ref_file)
		ref_file.close()

		# We write in the align statistics (not aligned) in the AllAlignmentStatistics.tsv
		AlignStat_File = open(temp_directory+"AllAlignmentStatistics.tsv", "a")
		what_to_print=considered_scaff+"\t"+str(len(Scaffols_dict[considered_scaff])) + "\t0\t0" 
		print(what_to_print, file=AlignStat_File)
		AlignStat_File.close()

# We take the references at the end as it is the onces we removed
final_command_cp = "cp " + temp_directory+"reference.tmp.fa " + out_prefix + ".RepresentativeScaffolds.fa"
print(final_command_cp)
os.system(final_command_cp)
final_command_cp2 = "cp " + temp_directory+"AllAlignmentStatistics.tsv " + out_prefix + ".AllAlignmentStatistics.tsv"
print(final_command_cp2)
os.system(final_command_cp2)
final_command_cp3 = "cp " + temp_directory+"RemovedScaffolds.fa " + out_prefix + ".RemovedScaffolds.fa"
print(final_command_cp3)
os.system(final_command_cp3)
final_command_cp4 = "cp " + temp_directory+"AllAlignments.Coord.tsv " + out_prefix + ".AllAlignments.Coord.tsv"
print(final_command_cp4)
os.system(final_command_cp4)






