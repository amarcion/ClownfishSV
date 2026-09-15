from Bio import SeqIO
import sys
import os

species = sys.argv[1]
outfile_name = sys.argv[2]
input_file_suffix = sys.argv[3]

os.chdir(species)


with open(outfile_name, "w") as outfile:

	# Not scaffolded
	notScaffolded = species+".ScaffoldsNotToScaffold.fa"
	for record in SeqIO.parse(notScaffolded, "fasta"):
		print(">"+record.id, file=outfile)
		print(record.seq, file=outfile)


	#Scaffolded
	folders = [f for f in os.listdir() if os.path.isdir(f)]
	for each_fold in folders:
		file_to_open = each_fold + "/"+species+".temp."+each_fold+input_file_suffix
		for record in SeqIO.parse(file_to_open, "fasta"):

			if "ntLin" in record.id:
				scaffs = record.description.split(" ")[-1].split(",")
				new_id = ""
				for el in scaffs:
					if ":" in el:
						el_level2 = el.split(":")
						for subel in el_level2:
							if "ptg" in subel:
								new_id += subel+"_"
					else:
						if "ptg" in el:
							new_id += el[:-1]+"_"

				new_id += "ntLink"
				print(new_id)
				print(">"+new_id, file=outfile)
				print(record.seq, file=outfile)					

			elif ":" in record.id:
				scaffs = record.id.split(":")
				new_id = ""
				for el in scaffs:
					if "ptg" in el:
						new_id += el+"_"
				print(new_id[:-1])
				print(">"+new_id[:-1], file=outfile)
				print(record.seq, file=outfile)		

			else:
				print(">"+record.id, file=outfile)
				print(record.seq, file=outfile)

