import sys
from Bio import SeqIO

assembly_file_name = sys.argv[1]
assembly_out_name = sys.argv[2]
splitting_point_file = sys.argv[3]

#scaffolds = sys.argv[3].split(",")
#coords = sys.argv[4].split(",")

# Retrieve Splitting Points
scaffolds_coords = {}

with open(splitting_point_file, "r") as infile:
	for line in infile:
		line = line.rstrip()
		line_list = line.split("\t")
		if line_list[0] not in scaffolds_coords.keys():
			scaffolds_coords[line_list[0]] = "-".join(line_list[1:3]) 
		else:
			scaffolds_coords[line_list[0]] += "-" + "-".join(line_list[1:3]) 

scaffolds = []
coords = []
for scaff in scaffolds_coords.keys():
	scaffolds.append(scaff)
	coords.append(scaffolds_coords[scaff])
# print(scaffolds)
# print(coords)


with open(assembly_out_name, "w") as out_file:
	for record in SeqIO.parse(assembly_file_name, "fasta"):
		if record.id not in scaffolds:
			print(">"+record.id[:-2], file=out_file)
			print(record.seq, file=out_file)

		else:
			for i in range(len(scaffolds)):
				if scaffolds[i] == record.id:
					if len(coords[i].split("-")) > 2:
						print("Double split ", scaffolds[i])
						subscaff_1 = scaffolds[i][:-2]+"a"
						subscaff_2 = scaffolds[i][:-2]+"b"
						subscaff_3 = scaffolds[i][:-2]+"c"

						start1, stop1, start2, stop2 = coords[i].split("-")
						seq1 = record.seq[:int(start1)+1]
						seq2 = record.seq[int(stop1):int(start2)+1]
						seq3 = record.seq[int(stop2):]

						print(subscaff_1, subscaff_2, subscaff_3)
						print(">"+subscaff_1, file=out_file)
						print(seq1, file=out_file)
						print(">"+subscaff_2, file=out_file)
						print(seq2, file=out_file)	
						print(">"+subscaff_3, file=out_file)
						print(seq3, file=out_file)
					else:
						print("Single cut ", scaffolds[i])
						start, stop = coords[i].split("-")
						subscaff1 = scaffolds[i][:-2]+"a"
						subscaff2 = scaffolds[i][:-2]+"b"	
						seq1 = record.seq[:int(start)+1]
						seq2 = record.seq[int(stop):]
						print(subscaff1, subscaff2)
						print(">"+subscaff1, file=out_file)
						print(seq1, file=out_file)
						print(">"+subscaff2, file=out_file)
						print(seq2, file=out_file)				
					
