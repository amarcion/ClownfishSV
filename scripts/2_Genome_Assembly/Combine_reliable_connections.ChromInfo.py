import sys
import os

path_files = sys.argv[1]
species=sys.argv[2]
out_path = sys.argv[3]
min_thresh =  float(sys.argv[4]) # min coverage threshold
max_thresh =  float(sys.argv[5]) # mac coverage threshold
chromFile =  sys.argv[6]

pe_minLength = {"pe5000": 1000, "pe7000": 2800, "pe9000": 5400, "pe11000": 6600, "pe15000": 9000} # in set to 1/5, 2/5 and 3/5  for the others
pe_maxLength = {"pe5000": 8000, "pe7000": 9000, "pe9000": 11000, "pe11000": 13000, "pe15000": 17000}
#While the ranges are wide, they allow to account for potential assembly errors and the presence of gaps, both of which can affect the observed insert size

def filter_connection(line_connection, peType, min_thresh, max_thresh):
	line_connection = line_connection.rstrip()
	line_list = line_connection.split("\t")
	# check min coverage
	if float(line_list[2]) < min_thresh:
		return 
	# check max coverage
	if float(line_list[2]) > max_thresh:
		return False
	# check min length
	if float(line_list[3]) <  pe_minLength[peType]:
		return False
	# check max length
	if float(line_list[3]) > pe_maxLength[peType]:
		return False

	return line_list

def get_expected_chrom(expectedChrom_file):
	Dictionary_results = {}
	with open(expectedChrom_file, "r") as input_file:
		for line in input_file:
			line = line.rstrip()
			if line.startswith("Species"):
				continue

			line_list = line.split("\t")
			if line_list[0] not in Dictionary_results.keys():
				Dictionary_results[line_list[0]] = {}

			Dictionary_results[line_list[0]][line_list[3]] = line_list[1]

	return(Dictionary_results)


def Same_Chrom(line_connection, Dict_chrom, species):
	line_connection = line_connection.rstrip()
	line_rp = line_connection.replace("/r", "")
	line_list = line_rp.split("\t")
	if line_list[0] in Dict_chrom[species].keys():
		chr_1 = Dict_chrom[species][line_list[0]]
	else:
		chr_1 = "NA1"

	if line_list[1] in Dict_chrom[species].keys():
		chr_2 = Dict_chrom[species][line_list[1]]
	else:
		chr_2 = "NA2"	

	if chr_1 == chr_2:
		return True
	else:
		return False


# Make folder if not exist
if not os.path.exists(out_path):
    os.makedirs(out_path)

# Get potential chrom
Dictionary_chromosomes = get_expected_chrom(chromFile)


outfile = open(out_path+"/reliable.connections", "w")
outfile2 = open(out_path+"/NoExpectedChrom.connections", "w")
retained_connections = []

for pe_type in pe_minLength.keys():
	print(pe_type)
	print()

	input_file_name = path_files + "/"+species+"."+pe_type+".out/reliable.connections"

	with open(input_file_name, "r") as in_file:
		for line in in_file:
			connection = filter_connection(line, pe_type, min_thresh, max_thresh)

			if connection == False:
				continue

			# avoid multiple connection is the same in different files
			if connection[0:2] not in retained_connections:

				# Check if potential chromosomes
				sc = Same_Chrom(line, Dictionary_chromosomes, species)
				if sc:
					retained_connections.append(connection[0:2])
					print("\t".join(connection), file=outfile)
					print("\t".join(connection))
				else:
					print(line, ": Not expected chrom")
					print("\t".join(connection), file=outfile2)

	in_file.close()

outfile.close()
outfile2.close()
