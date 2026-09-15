import sys

infile_name = sys.argv[1]
outfile_name = sys.argv[2]

# To be conservative, we only consider:
# SYNAL
# INVAL
# INVTRAL
# TRANSAL
consid_align = ["SYNAL","INVAL", "INVTRAL","TRANSAL"]

with open(outfile_name, "w") as outfile: 
	with open(infile_name, "r") as infile:
		for line in infile:
			line = line.rstrip()
			line_list = line.split("\t")
			if line_list[10] in consid_align:
				print("\t".join(line_list[:3]), file=outfile)
			




