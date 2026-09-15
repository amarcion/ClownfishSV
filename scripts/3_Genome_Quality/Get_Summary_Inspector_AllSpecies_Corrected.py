import sys
import os

path_to_files = sys.argv[1]
out_file_name = sys.argv[2]


prefix_assemblies = {"Hifiasm": "a_Corr_", "DupPurged": "b_Corr_", "RepScaffNoMito": "c_Corr_", "Samba": "d_Corr_", "SambaNtLink": "e_Corr_", "Final": "f_Corr_"}

def get_statistics_Inspector(input_file_name):
	List_of_stats = []
	List_of_values = []

	with open(input_file_name, "r") as input_file:
		for line in input_file:
			line = line.rstrip()
			line_list = line.split("\t")

			if len(line_list) == 2:
				List_of_stats.append(line_list[0])
				List_of_values.append(line_list[1])
	
	return(List_of_stats, List_of_values)


list_of_dir = os.listdir(path_to_files)
count = 0

List_of_results = []

for el in list_of_dir:
	if not os.path.isdir(path_to_files+"/"+el):
		continue
	if not el.startswith("OutCorr_"):
		continue

	Species = el.split("_")[1]
	Assembly = el.split("_")[2]
	
	tmp_file_name = path_to_files + "/" + el + "/summary_statistics"
	if not os.path.isfile(tmp_file_name):
		print("File " + tmp_file_name + " for species " + Species + ", Assembly " + Assembly + " not found! ")
		continue

	Stats_order, temp_values = get_statistics_Inspector(tmp_file_name)
	Values_Sp_Ass = [Species, prefix_assemblies[Assembly]+Assembly] + temp_values
	List_of_results.append(Values_Sp_Ass)


# Format The results and Save

Col_names = ["Species", "Assembly" ] + Stats_order
with open(path_to_files+"/"+out_file_name, "w") as out_file:
	print("\t".join(Col_names), file=out_file)

	List_of_results.sort()
	for el in List_of_results:
		new_res = "\t".join(el)
		print(new_res, file=out_file)





