import sys 

blastp_input = sys.argv[1]
blastp_out = sys.argv[2]

with open(blastp_out, "w") as out_file: 
	with open(blastp_input, "r") as input_file:

		Analyzing_gene = ""

		count_no_XP = 0
		count_Annotated_XP = 0
		total = 0

		for line in input_file:
			line = line.rstrip()
			line_list = line.split("\t")
			temp_analyzing_gene = line_list[0]
			annotation = line_list[1]

			#First line
			if Analyzing_gene == "" :
				Analyzing_gene = temp_analyzing_gene
				temp_list_annotation = [line_list]
				continue

			#Following lines
			if Analyzing_gene == temp_analyzing_gene:
				temp_list_annotation.append(line_list)

			elif Analyzing_gene != temp_analyzing_gene:
				# save the annotation for previous gene
				#print(len(temp_list_annotation))
				foundAnnotation = False
				for element in temp_list_annotation:
					if element[1].startswith("XP_"):
						foundAnnotation = True
						right_annotation = "\t".join(element)
						print(right_annotation, file=out_file)
						count_Annotated_XP += 1
						break
				if not foundAnnotation:
					count_no_XP += 1
					print("\t".join(temp_list_annotation[0]), file=out_file)

				total += 1
				#Replace info for new gene
				Analyzing_gene = temp_analyzing_gene
				temp_list_annotation = [line_list]


		#Last gene
		foundAnnotation = False
		for element in temp_list_annotation:
			if element[1].startswith("XP_"):
				foundAnnotation = True
				right_annotation = "\t".join(element)
				count_Annotated_XP += 1
				print(right_annotation, file=out_file)
				break
		if not foundAnnotation:
			count_no_XP += 1
			print("\t".join(temp_list_annotation[0]), file=out_file)
		total += 1


print("Total annotated transcripts: ", total)
print("Annotated with XP_: ", count_Annotated_XP)
print("Annotated without XP_: ", count_no_XP)




