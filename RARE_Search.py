import os
def rev_comp(sequence):
	new_Seq = ""
	for base in reversed(sequence):
		if base == "A":
			new_Seq += "T"
		if base == "C":
			new_Seq += "G"
		if base == "G":
			new_Seq += "C"
		if base == "T":
			new_Seq += "A"
			
	return new_Seq


'''
bases = {"A":["A"], "C":["C"],"G":["G"], "T":["T"], "R":["A","G"], "Y":["C","T"], "S":["G","C"], "W":["A","T"], "K":["G","T"], "M":["A","C"], "B":["C","G","T"], "D":["A","G","T"], "H":["A","C","T"], "V":["A","C","G"], "N":["A","C","G","T"]}
'''

permutations = ['AGGTCA', 'AGGTGA', 'AGTTCA', 'AGTTGA', 'GGGTCA', 'GGGTGA', 'GGTTCA', 'GGTTGA']
rev_permutations = []

for p in permutations:
	rev_permutations.append(rev_comp(p))
	


gene_list = os.listdir("C:\Users\Thomas\Desktop\Dickinson_Lab\RARA_Promoter_Sequences\human_mouse")

headers = []
sequences = []

for file in gene_list:
	fasta = open("C:\Users\Thomas\Desktop\Dickinson_Lab\RARA_Promoter_Sequences\human_mouse\\" + file, "r")
	gene = fasta.readlines()
	fasta.close()
	
	gene_seq = ""
	
	for i in range(1,len(gene)):
		gene_seq += gene[i].strip().upper()

	headers.append(gene[0].strip())
	sequences.append(gene_seq)
	


RARA_Sites = {}
RARA_Sequences = {}



for i in range(0,len(sequences)):
	
	RARA_Sequences[headers[i]] = []
	RARA_Sites[headers[i]] = []

	if "+++" in headers[i]:
		upstream_TSS = sequences[i][0000:20000]
		downstream_TSS = sequences[i][20000:40000]
		upstream_End = sequences[i][-40000:-20000]
		downstream_End = sequences[i][-20000:-1]
		
	if "---" in headers[i]:
		upstream_TSS = rev_comp(sequences[i][-20000:-1])
		downstream_TSS = rev_comp(sequences[i][-40000:-20000])
		upstream_End = rev_comp(sequences[i][20000:40000])
		downstream_End = rev_comp(sequences[i][0000:20000])

		
	for j in range(0,len(upstream_TSS)-17):
		if upstream_TSS[j:j+2] == "AG" or upstream_TSS[j:j+2] == "GG":
			if upstream_TSS[j:j+6] in permutations:
				if upstream_TSS[j+11:j+17] in permutations:
					if "+++" in headers[i]:
						RARA_Sites[headers[i]].append(str(20000-j) + "..upstream_TSS")
						RARA_Sequences[headers[i]].append(upstream_TSS[j:j+17])
					if "---" in headers[i]:
						RARA_Sites[headers[i]].append(str(20000-j) + "..upstream_TSS")
						RARA_Sequences[headers[i]].append(upstream_TSS[j:j+17])
						
		if upstream_TSS[j:j+2] == "TG" or upstream_TSS[j:j+2] == "TC":
			if upstream_TSS[j:j+6] in rev_permutations:
				if upstream_TSS[j+11:j+17] in rev_permutations:
					if "+++" in headers[i]:
						RARA_Sites[headers[i]].append(str(20000-j) + "..upstream_TSS")
						RARA_Sequences[headers[i]].append(upstream_TSS[j:j+17])
					if "---" in headers[i]:
						RARA_Sites[headers[i]].append(str(20000-j) + "..upstream_TSS")
						RARA_Sequences[headers[i]].append(upstream_TSS[j:j+17])
						
	for j in range(0,len(downstream_TSS)-17):
		if downstream_TSS[j:j+2] == "AG" or downstream_TSS[j:j+2] == "GG":
			if downstream_TSS[j:j+6] in permutations:
				if downstream_TSS[j+11:j+17] in permutations:
					if "+++" in headers[i]:
						RARA_Sites[headers[i]].append(str(j) + "..downstream_TSS")
						RARA_Sequences[headers[i]].append(downstream_TSS[j:j+17])
					if "---" in headers[i]:
						RARA_Sites[headers[i]].append(str(j) + "..downstream_TSS")
						RARA_Sequences[headers[i]].append(downstream_TSS[j:j+17])
						
		if downstream_TSS[j:j+2] == "TC" or downstream_TSS[j:j+2] == "TG":
			if downstream_TSS[j:j+6] in rev_permutations:
				if downstream_TSS[j+11:j+17] in rev_permutations:
					if "+++" in headers[i]:
						RARA_Sites[headers[i]].append(str(j) + "..downstream_TSS")
						RARA_Sequences[headers[i]].append(downstream_TSS[j:j+17])
					if "---" in headers[i]:
						RARA_Sites[headers[i]].append(str(j) + "..downstream_TSS")
						RARA_Sequences[headers[i]].append(downstream_TSS[j:j+17])
						
	for j in range(0,len(upstream_End)-17):
		if upstream_End[j:j+2] == "AG" or upstream_End[j:j+2] == "GG":
			if upstream_End[j:j+6] in permutations:
				if upstream_End[j+11:j+17] in permutations:
					if "+++" in headers[i]:
						RARA_Sites[headers[i]].append(str(20000-j) + "..upstream_End")
						RARA_Sequences[headers[i]].append(upstream_End[j:j+17])
					if "---" in headers[i]:
						RARA_Sites[headers[i]].append(str(20000-j) + "..upstream_End")
						RARA_Sequences[headers[i]].append(upstream_End[j:j+17])
						
		if upstream_End[j:j+2] == "TC" or upstream_End[j:j+2] == "TG":
			if upstream_End[j:j+6] in rev_permutations:
				if upstream_End[j+11:j+17] in rev_permutations:
					if "+++" in headers[i]:
						RARA_Sites[headers[i]].append(str(20000-j) + "..upstream_End")
						RARA_Sequences[headers[i]].append(upstream_End[j:j+17])
					if "---" in headers[i]:
						RARA_Sites[headers[i]].append(str(20000-j) + "..upstream_End")
						RARA_Sequences[headers[i]].append(upstream_End[j:j+17])
						
	for j in range(0,len(downstream_End)-17):
		if downstream_End[j:j+2] == "AG" or downstream_End[j:j+2] == "GG":
			if downstream_End[j:j+6] in permutations:
				if downstream_End[j+11:j+17] in permutations:
					if "+++" in headers[i]:
						RARA_Sites[headers[i]].append(str(j) + "..downstream_End")
						RARA_Sequences[headers[i]].append(downstream_End[j:j+17])
					if "---" in headers[i]:
						RARA_Sites[headers[i]].append(str(str(j) + "..downstream_End"))
						RARA_Sequences[headers[i]].append(downstream_End[j:j+17])
						
		if downstream_End[j:j+2] == "TC" or downstream_End[j:j+2] == "TG":
			if downstream_End[j:j+6] in rev_permutations:
				if downstream_End[j+11:j+17] in rev_permutations:
					if "+++" in headers[i]:
						RARA_Sites[headers[i]].append(str(j) + "..downstream_End")
						RARA_Sequences[headers[i]].append(downstream_End[j:j+17])
					if "---" in headers[i]:
						RARA_Sites[headers[i]].append(str(j) + "..downstream_End")
						RARA_Sequences[headers[i]].append(downstream_End[j:j+17])
						
						
						


output = open("C:\Users\Thomas\Desktop\Dickinson_Lab\RARA_Promoter_Sequences\RARA_Output_hm.txt", "w")
output.writelines("Coordinates given are distances from the transcription start site.\n\n")
	
for head in headers:
	output.writelines(head + "\n" + "\n")
	for i in range(0,len(RARA_Sequences[head])):	
		output.writelines(RARA_Sites[head][i] + "\t")
		output.writelines(RARA_Sequences[head][i] + "\n" + "\n")
		
output.close()
				
				
