import os
import itertools


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


rare = raw_input("What is your RARE? Use IUPAC notation for positions where multiple bases are allowed. ")

sep = int(raw_input("How many bases separate your direct repeats? "))

search_len = int(raw_input("How large is your search area? ")

stream = raw_input("Are you searching upstream (up), downstream (down), or both (both)? ")

rare_split = []
for base in rare:
	rare_split.append(bases[base])
	
permutations = list(itertools.product(*rare_split))
'''
	
permutations = ['AGGTCA', 'AGGTGA', 'AGTTCA', 'AGTTGA', 'GGGTCA', 'GGGTGA', 'GGTTCA', 'GGTTGA']
rev_permutations = []

for p in permutations:
	rev_permutations.append(rev_comp(p))
	

gff = open("C:\Users\Thomas\Desktop\Dickinson_Lab\RARA_Promoter_Sequences\XENLA_9.2_Xenbase.gff3")
gene_info = gff.readlines()
gff.close()

del gene_info[0:8]

fasta = open("C:\Users\Thomas\Desktop\Dickinson_Lab\RARA_Promoter_Sequences\XL9_2.fa")
genome = fasta.readlines()
fasta.close()

genome_dict = {}
for i in range(0,len(genome)):
	if genome[i][0] == ">":
		genome_dict[genome[i][1:].strip()] = genome[i+1].strip().upper()

RARA_Sites = {}
RARA_Sequences = {}
gene_names = []


for gene in gene_info:
	gene1 = gene.split("\t")
	gene2 = gene1[-1]
	
	if "Xenbase\tgene" not in gene:
		continue 
	
	if "protein_coding" not in gene2:
		continue 
	
	gene_chr = gene1[0].strip()
	gene_start = int(gene1[3].strip())
	gene_end = int(gene1[4].strip())
	gene_strand = gene1[6].strip()
	gene_name = gene2.split(";")[1].strip()
	
	gene_names.append(gene_name)
	
	if gene_strand == "+":
		upstream_TSS = genome_dict[gene_chr][gene_start - 10000 : gene_start]
		downstream_TSS = genome_dict[gene_chr][gene_start : gene_start + 10000]
		upstream_End = genome_dict[gene_chr][gene_end - 10000 : gene_end]
		downstream_End = genome_dict[gene_chr][gene_end : gene_end + 10000]
		
	if gene_strand == "-":
		upstream_TSS = rev_comp(genome_dict[gene_chr][gene_end : gene_end + 10000])
		downstream_TSS = rev_comp(genome_dict[gene_chr][gene_end - 10000 : gene_end])
		upstream_End = rev_comp(genome_dict[gene_chr][gene_start : gene_start + 10000])
		downstream_End = rev_comp(genome_dict[gene_chr][gene_start - 10000 : gene_start])
	

	RARA_Sites[gene_name] = []
	RARA_Sequences[gene_name] = []
	
	for j in range(0,len(upstream_TSS)-17):
		if upstream_TSS[j:j+2] == "AG" or upstream_TSS[j:j+2] == "GG":
			if upstream_TSS[j:j+6] in permutations:
				if upstream_TSS[j+11:j+17] in permutations:
					if gene_strand == "+":
						RARA_Sites[gene_name].append(gene_chr + str((gene_start-10000)+j) + ".." + str((gene_start-10000)+j+17))
						RARA_Sequences[gene_name].append(upstream_TSS[j:j+17])
					if gene_strand == "-":
						RARA_Sites[gene_name].append(gene_chr + str((gene_end+10000)-j-17) + ".." + str((gene_end+10000)-j))
						RARA_Sequences[gene_name].append(upstream_TSS[j:j+17])
						
		if upstream_TSS[j:j+2] == "TG" or upstream_TSS[j:j+2] == "TC":
			if upstream_TSS[j:j+6] in rev_permutations:
				if upstream_TSS[j+11:j+17] in rev_permutations:
					if gene_strand == "+":
						RARA_Sites[gene_name].append(gene_chr + str((gene_start-10000)+j) + ".." + str((gene_start-10000)+j+17))
						RARA_Sequences[gene_name].append(upstream_TSS[j:j+17])
					if gene_strand == "-":
						RARA_Sites[gene_name].append(gene_chr + str((gene_end+10000)-j-17) + ".." + str((gene_end+10000)-j))
						RARA_Sequences[gene_name].append(upstream_TSS[j:j+17])
						
	for j in range(0,len(downstream_TSS)-17):
		if downstream_TSS[j:j+2] == "AG" or downstream_TSS[j:j+2] == "GG":
			if downstream_TSS[j:j+6] in permutations:
				if downstream_TSS[j+11:j+17] in permutations:
					if gene_strand == "+":
						RARA_Sites[gene_name].append(gene_chr + str((gene_start+j)) + ".." + str((gene_start+j+17)))
						RARA_Sequences[gene_name].append(downstream_TSS[j:j+17])
					if gene_strand == "-":
						RARA_Sites[gene_name].append(gene_chr + str(gene_end-j-17) + ".." + str((gene_end-j)))
						RARA_Sequences[gene_name].append(downstream_TSS[j:j+17])
						
		if downstream_TSS[j:j+2] == "TC" or downstream_TSS[j:j+2] == "TG":
			if downstream_TSS[j:j+6] in rev_permutations:
				if downstream_TSS[j+11:j+17] in rev_permutations:
					if gene_strand == "+":
						RARA_Sites[gene_name].append(gene_chr + str((gene_start+j)) + ".." + str((gene_start+j+17)))
						RARA_Sequences[gene_name].append(downstream_TSS[j:j+17])
					if gene_strand == "-":
						RARA_Sites[gene_name].append(gene_chr + str(gene_end-j-17) + ".." + str((gene_end-j)))
						RARA_Sequences[gene_name].append(downstream_TSS[j:j+17])
						
	for j in range(0,len(upstream_End)-17):
		if upstream_End[j:j+2] == "AG" or upstream_End[j:j+2] == "GG":
			if upstream_End[j:j+6] in permutations:
				if upstream_End[j+11:j+17] in permutations:
					if gene_strand == "+":
						RARA_Sites[gene_name].append(gene_chr + str(gene_end-j-17) + ".." + str(gene_end-j))
						RARA_Sequences[gene_name].append(upstream_End[j:j+17])
					if gene_strand == "-":
						RARA_Sites[gene_name].append(gene_chr + str(gene_start+j) + ".." + str(gene_start+j+17))
						RARA_Sequences[gene_name].append(upstream_End[j:j+17])
						
		if upstream_End[j:j+2] == "TC" or upstream_End[j:j+2] == "TG":
			if upstream_End[j:j+6] in rev_permutations:
				if upstream_End[j+11:j+17] in rev_permutations:
					if gene_strand == "+":
						RARA_Sites[gene_name].append(gene_chr + str(gene_end-j-17) + ".." + str(gene_end-j))
						RARA_Sequences[gene_name].append(upstream_End[j:j+17])
					if gene_strand == "-":
						RARA_Sites[gene_name].append(gene_chr + str(gene_start+j) + ".." + str(gene_start+j+17))
						RARA_Sequences[gene_name].append(upstream_End[j:j+17])
						
	for j in range(0,len(downstream_End)-17):
		if downstream_End[j:j+2] == "AG" or downstream_End[j:j+2] == "GG":
			if downstream_End[j:j+6] in permutations:
				if downstream_End[j+11:j+17] in permutations:
					if gene_strand == "+":
						RARA_Sites[gene_name].append(gene_chr + str((gene_end+j)) + ".." + str((gene_end+j+17)))
						RARA_Sequences[gene_name].append(downstream_End[j:j+17])
					if gene_strand == "-":
						RARA_Sites[gene_name].append(gene_chr + str(gene_start-j-17) + ".." + str(gene_start-j))
						RARA_Sequences[gene_name].append(downstream_End[j:j+17])
						
		if downstream_End[j:j+2] == "TC" or downstream_End[j:j+2] == "TG":
			if downstream_End[j:j+6] in rev_permutations:
				if downstream_End[j+11:j+17] in rev_permutations:
					if gene_strand == "+":
						RARA_Sites[gene_name].append(gene_chr + str((gene_end+j)) + ".." + str((gene_end+j+17)))
						RARA_Sequences[gene_name].append(downstream_End[j:j+17])
					if gene_strand == "-":
						RARA_Sites[gene_name].append(gene_chr + str(gene_start-j-17) + ".." + str(gene_start-j))
						RARA_Sequences[gene_name].append(downstream_End[j:j+17])
						
						
output = open("C:\Users\Thomas\Desktop\Dickinson_Lab\RARA_Promoter_Sequences\RARA_Output_Whole_Genome_All.txt", "w")
output.writelines("Coordinates given are position on the chromosome\n\n")

for gene in gene_names: 
	if len(RARA_Sites[gene])>0:
		output.writelines(gene + "\tNumber of sites: " + str(len(RARA_Sites[gene])) + "\t")
		for i in range(0,len(RARA_Sites[gene])-1):
				output.writelines(RARA_Sites[gene][i] + "\t" + RARA_Sequences[gene][i] + "\t")
		output.writelines(RARA_Sites[gene][-1] + "\t" + RARA_Sequences[gene][-1] + "\n")
		
output.close()