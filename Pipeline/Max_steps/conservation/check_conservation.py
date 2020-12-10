





def parse_aln(file):
	my_dict = {}
	for line in open(file).readlines()[1:]:
		line = line.strip().split()
		if not len(line) == 2:
			continue
		gene = line[0]
		seq = line[1]
		if "*" in seq:
			continue
		if not my_dict.has_key(gene):
			my_dict[gene] = ""
		my_dict[gene] += seq

	return my_dict

def min_diff(a,b):
	min_diffs = 999
	for i in a:
		for j in b:
			diff = abs(i-j)
			if diff < min_diffs:
				min_diffs = diff
	return min_diffs
def parse_fimo(file):
	my_dict = {}
	for line in open(file).readlines()[1:]:
		line = line.strip().split()
		motif = line[0]
		gene = line[1].split("|")[0]
		start = int(line[2])
		if not my_dict.has_key(motif):
			my_dict[motif] = {}
		if not my_dict[motif].has_key(gene):
			my_dict[motif][gene] = []
		my_dict[motif][gene].append(start)
	return my_dict
	
def adjust_pos(seq,pos):
	adjust_pos = pos
	count_original_pos = 0
	for i in range(len(seq)):
		if seq[i] == "-":
			adjust_pos += 1
		else:
			count_original_pos += 1
		if count_original_pos == pos:
			return adjust_pos

import glob
gene_list = [] 
for item in glob.glob("AT*.aln"):
	item = item.split(".")[0]
	gene_list.append(item)
import sys
T_fimo = 'Ath.fimo'
L_fimo = "lyrata.fimo"

T_fimo_pos = parse_fimo(T_fimo)
L_fimo_pos = parse_fimo(L_fimo)
motif_cons = {}		
dis_cut = 1
for motif in T_fimo_pos:
	# print motif
	motif_cons[motif]={}
	motif_cons[motif]["cons"] = []
	motif_cons[motif]["amg"] = []
	if not L_fimo_pos.has_key(motif):
		# print "not in"
		continue
	for gene in T_fimo_pos[motif]:
		gene_temp = gene.split("|")[0]
		if gene_temp in gene_list:
			# print gene
			aln = parse_aln(gene_temp+".orthologs.aln")
			
			orthologs = aln.keys()
			T_pos_array = T_fimo_pos[motif][gene]
			T_adj_pos_array = []
			L_adj_pos_array = []
			for T_pos in T_pos_array:
				T_adj_pos = adjust_pos(aln[gene_temp],T_pos)
				T_adj_pos_array.append(T_adj_pos)
			# print T_adj_pos_array
			orthologs.remove(gene_temp)
			ortholog1 = orthologs[0]
			
			# ortholog2 = orthologs[1]
			# print 
			if L_fimo_pos[motif].has_key(orthologs[0]):
				# print "asd"
				ortholog1 = orthologs[0]
			# if L_fimo_pos[motif].has_key(orthologs[1]):
				# ortholog1 = orthologs[1]
			# print ortholog1
			
			if L_fimo_pos[motif].has_key(ortholog1):
				# print ortholog1
				motif_cons[motif]["amg"].append("C_elegans,"+ortholog1)
				C_pos_array = L_fimo_pos[motif][ortholog1]
				for C_pos in C_pos_array:
					C_adj_pos = adjust_pos(aln[ortholog1],C_pos)
					L_adj_pos_array.append(C_adj_pos)
				motif_cons[motif]["amg"].append("thaliana: " + ",".join(map(lambda x:str(x),T_adj_pos_array)) + "; lyrata: " + ",".join(map(lambda x:str(x),L_adj_pos_array)))
			if min_diff(T_adj_pos_array,L_adj_pos_array) < dis_cut:
				motif_cons[motif]["cons"].append("thaliana " + gene_temp +" : " + ",".join(map(lambda x:str(x),T_adj_pos_array)) + "@" + ",".join(map(lambda x:str(x),T_pos_array)) + "; lyrata "+ortholog1+" : " + ",".join(map(lambda x:str(x),L_adj_pos_array)))

out = open("motif_conservation_analysis.tsv","wb")				
for motif in motif_cons:
	if not len(motif_cons[motif]["cons"]) == 0:
		out_line = [motif] + ["conserved"] + motif_cons[motif]["cons"]
		print >>out,"\t".join(out_line)
	else:
		continue
		# if not len(motif_cons[motif]["amg"]) == 0:
			# out_line = [motif] + ["ambigous"] + motif_cons[motif]["amg"]
			# print >>out,"\t".join(out_line)