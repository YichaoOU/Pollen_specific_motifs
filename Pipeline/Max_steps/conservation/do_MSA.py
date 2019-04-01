import os
def parse_list(file):
	temp = {}
	for line in open(file).readlines():
		temp[line.strip()] = ""
	return temp

def parse_bed_file(file):
	my_dict = {}
	flag = True
	with open(file) as f:
		for line in f:
			if flag:
				flag = False
				continue
			line = line.strip().split()
			seq_id = line[0]
			gene_id = line[5]
			my_dict[gene_id] = seq_id

	return my_dict	
	
def parse_orthologs(file):
	for line in open(file).readlines():
		line = line.strip().split(",")
		T = line[0]
		L = line[3]
		do_MSA_2_seqs(T,L)

from Bio import SeqIO as io
from Bio.Seq import Seq
def read_fasta(file):
	
	seq_hash = {}
	f=open(file,"rU")
	record = io.parse(f,"fasta")
	
	for r in record:
		seq_hash[str(r.id).split("|")[0]]=str(r.seq)
	return seq_hash		

T_seq = read_fasta("/home/working/ISMB_proceeding/HRGP/hrgp_fore.fasta")
L_seq = read_fasta("/home/working/HRGP/lyrata.fasta.txt")
def do_MSA_2_seqs(id1,id2):
	try:
		seq1 = T_seq[id1]
		seq2 = L_seq[id2]
	except:
		print "something wrong"
		return 1
	out = open(id1+".orthologs.fa","wb")
	print >>out,">"+id1
	print >>out,seq1
	print >>out,">"+id2
	print >>out,seq2
	out.close()
	command = "clustalw2 -gapopen=10 -gapext=0.1 -infile=" + id1+".orthologs.fa"
	os.system(command)
	
import sys		
parse_orthologs("ath-homologs.csv")



