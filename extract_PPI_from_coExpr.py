import sys

f = open(sys.argv[1])
inweb = open("/Users/jessicaxinhu/Documents/noa_project/InWeb/InWeb_IM/inwebIM_noCutoff.tsv")


interact = set()
for line in inweb:
	tmp = line.strip().split("\t")
	interact.add((tmp[0], tmp[1]))


for line in f:
	tmp = line.strip().split("\t")
	gene_pair = (tmp[0], tmp[1])
	gene_pair_swap = (tmp[1], tmp[0])

	if gene_pair in interact or gene_pair_swap in interact:
		print line.strip() 

	



