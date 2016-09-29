import sys


f = open(sys.argv[1])
f_sign1 = open(sys.argv[2])
f_sign2 = open(sys.argv[3])
mod = sys.argv[4]
inweb = open("/Users/jessicaxinhu/Documents/noa_project/InWeb/InWeb_IM/inwebIM_noCutoff.tsv")
out = "/Users/jessicaxinhu/Documents/noa_bc_pilot/diffCoExpress_analysis/diffCoEx/beta12/diff_modules_beta12/diffNetw_inweb/"


# inweb
interact = set()
for line in inweb:
	tmp = line.strip().split("\t")
	interact.add((tmp[0], tmp[1]))






# make file

sign1 = []
sign2 = []
for line in f_sign1:
	tmp = line.strip().split(" ")
	sign1.append(tmp[4])
for line in f_sign2:
	tmp = line.strip().split(" ")
	sign2.append(tmp[4])



my_set = set()
my_set_cytoAttrib = set()
for i,line in enumerate(f):
	tmp = line.strip().split(" ")
	gene_pair = (tmp[0], tmp[1])
	gene_pair_swap = (tmp[1], tmp[0])

	if gene_pair in interact or gene_pair_swap in interact:
		if sign1[i] == "-1" and sign2[i] == "1":
			my_set.add((tmp[0], tmp[1], tmp[4], "-+"))
			my_set_cytoAttrib.add((tmp[0] + " (interacts with) " + tmp[1], tmp[4], "-+"))
		elif sign1[i] == "1" and sign2[i] == "-1":
			print sign1[i], "\t", sign2[i]
			my_set.add((tmp[0], tmp[1], tmp[4], "+-"))
			my_set_cytoAttrib.add((tmp[0] + " (interacts with) " + tmp[1], tmp[4], "+-"))			


out1 = open(out + mod + "_diffNetw_inweb.tsv", "w")
for line in my_set:
	out1.write("\t".join(list(line)))
	out1.write("\n")
out1.close()

# make cytoscape file
out2 = open(out + mod + "_diffNetw_inweb_cytoAttrib.tsv", "w")
for line in my_set_cytoAttrib:
	out2.write("\t".join(list(line)))
	out2.write("\n")
out2.close()
