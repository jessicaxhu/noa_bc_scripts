import sys


maf = open(sys.argv[1])


maf.readline()
for line in maf:
	tmp = line.strip().split("\t")
	print "%s\t%s%s\t%s\t%s\t%s\t%s\t%s" % (tmp[0], "chr", tmp[4], tmp[5], tmp[7], tmp[11], tmp[12], tmp[16])
	



