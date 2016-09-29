import sys
import re


# INPUT
inweb = open(sys.argv[1])
cutoff = 0


for line in inweb:
	tmp = line.strip().split("\t")
	m_pattern = "uniprotkb:(.+?)\(gene\ name\)"
	score = float(tmp[14].split("|")[0])
	if score >= cutoff:
		try:	
			gene1 = re.match(m_pattern, tmp[4]).group(1)
			gene2 = re.match(m_pattern, tmp[5]).group(1)
			print gene1 + "\t" + gene2 + "\t" + str(score)
		except:
			pass






