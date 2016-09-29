# packages
import sys
from subprocess import call
import glob



#########################################################################################################################
### LOAD files
#########################################################################################################################


in_folder = sys.argv[1]



####################################
### 2. MUT files
f1 = in_folder + '/vest.tsv'
vest = open(f1)			# VEST3 file


pval_cutoff = 0.05			# vest3 pval cutoff



#####################################
### 3. CNV files

f2 = in_folder + '/gistic.tsv'
gistic = open(f2)


gistic_tr = f2 + "_tr"
call(["awk","-f", "/Users/jessicaxinhu/Documents/my_programs/transpose.awk", f2], stdout = open(gistic_tr,'w'))
gistic_genelist = set()
# load gistic file into genelist set() structure
for line in open(gistic_tr):
	tmp = line.strip().split("\t")
	gist_genes = tmp[4:]
	try:
		qval = float(tmp[1])
		if qval < 0.05:
			gistic_genelist.update(gist_genes)
	except:
		pass
call(["rm", gistic_tr])



f3 = in_folder + '/segfile.tsv'
segfile = open(f3)			# use segfile with gene ids
segmean_cutoff = -0.2			# segmean cutoff





#####################################
### 4. expr data
f4 = in_folder + '/rsem.tsv'
rsem = open(f4)				# logfc matrix --> if upregu --> save 0/1

logfc_cutoff = 0.2





######################################
### 5. modules data
modules = {}
f_modules = open("/Users/jessicaxinhu/Documents/noa_bc_pilot/mutEx_analysis/modules/clusterONE_modules.txt")		# clusterONE






#########################################################################################################################
### ClusterONE modules
#########################################################################################################################
pat_gene_info = {}



modules = {}
module_genes = set()
module_index = 0
for line in f_modules:
	module_index += 1 
	module = line.strip().split("\t")
	for elem in module:
		module_genes.add(elem)
	modules[module_index] = module

#########################################################################################################################
### if mutation in vest --> save 1
#########################################################################################################################

i = 0
mut_info  = {}
mut_patients = set()
for line in vest:
        i += 1
        if i > 12:
                tmp = line.strip().split("\t")
                if len(tmp[11]) > 1:
                        pval = float(tmp[11])
		pid = tmp[7][0:15]
		mut_patients.add(pid)
                gene = tmp[1]
                if pval <= pval_cutoff:
                        mut_info[(pid, gene)] = 1
#########################################################################################################################
### gistic genelist --> look up in segfile if patient has del for any of genes in gistic genelist --> save 1
#########################################################################################################################

cnv_info = {}
cnv_patients = set()
for line in segfile:
	tmp = line.strip().split("\t")
	pid = tmp[0][0:15]
	cnv_patients.add(pid)
	gene = tmp[1]
	segmean = float(tmp[5])

	if segmean <= segmean_cutoff and gene in gistic_genelist:
		cnv_info[(pid, gene)] = 1


#########################################################################################################################
### RSEM matrix --> look up if patient has logfc > 0.2 --> save 1
#########################################################################################################################

noa_info = {}

pid_header = rsem.readline().strip().split("\t")
pid_header.pop(0)
rsem_patients = set()

for line in rsem:
	for i in range(0, len(pid_header)):
		tmp = line.strip().split("\t")
		gene = tmp.pop(0)
		pid = pid_header[i][0:15]
		rsem_patients.add(pid)

		if tmp[i] > logfc_cutoff:
			noa_info[(pid, gene)] = 1






#########################################################################################################################
#########################################################################################################################
#########################################################################################################################


### common pids

common_pids = set()
tmp_intercept = mut_patients.intersection(cnv_patients)
commmon_pids = tmp_intercept.intersection(rsem_patients)
print common_pids



mod_stat = {}

# find common patients in 3 data structures

for pid in common_pids:
	for mod_index in modules.keys():
		mod_stat[mod_index] = [0,0,0,0]		# 0:tNOA_tMUT, 1:tNOA_fMUT, 2:fNOA_tMUT, 3:fNOA_fMUT
		for gene in modules[mod_index]:
			if (pid, gene) in noa_info:
				if (pid, gene) in mut_info or (pid, gene) in cnv_info:
					mod_stat[mod_index][0] += 1
				else:
					mod_stat[mod_index][1] += 1
			elif (pid, gene) in mut_info or (pid, gene) in cnv_info:
				mod_stat[mod_index][2] += 1
			else:
				mod_stat[mod_index][3] += 1



