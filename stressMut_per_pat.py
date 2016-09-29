# calculates a table of number of Stress-related SNV and focal CNV per patient


allStress_genes = set( line.strip() for line in open('/Users/jessicaxinhu/Documents/noa_project/reactome/GIN_NOA_analysis/genomic_stress_genes_reactome.tsv')) 
maf = open('/Users/jessicaxinhu/Documents/noa_bc_pilot/edgeR_DE_analysis/GIN_NOA/bc_maf_short.tsv')
commonpat = set( line.strip() for line in open('/Users/jessicaxinhu/Documents/noa_bc_pilot/scripts/common_pat.txt'))
segfile = open('/Users/jessicaxinhu/Documents/noa_project/TCGA_data/segfiles_combined_wGeneid.tsv')


'''
######################## SNV #######################################

stressMut_per_pat = {}

for line in maf:
        tmp = line.strip().split('\t')
        #print tmp[1][0:16]
        if tmp[1][0:16] in commonpat:
                if tmp[0] in allStress_genes: 
                        if tmp[1] in stressMut_per_pat:
                                stressMut_per_pat[tmp[1]] += 1
                        else:
                                stressMut_per_pat[tmp[1]] = 1
                else:
                        if tmp[1] not in stressMut_per_pat:
                                stressMut_per_pat[tmp[1]] = 0

for key in stressMut_per_pat.keys():
        print(key + "\t" + str(stressMut_per_pat[key]))
                
'''

####################### CNV ##########################################

stressCNV_per_pat = {}
segmean_cutoff = -0.2
for line in segfile:
        tmp = line.strip().split('\t')
        gene = tmp[1]
        pid = tmp[0]
        segmean = float(tmp[5])
        if (pid[0:16] in commonpat and segmean <= segmean_cutoff):
                if gene in allStress_genes:
                        if pid in stressCNV_per_pat:
                                stressCNV_per_pat[pid] += 1
                        else:
                                stressCNV_per_pat[pid] = 1
                else:
                        if pid not in stressCNV_per_pat:
                                stressCNV_per_pat[pid] = 0

for key in stressCNV_per_pat.keys():
        print(key + "\t" + str(stressCNV_per_pat[key]))








