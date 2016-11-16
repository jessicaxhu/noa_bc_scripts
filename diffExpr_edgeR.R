
#######################################
#### edgeR differential analysis ######
#######################################

library(edgeR)
library(RColorBrewer)
library(devtools)



# color
red = rgb(202,49,23, max = 255)   #SLnoa
blue = rgb(0,88,190, max = 255)   #SLmut
yellow="goldenrod2" # GIN NOA

grey="grey67"

####### LOAD DATA ############################################################

# all data
rawdata <- read.delim("~/Documents/noa_project/GSE62944/GSM1536837_06_01_15_TCGA_24.tumor_Rsubread_FeatureCounts.txt", check.names=FALSE, stringsAsFactors=FALSE)


# taking BRCA patients
brca_patlist = as.character(read.table("~/Documents/noa_bc_pilot/edgeR_DE_analysis/GIN_NOA/patid_BRCA_rnaseqV2.tsv")[,1])
bc <- rawdata[, colnames(rawdata) %in% brca_patlist]
rownames(bc) <- rawdata[,1]
#bc <- cbind(rawdata[,1], bc)
#colnames(bc)[1] <- "Genes"




#######################################################
############# FUNCTIONS ###############################
#######################################################

short_id <- function (x) {
	substr(x, 1, 16)
}




##################################################
####### chromosome MB data #######################
##################################################



# all CHROMOSOME data 
chromo_mb = read.table(file="~/Documents/noa_project/muttables/cnv_mb_uniq.tsv", sep="\t", stringsAsFactors=FALSE)

# only bc
mb_brca <- chromo_mb[ substr(as.character(chromo_mb[,1]), 1, 16) %in% substr(brca_patlist, 1, 16) , ]


#lower = quantile(mb_brca[,2])[2]
#upper = quantile(mb_brca[,2])[4]
lower = 5
upper = 20
plot(density(mb_brca[,2]), main="BRCA CNV mutation burden", xlab = "% Nucleotides with CNV loss or gain", cex.lab=1.5, cex.axis=1.5)
abline(v=upper, col = "grey", lwd=3, lty=3)
abline(v=lower, col = "grey", lwd=3, lty=3)
#median(mb_brca[,2])

# count number of TRUE when %-wise nucleotides -0.2 > segmean > 0.2
summary(mb_brca[,2] >= upper)
summary(mb_brca[,2] < lower)



#################################################
####### MAF mutation data #######################
#################################################

# all mutation types (including silent)
#maf <- read.table(file="/Users/jessicaxinhu/Documents/noa_bc_pilot/edgeR_DE_analysis/GIN_NOA/bc_maf_short.tsv", header=F)
vest <- read.table(file="/Users/jessicaxinhu/Documents/noa_bc_pilot/vest3/output/jessica_hu_20161109_072209/Variant.Result.tsv", skip=11, header=T, sep="\t")





# COMMON PATIENT LIST between rnaseqV2 and cnv (16 char TCGA id)
common_pat <- intersect(substr(colnames(bc), 1, 16), substr(as.character(mb_brca[,1]), 1, 16))
#common_pat <- intersect(common_pat, substr(as.character(maf[,2]), 1, 16))      # 960 pat left
common_pat <- intersect(common_pat, substr(as.character(vest[,8]), 1, 16))      # 960 pat left



# take only patients that INTERSECT 
# RNAseq
bc = bc[, substr(colnames(bc), 1, 16) %in% common_pat]
# chromo
mb_brca <- mb_brca[ substr(mb_brca[,1], 1, 16) %in% common_pat, ]
# mut
#maf <- maf[ substr(maf[,2], 1, 16) %in% common_pat, ]
vest <- vest[ substr(vest[,8], 1, 16) %in% common_pat, ]




################################################################
### make mutation stress vectors for stress phenotype
################################################################


dnaRepairStress_genes <- read.table(file="/Users/jessicaxinhu/Documents/noa_project/reactome/GIN_NOA_analysis/dnaRepairStress_genes_reactome.tsv", header=F)
mitoticStress_genes <- read.table(file="/Users/jessicaxinhu/Documents/noa_project/reactome/GIN_NOA_analysis/mitoticStress_genes_reactome.tsv", header=F)
allStress_genes <- read.table(file="/Users/jessicaxinhu/Documents/noa_project/reactome/GIN_NOA_analysis/genomic_stress_genes_reactome.tsv", header=F)



# maf
dnaRepairStress_patients = c()
mitoticStress_patients = c()
allStress_patients = c()
nStressMut_per_pat = c()
for (i in 1:dim(maf)[1]) {
  if ( maf[i,1] %in% allStress_genes[,1] ) {
    allStress_patients <- c(allStress_patients, as.character(maf[i,2]))
  }
  if ( maf[i,1] %in% dnaRepairStress_genes[,1] ) {
    dnaRepairStress_patients <- c(dnaRepairStress_patients, as.character(maf[i,2]))
  }
  if ( maf[i,1] %in% mitoticStress_genes[,1] ) {
    mitoticStress_patients <- c(mitoticStress_patients, as.character(maf[i,2]))
  }
}
length(unique(allStress_patients))



# vest
dnaRepairStress_patients = c()
mitoticStress_patients = c()
allStress_patients = c()
nStressMut_per_pat = c()
for (i in 1:dim(vest)[1]) {
  if (!is.na(vest[i,12])) {
    if ((vest[i,9] %in% allStress_genes[,1]) && (vest[i,12] <= 0.05)) {
      allStress_patients <- c(allStress_patients, as.character(vest[i,8]))
      }
    if ((vest[i,9] %in% dnaRepairStress_genes[,1]) && (vest[i,12] <= 0.05)) {
      dnaRepairStress_patients <- c(dnaRepairStress_patients, as.character(vest[i,8]))
      }
    if ((vest[i,9] %in% mitoticStress_genes[,1]) && (vest[i,12] <= 0.05)) {
      mitoticStress_patients <- c(mitoticStress_patients, as.character(vest[i,8]))
    }
  }
}
length(unique(allStress_patients))



allStress_patients <- sort(unique(allStress_patients))  # maf: n_mut = 764, n_noMut = 196, vest: n_mut=587, n_nomut=960-587=373
dnaRepairStress_patients <- sort(unique(dnaRepairStress_patients))
mitoticStress_patients <- sort(unique(mitoticStress_patients))
#length(intersect(dnaRepairStress_patients, mitoticStress_patients))




allStress_patients_df = data.frame("pid" = character(1), "mut_stat" = numeric(1), stringsAsFactors=FALSE)
for (i in 1:length(common_pat)) {
  if (common_pat[i] %in% substr(allStress_patients, 1, 16)) {
    allStress_patients_df <- rbind(allStress_patients_df, c(common_pat[i], 1))
  }
  else {
    allStress_patients_df <- rbind(allStress_patients_df, c(common_pat[i], 0))
  }
}
# remove first row
allStress_patients_df <- allStress_patients_df[-1, ]


# SORT to alphabetically to fit patid from DGE object
allStress_patients_df <- allStress_patients_df[order(allStress_patients_df[,1]),] 		# vector with 1 if mut in at least 1 stress pw



################################################################
### make CNV stress vectors for stress phenotype
################################################################


cnv_mb = vector()
for(i in 1:length(mb_brca[,2])) {
	if (mb_brca[i,2] >= 20) {
		cnv_mb  <- c(cnv_mb, "highMB")			# 1: high MB
	} else if( mb_brca[i,2] < 5) {
		cnv_mb <-c(cnv_mb, "lowMB")			# 0: low MB
	} else {
		cnv_mb <- c(cnv_mb, "mediumMB")			# 2: medium MB
	}
}
mb_class <- cbind(mb_brca, cnv_mb)
colnames(mb_class) <- c("patient_id", "cnv_mb", "cnv_groups")
#write.table(mb_class, file="~/Documents/noa_bc_pilot/cnv_mb_groups.tsv", quote=F, sep="\t", row.names=F)



# SORT to alphabetically to fit patid from DGE object
mb_class <- mb_class[order(mb_class[,1]),] 		






###############################################################################
############ archetype group 3 analysis
###############################################################################
arch3 <- scan(file="/Users/jessicaxinhu/Documents/noa_bc_pilot/archetype_analysis/arch3_pidlist.txt", what="character")
arch_lumB <- scan(file="/Users/jessicaxinhu/Documents/noa_bc_pilot/archetype_analysis/arch_lumB_pidlist.txt", what="character")
arch_basal <- scan(file="/Users/jessicaxinhu/Documents/noa_bc_pilot/archetype_analysis/arch_basal_pidlist.txt", what="character")

arch3_cnvValues <- mb_class[substr(mb_class[,1],1,15) %in% arch3,2]
arch_lumB_cnvValues <- mb_class[substr(mb_class[,1],1,15) %in% arch_lumB,2]
arch_basal_cnvValues <- mb_class[substr(mb_class[,1],1,15) %in% arch_basal,2]
highMB_cnvValues <- mb_class[mb_class[,3] == "highMB",2]


blue = rgb(56, 86, 165, maxColorValue = 265)

plot(density(mb_class[,2]), col=grey, ylim=c(0,0.1), lwd = 3, main="% CNV mutation burden distributions of archetypes")
lines(density(highMB_cnvValues), col="palegreen3", lwd = 3)
lines(density(arch_lumB_cnvValues), col="darkorange", lwd = 3)
lines(density(arch_basal_cnvValues), col=blue, lwd = 3)
lines(density(arch3_cnvValues), col="red2", lwd = 3)

legend("topleft", c("All BC", "HighMB group", "Archetype LumB", "Archetype Basal", "Archetype 3"), lty=c(1,1,1,1,1,1),lwd=c(3,3,3,3,3,3),col=c(grey,"palegreen3","darkorange",blue,"red2"), bty="n")
###############################################################################
###############################################################################







###############################################################################
# DE analysis begins
###############################################################################

y <- DGEList(counts=bc, genes=rownames(bc))		# DGE object
colnames(bc) <- substr(colnames(bc), 1, 16)

# Filter out lowly expressed genes with low counts
keep <- rowSums(cpm(y)>1) >= 2
y <- y[keep, , keep.lib.sizes=FALSE]


y$samples$lib.size <- colSums(y$counts)		# recompute lib sizes
#rownames(y$counts) <- y$genes$genes


# TMM normalization
y <- calcNormFactors(y)		# norm differences between lib
head(y$samples)




### RAN TO HERE:

#################################################################################
### SIMPLE APPROACH (with one factor: CNV)
#################################################################################


#et <- exactTest(y, pair=c("A","B"))
#topTags(et)



#################################################################################
### ADVANCED APPROACH (with two factors: SNV and CNV)
#################################################################################

#setwd("/Users/jessicaxinhu/Documents/noa_bc_pilot/edgeR_DE_analysis/GIN_NOA")
#setwd("/Users/jessicaxinhu/Documents/noa_bc_pilot/edgeR_DE_analysis/GIN_NOA/with_vest")

# design matrix
targets <- as.data.frame(cbind(as.character(mb_class[,3]), allStress_patients_df[,2]))
rownames(targets) <- allStress_patients_df[,1]
colnames(targets) <- c("cnv", "stress_pw_mut")
Group <- factor(paste(targets$cnv,targets$stress_pw_mut,sep="."))
targets <- cbind(targets,Group=Group)
#write.table(targets, file="targets.tsv", quote=F, sep="\t")
#save(targets, file="targets.RData")
design <- model.matrix(~0+Group)
colnames(design) <- levels(Group)
# estimate dispersion
y <- estimateDisp(y, design)		# dispersion




fit <- glmFit(y, design)



# Make contrasts
my.contrasts <- makeContrasts(
  highMB.1vs0 = highMB.1 - highMB.0,
  highMB.1vslowMB.0 = highMB.1 - lowMB.0, # all high versus low
  mediumMB.1vslowMB.0 = mediumMB.1 - lowMB.0,
  highMB.0vslowMB.1 = highMB.0 - lowMB.1,
  
  lowMB.1vs0 = lowMB.1 - lowMB.0,  # pw specific NOAs
  highMB.1vslowMB.1 = highMB.1 - lowMB.1, # GIN NOA
  levels = design
)

lrt <- glmLRT(fit, contrast=my.contrasts[,"highMB.1vslowMB.0"])
lrt <- glmLRT(fit, contrast=my.contrasts[,"highMB.1vs0"])
lrt <- glmLRT(fit, contrast=my.contrasts[,"highMB.0vslowMB.1"])
lrt <- glmLRT(fit, contrast=my.contrasts[,"mediumMB.1vslowMB.0"])

lrt <- glmLRT(fit, contrast=my.contrasts[,"lowMB.1vs0"])
lrt <- glmLRT(fit, contrast=my.contrasts[,"highMB.1vslowMB.1"])


out <- topTags(lrt, n=Inf, adjust.method="BH")
#keep <- out$table$FDR <= 0.1 & abs(out$table$logFC) >= 1
keep <- out$table$FDR <= 0.001 & out$table$logFC >= 1          # THRESHOLD values 
DE_table <- as.data.frame(out[keep,])
DE_table_out <- DE_table[order(DE_table[,2], decreasing=T),]
dim(DE_table_out)
write.table(DE_table_out, file="highMB.1vslowMB.0.tsv", quote=F, sep="\t", row.names=F)



###############
# Get DE genes in same pathway as stress pw genes

DEup_stressPw_genes <- DE_table_out[DE_table_out[,1] %in% allStress_genes[,1], ]
head(DEup_stressPw_genes)
dim(DEup_stressPw_genes)
write.table(DEup_stressPw_genes[,1], file="highMB.1vslowMB.0_inStressPw.txt", quote=F, row.names=F, col.names=F, sep="\t")

###############################################

# PLOTS

plot(DEall_highMB_o$logFC)



# plot correlation between Stress pw mutations and global CNV burden
# SNV
stressMut_per_pat = read.table(file="/Users/jessicaxinhu/Documents/noa_bc_pilot/scripts/stressMut_per_pat_out.tsv", header=F)
stressMut_per_pat <- stressMut_per_pat[order(stressMut_per_pat[,1]),]
df <- cbind(stressMut_per_pat[,2], mb_class[,2])
rownames(df) <- stressMut_per_pat[,1]
colnames(df) <- c("mut", "cnv")
cutoff_nMut = 25
df1 <- df[df[,1] <= cutoff_nMut ,]
plot(df1, pch=16, xlab="n mut in stress pathways", ylab = "% cnv burden")
abline(h=5, col="blue")
abline(h=20, col="blue")
abline(v=1, col="red")
blue = rgb(0,88,190, max = 255)
boxplot(cnv ~ mut, df, pch=16, xlab="sqrt( n snv in stress pathways )", ylab = "sqrt(  % cnv burden )", main="Correlation between stress pathway SNV and CNV burden", border=blue)




########### PCA ################################################# 

de_genes <- DE_table_out[,1]
df <- bc[rownames(bc) %in% de_genes,]
df <- t(df)
df <- log(df+(df==0))
plot(density(df))
pc <- prcomp(df, center = TRUE, scale = TRUE)
#plot(df_pca, type = "l")  



group_mb = c()
for (i in 1:dim(targets)[1]) {
  if (targets[i,3] == "highMB.1") {
    group_mb <- c(group_mb, "firebrick2")
  }
  else if (targets[i,3] == "lowMB.0") {
    group_mb <- c(group_mb, "royalblue1")
  }
  else {
    group_mb <- c(group_mb, "grey")
  }
}


clin <- read.table("/Users/jessicaxinhu/Documents/noa_bc_pilot/tcga_data/BRCA.clin.merged.picked.txt", sep="\t", header=T, row.names=1)
stage <- t(clin[rownames(clin) == "pathologic_stage",])
rownames(stage) <- gsub("\\.", "-", toupper(rownames(stage)))
stage <- stage[rownames(stage) %in% substr(rownames(df),1,12),]
stage <- stage[order(names(stage)),]
group_stage = c()
for (i in 1:length(stage)) {
  if (is.na(stage[i])) {
    group_stage <- c(group_stage, "grey")
  }
  else if ( (stage[i] == "stage i") || (stage[i] == "stage iia") || (stage[i] == "stage iib") ) {
    group_stage <- c(group_stage, "firebrick2")
  }
  else if ( (stage[i] == "stage iiia") || (stage[i] == "stage iiib") || (stage[i] == "stage iiic") ) {
    group_stage <- c(group_stage, "royalblue1")
  }
  else {
    group_stage <- c(group_stage, "grey")
  }
}

par(mfrow=c(1,2))
# pca with groups as mutburden 
plot(pc$x[, 1], pc$x[, 2], col = group_mb, main = "BC - DE genes", xlab = "PC1", ylab = "PC2", pch=16)
legend(25,-25, legend=c("MB", "non-MB"), col=c("red","blue"), pch=16, bty='n')
# pca with groups as tumor stage
plot(pc$x[, 1], pc$x[, 2], col = group_stage, main = "BC - DE genes", xlab = "PC1", ylab = "PC2", pch=16)
legend(25,-25, legend=c("high stage", "low stage"), col=c("red","blue"), pch=16, bty='n')
