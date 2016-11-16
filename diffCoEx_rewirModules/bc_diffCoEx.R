#required libraries
library(WGCNA)          ###used for topological overlap calculation and clustering steps
library(RColorBrewer)   ###used to create nicer colour palettes
library(preprocessCore) ###used by the quantile normalization function
library(flashClust)
library(DESeq2)
library(flashClust)


# color
red = rgb(202,49,23, max = 255)   #SLnoa
blue = rgb(0,88,190, max = 255)   #SLmut
yellow="goldenrod2" # GIN NOA
grey="grey67"


# working dir
#setwd("/Users/jessicaxinhu/Documents/noa_bc_pilot/diffCoExpress_analysis/diffCoEx/")
setwd("/Users/jessicaxinhu/Documents/noa_bc_pilot/diffCoExpress_analysis/diffCoEx/with_vest")

# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);
# Allow multi-threading within WGCNA. At present this call is necessary.
# Any error here may be ignored but you may want to update WGCNA if you see one.
# Caution: skip this line if you run RStudio or other third-party R environments.
# See note above.
enableWGCNAThreads()



########################## DATA
rawData<-as.matrix(read.table(file="/Users/jessicaxinhu/Documents/noa_bc_pilot/bc_rnnaseq.tsv",row.names=1,sep="\t",header=T))



########################## PRE-PROCESS

## norm rnaseq data (recommended by WGCNA)
#removing all features that have a count of less than say 10 in more than 90% of the samples
n_samples = 1100
rowsToOmit = c()
rowsToKeep = c()
for (i in 1:dim(rawData)[1]) {
  SamplesLessThanTen = 0
  for (j in 1:dim(rawData)[2]) {
    # select min count:
    if (rawData[i,j] < 10) {    
      SamplesLessThanTen = SamplesLessThanTen + 1
    }
  }
  if (SamplesLessThanTen > n_samples * 0.9) {
    rowsToOmit <- c(rowsToOmit, i)
  } else {
    rowsToKeep <- c(rowsToKeep, i)
  }
}
rawData <- rawData[rowsToKeep,]





### 1. NORMALIZATION method: deseq2 step (global adjustment)
# varianceStabilizingTransformation: observed variability in global properties 
# are due only to technical reasons and are unrelated to the biology of interest
normData <- varianceStabilizingTransformation(rawData)
dimnames(normData)<-dimnames(rawData)
#normData0 <- varianceStabilizingTransformation(rawData)
#save(normData, file="normDat_n1100_brca.RData")
#load("/Users/jessicaxinhu/Documents/noa_bc_pilot/diffCoExpress_analysis/diffCoEx/normDat_n1100_brca.RData")
#write.table(normData0, file="/Users/jessicaxinhu/Documents/noa_bc_pilot/diffCoExpress_analysis/normDat_n1100_brca.tsv", sep="\t", quote=F)


### 2. NORMALIZATION method - remove unwanted technical variation (can use this instead of previous normalization)
# Quantile normalization is a global adjustment method that assumes the statistical distribution of each sample is the same.
#normData<-normalize.quantiles(log2(rawData))
#dimnames(normData)<-dimnames(rawData)



#we remove the probeset of zero variance
RowVar <- function(x) {
  rowSums((x - rowMeans(x))^2)/(dim(x)[2] - 1)
}
any( RowVar(normData) == 0 )


# transpose: pid:rows, genes:cols
normData <- t(normData)
# make pid matchable
rownames(normData) <- gsub("\\.", "-", rownames(normData))



# groups from diffExpr_incl_mutDat.R script:
#targets <- read.table(file="/Users/jessicaxinhu/Documents/noa_bc_pilot/edgeR_DE_analysis/GIN_NOA/targets.tsv", header=T, row.names=1)
load("targets.RData")

highMB1pids <- rownames(targets[targets[,3] == "highMB.1",])    # maf: 251,   vest: 214
lowMB0pids <- rownames(targets[targets[,3] == "lowMB.0",])      # maf: 69,    vest: 117


# subset normdata only for collected pids
collected_pids_chosen <- c(highMB1pids, lowMB0pids)
normData <- normData[substr(rownames(normData), 1, 16) %in% collected_pids_chosen,]



# trait vector
mutstatus = targets[rownames(targets) %in% collected_pids_chosen, 3, drop=FALSE]  # 1 is highMB.1,   0 is lowMB.0
a <- as.data.frame(sapply(mutstatus, function(x) gsub("highMB.1", "1", x)))
a <- as.data.frame(sapply(a, function(x) gsub("lowMB.0", "0", x)))
rownames(a) <- rownames(mutstatus)
colnames(a) <- "status"
mutstatus <- a


# groups of patients with spec trait
datC1 <- normData[substr(rownames(normData), 1, 16) %in% highMB1pids, ] ### these samples correspond to highMB (stress phenotype) group (n=251)
datC2 <- normData[substr(rownames(normData), 1, 16) %in% lowMB0pids,] ###those samples correspond to non-stress phenotype group (n=69)




#####################################################################
######################### Applying DiffCoEx #########################
#####################################################################



############# 1. Network construction (WGCNA)



#########################################################################################################
##### Choose a set of soft-thresholding powers (CHOOSE BETA) --> DO it in R (not Rstudio), or else error
#########################################################################################################
powers = c(c(1:10), seq(from = 12, to=20, by=2))
# Call the network topology analysis function
#sft = pickSoftThreshold(datC1, powerVector = powers, verbose = 5)
bsize = 5000 # only useful for big data
sft = pickSoftThreshold(normData, networkType="signed", corFnc="bicor", verbose=3, powerVector=c(seq(8, 30, by = 2)), blockSize=bsize)   # choose the power at SFT.R.sq = 0.8

# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
################################################################################################
################################################################################################




# parameter for soft thresholding
# NB!!! If settings are user defined parameter for soft thresholding -> need to do statistical significance testing

#softPower=6
softPower=12    #(when using VST normalization on 320 samples)  



# Fast calc if pearson correlation
# Calc gene-gene pair similarities
AdjMatC1<-sign(cor(datC1,method="spearman")) * (cor(datC1,method="spearman"))^2
AdjMatC2<-sign(cor(datC2,method="spearman"))*(cor(datC2,method="spearman"))^2
#save(AdjMatC1, file="AdjMatC1.RData")
#save(AdjMatC2, file="AdjMatC2.RData")
#load("AdjMatC1.RData")
#load("AdjMatC2.RData")
diag(AdjMatC1)<-0
diag(AdjMatC2)<-0
collectGarbage()

AdjDiff = (abs(AdjMatC1-AdjMatC2)/2)^(softPower/2)
AdjDiff_subset <- AdjDiff[substr(rownames(AdjDiff), 1, 16) %in% collected_pids_chosen,]
#save(AdjDiff, file="AdjDiff.RData")
#load("AdjDiff.RData")
#AdjDiff = AdjMatC1-AdjMatC2


# signed co-expr netw -> modules correspond to pos corr genes
# Make differential matrix
#dissTOMC1C2=TOMdist((abs(AdjMatC1-AdjMatC2)/2)^(softPower/2))         # unsigned
dissTOMC1C2=TOMdist((abs(AdjMatC1-AdjMatC2)/2)^(softPower/2), TOMType = "signed")     # Calc TOM dissimilarity network, UNSIGNED is default. Low value -> high similarity -> gene A and B both have sign corr changes with same group of genes
collectGarbage()
#save(dissTOMC1C2, file="dissTOMC1C2.RData")
#load("dissTOMC1C2.RData")







############## 2. Module detection (diffCoEx)

#hierarchical clustering of the genes based on the TOM dissimilarity measure
#Hierarchical clustering is performed using the Topological Overlap of the adjacency difference as input distance matrix
geneTreeC1C2 = flashClust(as.dist(dissTOMC1C2), method = "average");

# Plot the resulting clustering tree (dendrogram)
#png(file="hierarchicalTree.png",height=1000,width=1000)
plot(geneTreeC1C2, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity - BRCA",labels = FALSE, hang = 0.04);
dev.off()


# We now extract modules from the hierarchical tree. This is done using cutreeDynamic. Please refer to WGCNA package documentation for details
# NB!! If settings are user-defined -> need to do statistical significance testing
minModuleSize = 30
dynamicModsHybridC1C2 = cutreeDynamic(dendro = geneTreeC1C2, 
                                      distM = dissTOMC1C2,
                                      method="hybrid", deepSplit = 2,
                                      pamRespectsDendro = FALSE, 
                                      minClusterSize = minModuleSize);



#Every module is assigned a color. Note that GREY is reserved for genes which do not belong to any differentially coexpressed module
dynamicColorsHybridC1C2 = labels2colors(dynamicModsHybridC1C2)

#the next step merges clusters which are close (see WGCNA package documentation)
mergedColorC1C2<-mergeCloseModules(rbind(datC1,datC2),dynamicColorsHybridC1C2,cutHeight=.2)$colors
colorh1C1C2<-mergedColorC1C2
#reassign better colors
colorh1C1C2[which(colorh1C1C2 =="midnightblue")]<-"red"
colorh1C1C2[which(colorh1C1C2 =="lightgreen")]<-"yellow"
colorh1C1C2[which(colorh1C1C2 =="cyan")]<-"orange"
colorh1C1C2[which(colorh1C1C2 =="lightcyan")]<-"green"



# Plot the dendrogram and colors underneath
#png(file="module_assignment.png",width=1000,height=1000)
plotDendroAndColors(geneTreeC1C2, cbind(dynamicColorsHybridC1C2, colorh1C1C2), "Hybrid Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors cells - BRCA")
dev.off()




#We write each module to an individual file containing affymetrix probeset IDs
modulesC1C2Merged<-extractModules(colorh1C1C2,datC1,dir="modules_beta12",file_prefix=paste("Output","Specific_module",sep=''),write=T)
#write.table(colorh1C1C2,file="module_assignment.txt",row.names=F,col.names=F,quote=F)


#We plot to a file the comparative heatmap showing correlation changes in the modules
#The code for the function plotC1C2Heatmap and others can be found below under the Supporting Functions section
plotC1C2Heatmap(colorh1C1C2,AdjMatC1,AdjMatC2, datC1, datC2, brewercol="RdGy")
#png(file="exprChange.png",height=500,width=500)
plotExprChange(datC1,datC2,colorh1C1C2)
#png(file="exprChange_bar_beta12.png",height=500,width=500)
dev.off()








##########################################################################################################################################################################
##########################################################################################################################################################################
################ 3. Gene and module significance (WGCNA)

nGenes = ncol(normData)
nSamples = nrow(normData)
# trait
status <- mutstatus
# merged module colors
moduleColors = colorh1C1C2
mergedColors <- moduleColors
# unmerged module colors
dynamicColors = dynamicColorsHybridC1C2




### recalc MEs with color labels
MEs0 = moduleEigengenes(AdjDiff, moduleColors)$eigengenes  #should use AdjDiff instead of normData since our significant modules are those that are correlated in respect to the CHANGE in coexpr. 
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, as.numeric(status[,1]), use='p')
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)
sizeGrWindow(10,6)
# Will display correlations and their p-values
textMatrix =  paste(signif(moduleTraitCor, 2), "\n(",
                    signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(status),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))



############# MODULE significance ####################################
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(normData, MEs, use = "p"));
#write.table(geneModuleMembership, "MMCorr.txt", sep="\t", quote=F)
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");
#write.table(MMPvalue, "MMPvalue.txt", sep="\t", quote=F)




############### GENE significance ####################################
geneTraitSignificance = as.data.frame(cor(normData, status, use = "p"));
#write.table(geneTraitSignificance, "GSCorr.txt", sep="\t", quote=F)
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
names(geneTraitSignificance) = paste("GS.", names(status), sep="");
names(GSPvalue) = paste("p.GS.", names(status), sep="");
#write.table(GSPvalue, "GSPvalue.txt", sep="\t", quote=F)





# Convert numeric lables into colors. Convert a vector or array of numerical labels in colors.
#In this case we want to convert the geneTraitSignificance (GS) into colors and then add this
#colors to the geneTree.
#GScolors = labels2colors(geneTraitSignificance, zeroIsGrey=FALSE)
#GScolors = numbers2colors(as.matrix(geneTraitSignificance), signed = TRUE,
#                lim = NULL, commonLim = FALSE, naColor = "grey")

#If signed is TRUE, the default setting is to use a palette that starts with green for the most negative values, continues
#with white for values around zero and turns red for positive values. 
signed=TRUE
GScolors= numbers2colors(
  as.matrix(geneTraitSignificance),
  signed,
  centered = signed,
  lim = NULL,
  colors = if (signed) greenWhiteRed(100) else greenWhiteRed(100)[50:100],
  naColor = "grey")
# Plot the dendrogram and colors underneath
pdf(file = "geneTree_moduleSize40_mergedColors_vs_GScolors.pdf", width = 12, height = 6);
sizeGrWindow(12,6)
plotDendroAndColors(geneTreeC1C2, cbind(mergedColors, GScolors),
                    c("Module Colors", "Gene Significance"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()







####Measure of module significance as average gene significance#####
###One can also define a measure of module significance as the average gene 
#significance of all genes in the module. We use the absolute value for 
#defining a correlation based gene significance measure
#GS1=as.numeric(cor(y,datExpr, use="p"))
GeneSignificance=as.matrix(abs(geneTraitSignificance))
colorGS <- as.matrix(cbind(colorh1C1C2,GeneSignificance))
colorGS1 = colorGS
colorGS1 = as.matrix(colorGS[,-2])
# Next module significance is defined as average gene significance.
ModuleSignificance=tapply(GeneSignificance, colorGS1, mean, na.rm=T)
#dim(GeneSignificance)
#dim(colorGS1)

sizeGrWindow(8,7)
par(mfrow = c(1,1))
plotModuleSignificance(GeneSignificance, colorGS1)







#We now create a plot that explains the relationships between modules (heatmap) 
#and the corresponding module eigengene (barplot):
sizeGrWindow(8,7);
#which.module="saddlebrown"
which.module="darkgreen"
MEbarplot = MEs[, paste("ME",which.module, sep="")]
par(mfrow=c(2,1), mar=c(0.3, 5.5, 3, 2))
plotMat(t(scale(normData[,mergedColors==which.module ]) ),
        nrgcols=30,rlabels=F,clabels=F,rcols=which.module,
        main=which.module, cex.main=2)
par(mar=c(5, 4.2, 0, 0.7))
barplot(MEbarplot, col=which.module, main="", cex.main=2,
        ylab="eigengene expression",xlab="array sample")
#Note that the module eigengene takes on low values in arrays where a lot of module genes
#are under-expressed (green color in the heatmap). The ME takes on high values for arrays where a lot of module
#genes are over-expressed (red in the heatmap). ME can be considered the most representative gene expression profile
#of the module.












##########################################################################
## Write network file
##########################################################################
AdjMatC1_sign = sign(AdjMatC1)
AdjMatC2_sign = sign(AdjMatC2)


# Select module
#setwd = "/Users/jessicaxinhu/Documents/noa_bc_pilot/diffCoExpress_analysis/diffCoEx/beta12/diff_modules_beta12_all"
setwd("/Users/jessicaxinhu/Documents/noa_bc_pilot/diffCoExpress_analysis/diffCoEx/with_vest/diffmodules/")


# loop to take module of TOM
for (i in 1:length(unique(colorh1C1C2))) {
  module = unique(colorh1C1C2)[i];
  genes = colnames(normData)
  inModule = (colorh1C1C2==module);
  modGenes = genes[inModule];
  
  #Take module of TOM
  modTOM = dissTOMC1C2[inModule, inModule];
  dimnames(modTOM) = list(modGenes, modGenes)
  # Export the network into an edge list file VisANT can read
  vis1 = exportNetworkToVisANT(modTOM,
                               file = paste(module, "-TOM-beta12", ".txt", sep=""),
                               weighted = TRUE,
                               threshold = 0 )
  
  signMat1 = AdjMatC1_sign[inModule, inModule]
  dimnames(signMat1) = list(modGenes, modGenes)
  signMat2 = AdjMatC2_sign[inModule, inModule]
  dimnames(signMat2) = list(modGenes, modGenes)
  vis2 = exportNetworkToVisANT(signMat1,
                               file = paste(module, "-C1sign-TOM-beta12", ".txt", sep=""),
                               weighted = TRUE,
                               threshold = 0 )
  vis3 = exportNetworkToVisANT(signMat2,
                               file = paste(module, "-C2sign-TOM-beta12", ".txt", sep=""),
                               weighted = TRUE,
                               threshold = 0 )
  
}


# NB!! turquoise module prints vis2 and vis3 in two different lengths. Don't know why.. Have used diff to find the interaction differences in two files. 

module = "turquoise";
genes = colnames(normData)
inModule = (colorh1C1C2==module);
modGenes = genes[inModule];
# Select the corresponding Topological Overlap


#Take module of differential Adj matrix, co-expression value is DIFFERENCE between co-expr values for two conditions
#modTOM = AdjDiff[inModule, inModule];

  






#Take module of TOM
modTOM = dissTOMC1C2[inModule, inModule];
dimnames(modTOM) = list(modGenes, modGenes)
# Export the network into an edge list file VisANT can read
vis1 = exportNetworkToVisANT(modTOM,
                             file = paste(module, "-TOM-beta12", ".txt", sep=""),
                             weighted = TRUE,
                             threshold = 0 )

signMat1 = AdjMatC1_sign[inModule, inModule]
dimnames(signMat1) = list(modGenes, modGenes)
signMat2 = AdjMatC2_sign[inModule, inModule]
dimnames(signMat2) = list(modGenes, modGenes)
vis2 = exportNetworkToVisANT(signMat1,
                             file = paste(module, "-C1sign-TOM-beta12", ".txt", sep=""),
                             weighted = TRUE,
                             threshold = 0 )
vis3 = exportNetworkToVisANT(signMat2,
                             file = paste(module, "-C2sign-TOM-beta12", ".txt", sep=""),
                             weighted = TRUE,
                             threshold = 0 )

##########################################################################
##########################################################################



















#############################################################################
############ EXTRA FUNCTIONS ################################################
#############################################################################

##extractModules: a function which uses the module assignment list as input and writes individual files with the probeset ids for each module
extractModules<-function(colorh1,datExpr,write=F,file_prefix="",dir=NULL)
{
  module<-list()
  if (!is.null(dir))
  {
    dir.create(dir)
    file_prefix=paste(dir,"/",file_prefix,sep="")
  }
  i<-1
  for (c in unique(colorh1))
  {
    module[[i]]<-colnames(datExpr)[which(colorh1==c)]
    if (write) {write.table( module[[i]],file=paste(file_prefix,"_",c,".txt",sep=""),quote=F,row.names=F,col.names=F)}
    i<-i+1
  }
  names(module)<-unique(colorh1)
  module
}


##EigenGenes : this is used by the plotting function to display close together similar modules based on their eigen values
getEigenGeneValues<-function(datRef,colorh1,datAll)
{
  eigenGenesCoef<-list()
  i<-0
  for (c in unique(colorh1))
  {
    i<-i+1
    eigenGenesCoef[[i]]<-prcomp(scale(datRef[,which(colorh1 == c)]))$rotation[,1]
  }
  names(eigenGenesCoef)<-unique(colorh1)
  values<-NULL
  for( c in unique(colorh1))
  {
    v<-rbind(datAll)[,which(colorh1 == c)] %*%  eigenGenesCoef[[c]]
    values<-cbind(values,sign(mean(v))*v)
  }
  colnames(values)<-unique(colorh1)
  values
}



####plotting function for comparative heatmap
plotC1C2Heatmap<-function(colorh1C1C2,AdjMat1C1,AdjMat1C2, datC1, datC2,ordering=NULL,file="DifferentialPlot.png", brewercol)
{
  if (is.null(ordering))
  {
    h<-hclust(as.dist(1-abs(cor(getEigenGeneValues(datC1[,which(colorh1C1C2!="grey")],colorh1C1C2[which(colorh1C1C2!="grey")],rbind(datC1,datC2)[,which(colorh1C1C2!="grey")])))))
    for (c in h$label[h$order])
    {
      ordering<-c(ordering,which(colorh1C1C2 ==c))
    }
  }
  mat_tmp<-(AdjMat1C1[ordering,ordering])
  mat_tmp[which(row(mat_tmp)>col(mat_tmp))]<-(AdjMat1C2[ordering,ordering][which(row(mat_tmp)>col(mat_tmp))])
  diag(mat_tmp)<-0
  mat_tmp<-sign(mat_tmp)*abs(mat_tmp)^(1/2)
  png(file=file,height=1000,width=1000)
  image(mat_tmp,col=rev(brewer.pal(11,brewercol)),axes=F,asp=1,breaks=seq(-1,1,length.out=12))     # CHANGE THIS IF CHANGE R BREWER COLOR!!!!!!!!!
  dev.off()
  unique(colorh1C1C2[ordering])
}

##This function plots side by side the color bar of module assignments, and the change in mean expression of the modules between the two conditions.
plotExprChange<-function(datC1,datC2, colorhC1C2,ordering=NULL)
{
  if (is.null(ordering))
  {
    h<-hclust(as.dist(1-abs(cor(getEigenGeneValues(datC1[,which(colorh1C1C2!="grey")],colorh1C1C2[which(colorh1C1C2!="grey")],rbind(datC1,datC2)[,which(colorh1C1C2!="grey")])))))
    for (c in h$label[h$order])
    {
      ordering<-c(ordering,which(colorh1C1C2 ==c))
    }
  }
  mycolors<-colorh1C1C2[ordering]
  plot(x=0:length(which(mycolors!="grey")),y=rep(1,length(which(mycolors!="grey"))+1),col="white",axes=F,xlab="",ylab="",ylim=c(0,1))
  rr=c(244,239,225,215,209,193,181,166,151,130,110)
  gg=c(228,204,174,160,146,117,94,58,44,45,45)
  bb=c(176,140,109,105,102,91,84,74,70,68,66)
  MyColours<-NULL
  for ( i in 1:11)
  {
    MyColours=c(MyColours,rgb(rr[i],gg[i],bb[i],maxColorValue=255)  )
  }
  exprDiff<-NULL
  l<-0
  for (c in setdiff(unique(mycolors),"grey"))
  {
    meanC1<-mean(t(datC1)[colnames(datC1)[which(colorh1C1C2 == c)],])
    meanC2<-mean(t(datC2)[colnames(datC2)[which(colorh1C1C2 == c)],])
    exprDiff<-rbind(exprDiff,c(meanC1,meanC2))
    r<-l+length(which(mycolors==c))
    rect(l,0.85,r,1,col=c,border=F)
    rect(l,0,r,.4,col=MyColours[floor(meanC2*2)-10],border="white",lwd=2)
    rect(l,0.4,r,.8,col=MyColours[floor(meanC1*2)-10],border="white",lwd=2)
    l<-r
  }
  exprDiff
}

#plotMatrix is a function used to make Additional File 2: Figure S1 plot displaying the
# permutation results.
plotMatrix<-function(mat)
{
  mat[which(row(mat)>col(mat))]<-1001
  image(mat,col=c(gray.colors(4),"white"),breaks=c(0,0.1,50,100,1000,1001),xaxt='n',yaxt='n',xlim=c(-0.2,1.2),ylim=c(-0.2,1.2),bty='n',asp=1)
  text(0:(nrow(mat)-1)/(nrow(mat)-1),1.1,rownames(mat),cex=1,col=rownames(mat))
  text(-0.15,0:(ncol(mat)-1)/(ncol(mat)-1),colnames(mat),cex=1,col=colnames(mat))
  text(apply(matrix(0:(nrow(mat)-1)/(nrow(mat)-1)),1,rep,ncol(mat)),rep(0:(ncol(mat)-1)/(ncol(mat)-1),nrow(mat)),as.numeric(t(mat)),col="white",cex=1.5)
}



