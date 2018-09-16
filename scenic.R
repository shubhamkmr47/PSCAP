library(scran)        #bioconductor
library(scater)       #bioconductor
library(Biobase)      #bioconductor
library(EBSeq)        #bioconductor
library(scDD)         #bioconductor
library(SummarizedExperiment) # bioconductor
library(ggplot2)      #cran
library(RColorBrewer) #cran
library(SingleCellExperiment)

load("./data/data.RData")

cell_type <- as.vector(cells$Characteristics.cell.type.)
cell_type[cell_type == "alpha cell" | cell_type == "beta cell" |
            cell_type == "gamma cell" | cell_type == "delta cell" |
            cell_type == "epsilon cell"] <- "endocrine cells"
cell_type[cell_type == "acinar cell" |
            cell_type == "ductal cell"] <- "exocrine cells"

cells <- cbind(cells, cell_type)
cells <- cells[order(cells$Characteristics.disease.), ]
cells <- cells[order(cells$Characteristics.cell.type.), ]

# convert data to numeric matrix
ids <- as.vector(cells$Source.Name)
data.count <- data.count[, ids]
data.rpkm <- data.rpkm[, ids]

rm(ids, cell_type)
rownames(cells) <- as.vector(cells$Source.Name)

row <- rownames(data.count)
data.count <- apply(data.count, 2, as.numeric)
rownames(data.count) <- row

row <- rownames(data.rpkm)
data.rpkm <- apply(data.rpkm, 2, as.numeric)
data.rpkm <- log(data.rpkm + 1, 2)
rownames(data.rpkm) <- row
rm(row)

gene.var <- apply(data.rpkm, 1, function(x) IQR(x[x > 0]))
gene.var[which(is.na(gene.var))] <- 0
hist(gene.var[gene.var > 0], main = "gene IQR freq", xlab = "range", breaks = 50)

data.genes <- cbind(data.genes, iqr_rpkm = as.vector(gene.var))

data.rpkm <- data.rpkm[-c(26179:26271), ]

# select genes to do PCA on rpkm data
index <- which(gene.var >= 7)
index <- index [!index %in% c(26179:26271)]
topvar <- data.rpkm[index, ]
dim(topvar)

topgene <- data.genes[which(data.genes$G %in% rownames(topvar)), ]
rpkm.pca <- prcomp(topvar, center = TRUE, scale. = TRUE)
summary(rpkm.pca)$importance[, 1:5]

#Plot the first PCs
df <- data.frame(PC1 = rpkm.pca$rotation[,1], PC2 = rpkm.pca$rotation[,2], cells)

rm(gene.var, rpkm.pca, topvar, index, topgene)
gc()  # clean up workspace

rm(data.rpkm)

name <- as.vector(data.genes$name)
name[c(26179:26270)] <- as.vector(data.genes$id[c(26179:26270)])
data.genes$name <- name
rm(name)

all.equal(as.vector(data.genes$G), rownames(data.count))

dim(data.count)

data.count <- data.count[-c(26271), ]
data.genes <- data.genes[-c(26271), ]

gene.var <- apply(log(data.count + 1, 2), 1, function(x) iqr(x[x > 0]))

# remove duplicates of count data
data.genes <- cbind(data.genes, iqr_count = gene.var)
data.genes <- data.genes[order(data.genes$iqr_count, decreasing = TRUE),]
data.genes <- data.genes[!duplicated(data.genes$name), ]
data.genes <- data.genes[order(as.vector(data.genes$G)), ]

data.count <- data.count[as.vector(data.genes$G), ]
data.count <- data.count[order(rownames(data.count)), ]

all.equal(as.vector(data.genes$G), rownames(data.count))

dim(data.count)
rownames(data.count) <- as.vector(data.genes$name)

# construct a SCESet that also contains the gene counts, ercc counts, and metadata
eset <- newSCESet(countData = data.count, phenoData = AnnotatedDataFrame(cells))
isSpike(eset, "ERCC") <- grepl("^ERCC-", rownames(eset))

which(isSpike(eset))

rm(gene.var, data.count); gc()

eset <- calculateQCMetrics(eset,
        feature_controls = list(ERCC = which(isSpike(eset))))

which(isSpike(eset))

# QC to compare the level of dropout in endogeneous genes to ERCC spike ins in raw data
plotQC(eset, type = "exprs-freq-vs-mean")

# first, filter out genes that are almost always zero 
keep <- rowSums(counts(eset) > 0) >= round(dim(cells)[1]/4)
sum(keep)
eset <- eset[keep,] 

# filter out genes
keep <- rowMeans(counts(eset)) >= 5
sum(keep)
eset <- eset[keep,]

rm(keep)

# normalize counts for library size using the pool & deconvolve method of Lun et al. (2016)
eset <- computeSumFactors(eset, sizes = seq(10, 200, 10))
summary(sizeFactors(eset))

# use the size factors calculated above to normalize the counts
eset <- normalizeSCE(eset)

plot(sizeFactors(eset), colSums(counts(eset))/1e6, log="xy",
     ylab="Library Size (Total Counts in Millions)", xlab="Pooled Size Factor Estimate",
     main="Normalization factor versus library size")

detectionRate <- apply(counts(eset), 2, function(x) sum(x > 0) / length(x))
hist(detectionRate)

plot(detectionRate, df$PC1, pch=20, col="grey", ylab="PC 1")

rm(detectionRate, df); gc()  #clean up workspace

# save.image("./data/normalize.RData")

### scenic

library(SCENIC)
library(GENIE3)
library(AUCell)
library(RcisTarget)

library(mixtools)
library(doMC)
library(doRNG)
library(NMF)
library(Rtsne)
library(data.table)
library(reshape2)

library(scater)

load("./data/normalize.RData")

cells <- cells[cells$Characteristics.cell.type. != "unclassified cell" & 
                 cells$Characteristics.cell.type. != "unclassified endocrine cell", ]
id <- cells$Source.Name

cellLabels <- cells
exprMatrix <- exprs(eset)[, id]
geneNames <- rownames(exprMatrix)

rm(id, cells, eset)

dim(exprMatrix)
eset <- new("ExpressionSet", exprs=exprMatrix,
            phenoData=new("AnnotatedDataFrame",
                          data=data.frame(cellLabels[colnames(exprMatrix),, drop=FALSE])))

unique(cellLabels$Characteristics.cell.type.)

rm(geneNames, exprMatrix, cellLabels)

options(stringsAsFactors=FALSE)

# Infer potential transcription factor targets based on the expression data
exprMat <- exprs(eset)
dim(exprMat)

cellInfo <- pData(eset)[colnames(exprMat), "Characteristics.cell.type.", drop=F]
paste(unique(cellInfo$Characteristics.cell.type.), collapse = "', '")
cells <- c('acinar cell', 'alpha cell', 'beta cell', 'co-expression cell', 
           'delta cell', 'ductal cell', 'endothelial cell', 'epsilon cell', 
           'gamma cell', 'mast cell', 'MHC class II cell', 'PSC cell')
cellsCol <- c("forestgreen", "red3", "darkorange", "magenta4",
              "hotpink", "chartreuse", "skyblue", "darkblue", 
              "brown", "coral", "darkorchid", "darkcyan")

# Color to assign to the variables
colVars <- list(Characteristics.cell.type.=setNames(cellsCol, cells))

cellInfo <- pData(eset)[colnames(exprMat), "Characteristics.individual.", drop=F]
cells <- unique(data.info$Characteristics.individual.)

cellsCol <- c("forestgreen", "red3", "darkorange", "magenta4",
              "hotpink", "chartreuse", "skyblue", "darkblue",
              "brown", "coral")
colVars <- list(Characteristics.cell.type.=setNames(cellsCol, cells))
plot.new(); legend(0,1, fill=colVars$Characteristics.cell.type.,
                   legend=names(colVars$Characteristics.cell.type.))

# http://scenic.aertslab.org/downloads/databases/RcisTarget.hg19.motifDatabases.20k_0.1.1.tar.gz
# http://scenic.aertslab.org/downloads/databases/RcisTarget.mm9.motifDatabases.20k_0.1.1.tar.gz
# install.packages('./data/RcisTarget.mm9.motifDatabases.20k_0.1.1.tar.gz', repos = NULL)
# install.packages('./data/RcisTarget.hg19.motifDatabases.20k_0.1.1.tar.gz', repos = NULL)

org <- "hg19"

if(org=="hg19")
{
  library(RcisTarget.hg19.motifDatabases.20k)
  
  # Get genes in databases:
  data(hg19_500bpUpstream_motifRanking) # or 10kbp, they should have the same genes
  genesInDatabase <- hg19_500bpUpstream_motifRanking@rankings$rn
  rm(hg19_500bpUpstream_motifRanking)
  
  # Get TFS in databases:
  data(hg19_direct_motifAnnotation)
  allTFs <- hg19_direct_motifAnnotation$allTFs
  rm(hg19_direct_motifAnnotation)
}

nCellsPerGene <- apply(exprMat, 1, function(x) sum(x>0))
nCountsPerGene <- apply(exprMat, 1, sum)

summary(nCellsPerGene)
summary(nCountsPerGene)

max(exprMat)

sum(exprMat>0) / sum(exprMat==0)

minReads <- 3*.01*ncol(exprMat)
genesLeft_minReads <- names(nCountsPerGene)[which(nCountsPerGene > minReads)]
length(genesLeft_minReads)

minSamples <- ncol(exprMat)*.01
nCellsPerGene2 <- nCellsPerGene[genesLeft_minReads]
genesLeft_minCells <- names(nCellsPerGene2)[which(nCellsPerGene2 > minSamples)]
length(genesLeft_minCells)

# The genes available on RcisTarget databases will be used
genesLeft_minCells_inDatabases <- genesLeft_minCells[which(genesLeft_minCells %in% genesInDatabase)]
length(genesLeft_minCells_inDatabases)

exprMatrix_filtered <- exprMat[genesLeft_minCells_inDatabases, ]

rm(nCellsPerGene, nCountsPerGene, nCellsPerGene2)
rm(minReads, minSamples, genesInDatabase, genesLeft_minCells, 
   genesLeft_minCells_inDatabases, genesLeft_minReads)

# Potential regulators: List of transcription factors
inputTFs <- allTFs[allTFs%in% rownames(exprMatrix_filtered)]

c(allTFs=length(allTFs), inputTFs=length(inputTFs))

hist(exprMatrix_filtered[exprMatrix_filtered != 0])
dim(exprMatrix_filtered)

# Run GENIE3
set.seed(123)

weightMatrix <- GENIE3(exprMatrix_filtered, regulators=inputTFs, nCores=3, verbose = TRUE)
dim(weightMatrix)

corrMat <- cor(t(exprMatrix_filtered), method="spearman")

# Create co-expression modules
linkList <- getLinkList(weightMatrix, threshold=0.001) # (slighly faster)
linkList <- getLinkList(weightMatrix)
colnames(linkList) <- c("TF", "Target", "weight")
linkList <- linkList[order(linkList[,"weight"], decreasing=TRUE),]

dim(linkList)
head(linkList)
length(unique(linkList$TF))
length(unique(linkList$Target))

# Creating TF modules (potential TF-targets)
quantile(linkList$weight, probs=c(0.5, 0.75, 0.8, 0.90))

plot(linkList$weight[1:1000000], type="l", ylim=c(0, max(linkList$weight)), main="Weight of the links",
     ylab="Weight", xlab="Links sorted decreasingly")
abline(h=0.001, col="blue") # Threshold

sum(linkList$weight>0.001)/nrow(linkList)

# Keep only the links over 0.001:
linkList_001 <- linkList[which(linkList[,"weight"]>0.001),]
nrow(linkList_001) 

tfModules <- list()

linkList_001$TF <- as.character(linkList_001$TF)
linkList_001$Target <- as.character(linkList_001$Target)

# Create TF-modules:
# Weight > 0.001 
tfModules[["w001"]] <- split(linkList_001$Target, factor(linkList_001$TF))

# 2: Weight > 0.005
llminW <- linkList_001[which(linkList_001[,"weight"]>0.005),]
tfModules[["w005"]] <- split(llminW$Target, factor(llminW$TF))

# 3: Top 50 targets for each TF
tfModules[["top50"]] <- lapply(tfModules[["w001"]], function(x) x[1:(min(length(x), 50))])

# 4-6: Top regulators per target 
linkList_001_byTarget <- split(linkList_001, factor(linkList_001$Target))

nTopTfs <- c(5, 10, 50)
nTopTfs <- setNames(nTopTfs, paste("top", nTopTfs, "perTarget", sep=""))

topTFsperTarget <- lapply(linkList_001_byTarget, function(llt) {
  nTFs <- nTopTfs[which(nTopTfs <= nrow(llt))]
  melt(lapply(nTFs, function(x) llt[1:x,"TF"]))
})
topTFsperTarget <- topTFsperTarget[which(!sapply(sapply(topTFsperTarget, nrow), is.null))]
topTFsperTarget.asDf <-  data.frame(rbindlist(topTFsperTarget, idcol=TRUE))
head(topTFsperTarget.asDf)
colnames(topTFsperTarget.asDf) <- c("Target", "TF", "method")

# Merge the all the gene-sets:
tfModules.melted <- melt(tfModules)
colnames(tfModules.melted) <- c("Target", "TF", "method")
tfModules <- rbind(tfModules.melted, topTFsperTarget.asDf)

table(as.vector(tfModules$method))

# corr = 1 if weight >= 0.03
# corr = -1 if weight <= -0.03
# otherwise weigt = 0
tfs <- unique(tfModules$TF)
corrMat <- corrMat[tfs,]

tfModules_byTF <- split(tfModules, factor(tfModules$TF))
tfModules_withCorr_byTF <- lapply(tfModules_byTF, function(tfGeneSets)
{
  tf <- unique(tfGeneSets$TF)
  targets <- tfGeneSets$Target
  cbind(tfGeneSets, corr=c(as.numeric(corrMat[tf,targets] > 0.03) - as.numeric(corrMat[tf,targets] < -0.03)))
})
tfModules_withCorr <- data.frame(rbindlist(tfModules_withCorr_byTF))
head(tfModules_withCorr)
dim(tfModules_withCorr)

rm(corrMat, exprMatrix_filtered, linkList_001, llminW, tfModules.melted, 
   topTFsperTarget.asDf, weightMatrix, linkList_001_byTarget)

# cis-regulatory motif analysis on each of the TF regulons with RcisTarget.

# Motif enrichment analysis & identifying direct targets
if(org=="hg19")
{
  library(RcisTarget.hg19.motifDatabases.20k)
  
  # Motif rankings (genes x motifs)
  data(hg19_500bpUpstream_motifRanking)
  data(hg19_10kbpAroundTss_motifRanking)
  motifRankings <- list()
  motifRankings[["500bp"]] <- hg19_500bpUpstream_motifRanking
  motifRankings[["10kbp"]] <- hg19_10kbpAroundTss_motifRanking
  rm(hg19_500bpUpstream_motifRanking, hg19_10kbpAroundTss_motifRanking)
  
  # Motif annotation (TFs)
  data(hg19_direct_motifAnnotation)
  direct_motifAnnotation <- hg19_direct_motifAnnotation
  data(hg19_inferred_motifAnnotation) # optional
  inferred_motifAnnotation <- hg19_inferred_motifAnnotation
  
  rm(hg19_direct_motifAnnotation, hg19_inferred_motifAnnotation)
}

tfModules_withCorr <- tfModules_withCorr[which(as.character(tfModules_withCorr$TF) %in% allTFs),]
geneInDb <- tfModules_withCorr$Target %in% motifRankings[["500bp"]]@rankings$rn

missingGenes <- sort(unique(tfModules_withCorr[which(!geneInDb),"Target"]))
tfModules_withCorr <- tfModules_withCorr[which(geneInDb),]

# Targets with positive correlation
tfModules_Selected <- tfModules_withCorr[which(tfModules_withCorr$corr==1),]
tfModules_Selected <- cbind(tfModules_Selected, geneSetName=paste(tfModules_Selected$TF, tfModules_Selected$method, sep="_"))
head(tfModules_Selected)

tfModules <- split(tfModules_Selected$Target, tfModules_Selected$geneSetName)

# Keep gene sets with at least 20 genes
tfModules <- tfModules[which(lengths(tfModules)>=20)]

# Add TF to the gene set
tfModules <- setNames(lapply(names(tfModules), function(gsn) {
  tf <- strsplit(gsn, "_")[[1]][1]
  unique(c(tf, tfModules[[gsn]]))
}), names(tfModules))

tfModulesSummary <- t(sapply(strsplit(names(tfModules), "_"), function(x) x[1:2]))
sort(table(tfModulesSummary[,2]))

rm(tfModules_byTF, tfModules_withCorr_byTF)

# Calculate motif enrichment for each TF-module
motifs_AUC <- lapply(motifRankings, function(ranking) calcAUC(tfModules, ranking,
                                                              aucMaxRank=0.01*nrow(ranking@rankings), nCores=1, verbose=TRUE))

motifEnrichment <- lapply(motifs_AUC, function(aucOutput)
{
  tf <- sapply(setNames(strsplit(rownames(aucOutput), "_"), rownames(aucOutput)), function(x) x[[1]])
  addMotifAnnotation(aucOutput, highlightTFs=tf, nesThreshold=3.0, digits=3,
                     motifAnnot_direct=direct_motifAnnotation,
                     motifAnnot_inferred=inferred_motifAnnotation)
})

motifEnrichment <- do.call(rbind, lapply(names(motifEnrichment), function(dbName){
  cbind(motifDb=dbName, motifEnrichment[[dbName]])
}))
head(motifEnrichment)

# TFinDB: direct annotation is ** and inferred annotation is *

# keep only the motifs annotated
motifEnrichment_selfMotifs <- motifEnrichment[which(motifEnrichment$TFinDB != ""),, drop=FALSE]

# Prune targets
motifEnrichment_selfMotifs_wGenes <- lapply(names(motifRankings), function(motifDbName){
  addSignificantGenes(resultsTable=motifEnrichment_selfMotifs[motifDb==motifDbName],
                      geneSets=tfModules,
                      rankings=motifRankings[[motifDbName]],
                      maxRank=5000, method="aprox", nCores=1)
})
motifEnrichment_selfMotifs_wGenes <- rbindlist(motifEnrichment_selfMotifs_wGenes)
motifEnrichment_selfMotifs_wGenes[order(NES,decreasing=TRUE)][1:5,-"enrichedGenes", with=F]

motifEnrichment.asIncidList <- apply(motifEnrichment_selfMotifs_wGenes, 1, function(oneMotifRow) {
  genes <- strsplit(oneMotifRow["enrichedGenes"], ";")[[1]]
  oneMotifRow <- data.frame(rbind(oneMotifRow), stringsAsFactors=FALSE)
  data.frame(oneMotifRow[rep(1, length(genes)),c("NES", "motif", "highlightedTFs", "TFinDB")], genes, stringsAsFactors = FALSE)
})
motifEnrichment.asIncidList <- rbindlist(motifEnrichment.asIncidList)
colnames(motifEnrichment.asIncidList) <- c("NES", "motif", "TF", "annot", "gene")
motifEnrichment.asIncidList <- data.frame(motifEnrichment.asIncidList, stringsAsFactors = FALSE)
head(motifEnrichment.asIncidList)

regulonTargetsInfo <- lapply(split(motifEnrichment.asIncidList, motifEnrichment.asIncidList$TF), function(tfTargets){
  tfTable <- as.data.frame(do.call(rbind, lapply(split(tfTargets, tfTargets$gene), function(enrOneGene){
    directAnnot <- "**" %in% enrOneGene$annot
    enrOneGeneByAnnot <- enrOneGene
    if(directAnnot) enrOneGeneByAnnot <- enrOneGeneByAnnot[which(enrOneGene$annot == "**"),]
    bestMotif <- which.max(enrOneGeneByAnnot$NES)
    
    cbind(TF=unique(enrOneGene$TF), gene=unique(enrOneGene$gene), nMotifs=nrow(enrOneGene), 
          bestMotif=as.character(enrOneGeneByAnnot[bestMotif,"motif"]), NES=as.numeric(enrOneGeneByAnnot[bestMotif,"NES"]), 
          directAnnot=directAnnot)
  })), stringsAsFactors=FALSE)
  tfTable[order(tfTable$NES, decreasing = TRUE),]
})
regulonTargetsInfo <- rbindlist(regulonTargetsInfo)
colnames(regulonTargetsInfo) <- c("TF", "gene", "nMotifs", "bestMotif", "NES", "directAnnot")

linkList <- linkList[which(linkList$weight>=0.001),]
rownames(linkList) <- paste(linkList$TF, linkList$Target,sep="__")

head(regulonTargetsInfo)
regulonTargetsInfo <- cbind(regulonTargetsInfo, Genie3Weight=linkList[paste(regulonTargetsInfo$TF, regulonTargetsInfo$gene,sep="__"),"weight"])
head(regulonTargetsInfo)

regulonTargetsInfo_splitByAnnot <- split(regulonTargetsInfo, regulonTargetsInfo$directAnnot)
regulons <- sapply(split(regulonTargetsInfo_splitByAnnot[["TRUE"]], regulonTargetsInfo_splitByAnnot[["TRUE"]][,"TF"]), function(x) sort(as.character(unlist(x[,"gene"]))))
regulons_extended <- sapply(split(regulonTargetsInfo_splitByAnnot[["FALSE"]],regulonTargetsInfo_splitByAnnot[["FALSE"]][,"TF"]), function(x) unname(x[,"gene"]))
regulons_extended <- sapply(names(regulons_extended), function(tf) sort(unique(c(regulons[[tf]], regulons_extended[[tf]]))))
names(regulons_extended) <- paste(names(regulons_extended), "_extended", sep="")

regulons <- c(regulons, regulons_extended)

incidList <- melt(regulons)
incidMat <- table(incidList[,2], incidList[,1])
dim(incidMat)

table(sapply(names(regulons), function(x) x %in% regulons[[x]]))

# network activity in each cell
regulons <- regulons[order(lengths(regulons), decreasing=TRUE)]
regulons <- regulons[lengths(regulons)>=10]

# Add the TF & rename
regulons <- setNames(lapply(names(regulons), function(tf) sort(unique(c(gsub("_extended", "", tf), regulons[[tf]])))), names(regulons))
names(regulons) <- paste(names(regulons), " (",lengths(regulons), "g)", sep="")

# create gene rankings for each cell
aucellRankings <- AUCell.buildRankings(exprMat, nCores=1, plotStats=TRUE)
abline(v=aucellRankings@nGenesDetected["1%"], col="skyblue3", lwd=5, lty=3)

# calculate AUC
regulonAUC <- AUCell.calcAUC(regulons, aucellRankings, aucMaxRank=aucellRankings@nGenesDetected["1%"], nCores=3)
load("int/3.2_regulonAUC.RData")

variableRegulons <- names(which(apply(getAuc(regulonAUC), 1, sd) > 0))
reguDist <-as.dist(1-cor(t(getAuc(regulonAUC)[variableRegulons,]), method="spear"))
reguClust <- hclust(reguDist, method="ward.D2")
regulonClusters <- setNames(dynamicTreeCut::cutreeDynamic(reguClust, distM=as.matrix(reguDist), verbose = FALSE), reguClust$labels)
regulonOrder <- reguClust$labels[reguClust$order]
regulonOrder <- regulonOrder[order(regulonClusters[regulonOrder], decreasing = TRUE)]
regulonAUC@matrix <- regulonAUC@matrix[regulonOrder,]

regulonAUC_subset <- subset(regulonAUC, unlist(onlyNonDirectExtended(rownames(regulonAUC))))

# PCA-based t-SNE
set.seed(123)
tsneAUC <- Rtsne::Rtsne(t(getAuc(regulonAUC_subset)), initial_dims=10, perplexity=10)
rownames(tsneAUC$Y) <- colnames(regulonAUC_subset)
colnames(tsneAUC$Y) <- c("tsne1", "tsne2")
load("int/3.3_tsneRegulonAUC_PCA.RData")

# Distance-based t-SNE:
corDist <- as.dist(1-cor(getAuc(regulonAUC_subset)))
set.seed(123)
tsneAUC <- Rtsne::Rtsne(corDist, is_distance=TRUE, perplexity=10)
rownames(tsneAUC$Y) <- labels(corDist)
colnames(tsneAUC$Y) <- c("tsne1", "tsne2")

tSNE <- tsneAUC$Y
par(mfrow=c(2,2))

# Number of genes detected:
nGenesPerCell <- apply(exprMat, 2, function(x) sum(x>0))
colorPal <- grDevices::colorRampPalette(c("darkgreen", "yellow","red"))
cellColorNgenes <- setNames(adjustcolor(colorPal(10), alpha=.8)[as.numeric(cut(nGenesPerCell,breaks=10, right=F,include.lowest=T))], names(nGenesPerCell))

plot(tSNE, col=cellColorNgenes[rownames(tSNE)], pch=16, main="nGenes", sub="t-SNE on regulon activity (AUC)")

for(varName in names(colVars))
{
  cellColor <- setNames(colVars[[varName]][cellInfo[,varName]], rownames(cellInfo))
  plot(tSNE, col=cellColor[rownames(tSNE)], pch=16, main=varName, sub="t-SNE on regulon activity (AUC)")
}

load("./data/data.RData")
rm(data.rpkm, data.count, data.genes)

rownames(data.info) <- data.info$Source.Name
data.info <- data.info[rownames(cellInfo), ]

### disease
cellColor <- rep("#B3B6B7", 2166)
cellColor[which(data.info$Characteristics.disease. != "normal")] = "#E74C3C"
cellColor <- setNames(cellColor, rownames(cellInfo))

for(varName in names(colVars))
{
  plot(tSNE, col=cellColor[rownames(tSNE)], pch=16, main=varName, sub="t-SNE on regulon activity (AUC)")
}

### bmi
colors <- rep("", 10)
colors[4] <- "#EAF2F8"
colors[3] <- "#AED6F1"
colors[2] <- "#AED6F1"
colors[10] <- "#7FB3D5"
colors[1] <- "#5499C7"
colors[8] <- "#2980B9"
colors[6] <- "#2471A3"
colors[9] <- "#1A5276"
colors[5] <- "#154360"
colors[7] <- "#17202A"
ind <- as.vector(unique(data.info$Characteristics.individual.))
for(i in 1:10){
  cellColor[which(data.info$Characteristics.individual. == ind[i])] = colors[i]
}
cellColor <- setNames(cellColor, rownames(cellInfo))

for(varName in names(colVars))
{
  plot(tSNE, col=cellColor[rownames(tSNE)], pch=16, main=varName, sub="t-SNE on regulon activity (AUC)")
}

rm(data.info)

Cairo::CairoPDF("scenic/RegulonActivity_AUCtSNE.pdf", width=20, height=5)
par(mfrow=c(1,4))

# tSNE (colored by number of genes detected per cell)
plot(tSNE, col=cellColorNgenes[rownames(tSNE)], pch=16, main="nGenes", sub="t-SNE on regulon activity (AUC)")
plot(tSNE, col=cellColor[rownames(tSNE)], pch=16, main=varName, sub="t-SNE on regulon activity (AUC)")
plot.new(); plot.new()

# Plot module activity, thresholds & assignment:
cells_AUCellThresholds <- plot_aucTsne(tSNE=tSNE, exprMat=exprMat, regulonAUC=regulonAUC, alphaOff=0.1)
dev.off()

# Get cells assigned to each regulon
regulonsCells <- lapply(cells_AUCellThresholds, function(x) x$assignment)

trhAssignment <- sapply(cells_AUCellThresholds, function(x) unname(x$aucThr$selected))
commentsThresholds <- sapply(cells_AUCellThresholds, function(x) unname(x$aucThr$comment))

table2edit <- cbind(regulon=names(trhAssignment), 
                    threshold=trhAssignment, 
                    nCellsAssigned=lengths(regulonsCells)[names(trhAssignment)],
                    AUCellComment=commentsThresholds, 
                    nGenes=gsub("[\\(g\\)]", "", regmatches(names(cells_AUCellThresholds), gregexpr("\\(.*?\\)", names(cells_AUCellThresholds)))),
                    clusteringOrder=1:length(trhAssignment), 
                    clusterGroup=regulonClusters[names(trhAssignment)], 
                    onlyNonDirectExtended=(names(trhAssignment) %in% onlyNonDirectExtended(names(trhAssignment))),
                    personalNotes="")
write.table(table2edit, file="scenic/AUCellThresholds.txt", row.names=F, quote=F, sep="\t")

# Binary network activity
manualThresholds <- read.table("scenic/AUCellThresholds.txt", sep="\t", header=T)
newThresholds <- setNames(manualThresholds[,"threshold"], manualThresholds[,"regulon"])

convert <- newThresholds[which(newThresholds %in% c("Global_k1", "L_k2", "R_k3", "minimumDens"))]
newThresholds[names(convert)] <- sapply(names(convert), function(reg) cells_AUCellThresholds[[reg]]$aucThr$thresholds[convert[reg],"threshold"])
cbind(newThresholds[names(convert)])

secondHigherRegs <- manualThresholds[grep("Bimodal", manualThresholds[,"threshold"]), "regulon"]

# new thresholds
par(mfrow=c(3,5))
for(modName in secondHigherRegs)
{
  densCurve <- density(regulonAUC[modName,], adjust=1)
  inflPoints <- diff(sign(diff(densCurve$y)))
  minimumDens <- densCurve$x[which(inflPoints==2)]
  minimumDens <- minimumDens[1]
  
  AUC.plot(regulonAUC[modName,], gSetName=paste(modName, "module")); abline(v=minimumDens, col="darkorange", lwd=2, lty=2)
  lines(densCurve, col="blue")
  
  newThresholds[modName] <- minimumDens
}

any(is.na(as.numeric(newThresholds)))
newThresholds  <- setNames(as.numeric(newThresholds), names(newThresholds))

Cairo::CairoPDF("scenic/RegulonActivity_AUCtSNE.pdf", width=20, height=5)
par(mfrow=c(1,4))
newAssignment <- plot_aucTsne(exprMat, regulonAUC=regulonAUC, tSNE=tSNE, thresholds=newThresholds)
dev.off()

cells_AUCellThresholds <- newAssignment
rm(newAssignment)

regulonsCells <- lapply(cells_AUCellThresholds, function(x) x$assignment)
length(regulonsCells)

regulonActivity <- reshape2::melt(regulonsCells)
binaryRegulonActivity <- t(table(regulonActivity[,1], regulonActivity[,2]))
class(binaryRegulonActivity) <- "matrix"
sum(rowSums(binaryRegulonActivity)>5)
dim(binaryRegulonActivity)
binaryRegulonActivity[1:10,1:3]
binaryRegulonActivity_nonDupl <- binaryRegulonActivity[which(rownames(binaryRegulonActivity) %in% onlyNonDirectExtended(rownames(binaryRegulonActivity))),]

par(mfrow=c(1,2))
boxplot(rowSums(binaryRegulonActivity_nonDupl), main="nCells per regulon", 
        sub='number of cells \nthat have the regulon active',
        col="darkolivegreen1", border="#001100", lwd=2, frame=FALSE)
boxplot(colSums(binaryRegulonActivity_nonDupl), main="nRegulons per Cell", 
        sub='number of regulons \nactive per cell',
        col="darkolivegreen1", border="#001100", lwd=2, frame=FALSE)

cellInfo <- pData(eset)[,names(colVars), drop=F]
minCells <- ncol(eset) * .01

regulonSelection <- list()

# All regulons.
regulonSelection[["All regulons \n (including duplicated regulons)"]] <- rownames(binaryRegulonActivity)

# Active in > 1% cells
regMinCells <- names(which(rowSums(binaryRegulonActivity_nonDupl) > minCells))
regulonSelection[["Regulons active in more than 1% of cells"]] <- regMinCells

# Correlation across regulons (based on binary cell activity)
reguCor <- cor(t(binaryRegulonActivity_nonDupl[regMinCells,]))
diag(reguCor) <- 0

# Regulons that co-ocurr in similar cells. If a regulon is relevant by itself it will not be shown, also check the regulons ignored.
corrRegs <- names(which(rowSums(abs(reguCor) > 0.30) > 0))
regulonSelection[["Regulons with any other regulon correlated\n with abs(cor)>0.30 \n(and active in at least 1% of cells)"]]  <- corrRegs

missingRegs <- rownames(binaryRegulonActivity_nonDupl)[which(!rownames(binaryRegulonActivity_nonDupl) %in% corrRegs)]
regulonSelection[["Regulons no other regulons correlated\n with abs(cor)>0.30 \n or active in fewer than 1% of cells"]]  <- missingRegs

binaryRegulonOrder <- hclust(as.dist(1-reguCor[corrRegs,corrRegs]))
binaryRegulonOrder <- binaryRegulonOrder$labels[binaryRegulonOrder$order]

for(i in seq_len(length(regulonSelection)))
{
  selRegs <- names(regulonSelection)[i]
  if(length(selRegs)>=1)
  {
    binaryMat <- binaryRegulonActivity[regulonSelection[[selRegs]],,drop=FALSE]
    NMF::aheatmap(binaryMat, scale="none", revC=TRUE, main=selRegs,
                  annCol=cellInfo[colnames(binaryMat),, drop=FALSE],
                  annColor=colVars,
                  color = c("white", "black"),
                  filename=paste0("scenic/binaryRegulonActivity_heatmap_",i,".pdf"))
  }
}

# GRN-based cell states
cellInfo <- pData(eset)[, names(colVars), drop=F]
tBinaryAct <- t(binaryRegulonActivity_nonDupl)
tBinaryAct <- t(binaryRegulonActivity_nonDupl[binaryRegulonOrder, ])

# PCA based t-SNE
set.seed(123)
tBinaryAct_jitter <- jitter(tBinaryAct, factor=1)
tsneBinaryActivity_PCA <- Rtsne(tBinaryAct_jitter, initial_dims=5, perplexity=30)
rownames(tsneBinaryActivity_PCA$Y) <- rownames(tBinaryAct_jitter)
colnames(tsneBinaryActivity_PCA$Y) <- c("tsne2", "tsne1")
tsneBinaryActivity_PCA$Y <- tsneBinaryActivity_PCA$Y[,c("tsne1", "tsne2")]

# PCA based t-SNE
set.seed(123)
tBinaryAct_jitter <- jitter(tBinaryAct, factor=1)
tsneBinaryActivity_PCA <- Rtsne(tBinaryAct_jitter, initial_dims=50, perplexity=30)
rownames(tsneBinaryActivity_PCA$Y) <- rownames(tBinaryAct_jitter)
colnames(tsneBinaryActivity_PCA$Y) <- c("tsne2", "tsne1")
tsneBinaryActivity_PCA$Y <- tsneBinaryActivity_PCA$Y[,c("tsne1", "tsne2")]

# Distance-based t-SNE
corDist <- as.dist(1-cor(t(tBinaryAct)))
set.seed(123)
tsneBinaryActivity_Dist <- Rtsne(corDist, is_distance=TRUE, perplexity=30)
rownames(tsneBinaryActivity_Dist$Y) <- labels(corDist)
colnames(tsneBinaryActivity_Dist$Y) <- c("tsne1", "tsne2")

tSNEs_binary <- list()
tSNEs_binary[["Dist"]] <- tsneBinaryActivity_Dist$Y
tSNEs_binary[["5PC"]] <- tsneBinaryActivity_PCA$Y
tSNEs_binary[["50PC"]] <- tsneBinaryActivity_PCA$Y

tsneName <- names(tSNEs_binary)[3]
tSNE_binary <- tSNEs_binary[[tsneName]]

# Choose a t-SNE
tSNE_binary <- tsneBinaryActivity_PCA$Y
tSNEname <- "tsneBinaryActivity_50PC"

regOrder <- binaryRegulonOrder[which(binaryRegulonOrder %in% rownames(tBinaryAct))]
Cairo::CairoPDF(paste0("scenic/",tSNEname,"_BinaryRegulons.pdf"), width=20, height=15)
par(mfrow=c(4,6))
cells_trhAssignment <- plot_aucTsne(tSNE=tSNE_binary, exprMat=exprMat,
                                    regulonAUC=t(tBinaryAct)[binaryRegulonOrder,], cex=1.5, plots="binary", thresholds=0)
dev.off()

regOrder<- binaryRegulonOrder[which(binaryRegulonOrder %in% rownames(regulonAUC))]
Cairo::CairoPDF(paste0("scenic/AUCRegulons.pdf"), width=8, height=5, family = "arial")
par(mfrow=c(2,3), mar=c(1.5,1.5,1.5,1.5))
cells_trhAssignment <- plot_aucTsne(tSNE=tSNE_binary, exprMat=exprMat, 
                                    regulonAUC=regulonAUC[regOrder,], cex=1, plots="AUC", thresholds=0)
dev.off()

Cairo::CairoPDF(paste0("scenic/",tSNEname,"_allPlots.pdf"), width=20, height=5)
par(mfrow=c(1,4))
cells_trhAssignment <- plot_aucTsne(tSNE=tSNE_binary, exprMat=exprMat, 
                                    regulonAUC=regulonAUC[regOrder,],
                                    alphaOff=0.1, thresholds=cells_AUCellThresholds[regOrder])
dev.off()
