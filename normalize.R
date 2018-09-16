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

# QC to compare the level of dropout in endogeneous genes to ERCC spike ins in 
# raw data designate which rows contain spikeins instead of genes (for HVG analysis)
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

plotExplanatoryVariables(eset, exprs_values = "counts", variables =
     c("Characteristics.age.", "Characteristics.body.mass.index."))

plot(sizeFactors(eset), colSums(counts(eset))/1e6, log="xy",
     ylab="Library Size (Total Counts in Millions)", xlab="Pooled Size Factor Estimate",
     main="Normalization factor versus library size")

detectionRate <- apply(counts(eset), 2, function(x) sum(x > 0) / length(x))
hist(detectionRate)

plot(detectionRate, df$PC1, pch=20, col="grey", ylab="PC 1")

rm(detectionRate, df); gc()  #clean up workspace

# save.image("./data/normalize.RData")
