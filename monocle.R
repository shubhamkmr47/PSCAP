library(monocle)
library(scran)
library(ReactomePA)
library(clusterProfiler)
library(erer)
library(reshape2)

# load data
load("./data/data.RData")
rm(data.rpkm)

# consider only good quality cells
cells <- data.info[data.info$Characteristics.single.cell.well.quality. == "OK", ]

cell_type <- as.vector(cells$Characteristics.cell.type.)
cell_type[cell_type == "alpha cell" | cell_type == "beta cell" | cell_type == "gamma cell" | cell_type == "delta cell" | cell_type == "epsilon cell"] <- "endocrine cells"
cell_type[cell_type == "acinar cell" | cell_type == "ductal cell"] <- "exocrine cells"

cells <- cbind(cells, cell_type)
cells <- cells[order(cells$Characteristics.disease.), ]
cells <- cells[order(cells$Characteristics.cell.type.), ]

table(cells$Characteristics.cell.type.)

ids <- as.vector(cells$Source.Name)
data.count <- data.count[, ids]
rownames(cells) <- as.vector(cells$Source.Name)

row <- rownames(data.count)
data.count <- apply(data.count, 2, as.numeric)
rownames(data.count) <- row

rm(ids, cell_type, row)

all.equal(as.vector(data.genes$G), rownames(data.count))
dim(data.count)

# remove ERCC genes
data.count <- data.count[-c(26179:26271), ]
data.genes <- data.genes[-c(26179:26271), ]

# remove duplicates
gene.var <- apply(log(data.count + 1, 2), 1, function(x) iqr(x[x > 0]))
data.genes <- cbind(data.genes, iqr_count = gene.var)
data.genes <- data.genes[order(data.genes$iqr_count, decreasing = TRUE),]
data.genes <- data.genes[!duplicated(data.genes$name), ]
data.genes <- data.genes[order(as.vector(data.genes$G)), ]

data.count <- data.count[as.vector(data.genes$G), ]
data.count <- data.count[order(rownames(data.count)), ]

all.equal(as.vector(data.genes$G), rownames(data.count))

dim(data.count)

# filter out genes with zero expression across 25% of cells
keep <- rowSums(data.count > 0) >= round(dim(cells)[1]/4)
sum(keep)
data.count <- data.count[keep,] 

# filter out genes with very low average nonzero expression across all cells
keep <- rowMeans(data.count) >= 5
sum(keep)
data.count <- data.count[keep,]

rm(keep, gene.var)

rownames(data.genes) <- as.vector(data.genes$G)
data.genes <- data.genes[rownames(data.count), ]
colnames(data.genes)[2] <- "gene_short_name"

### Monocle
pd <- new("AnnotatedDataFrame", data = cells)
fd <- new("AnnotatedDataFrame", data = data.genes)

# create celldataset (negative binomial distribution with fixed variance)
HSMM <- newCellDataSet(as(data.count, "sparseMatrix"),
                       phenoData = pd,
                       featureData = fd,
                       lowerDetectionLimit = 1,
                       expressionFamily = negbinomial.size())

# size factors to normalize for differences in mRNA recovered across cells
# dispersion values will help us perform differential expression analysis.
HSMM <- estimateSizeFactors(HSMM)
HSMM <- estimateDispersions(HSMM)
HSMM <- detectGenes(HSMM, min_expr = 1)

round(dim(cells)[1]/4)
expressed_genes <- row.names(fData(HSMM))
print(head(pData(HSMM)))

# vVerifying normalization
L <- log(exprs(HSMM[expressed_genes,]))
melted_dens_df <- melt(Matrix::t(scale(Matrix::t(L))))

qplot(value, geom = "density", data = melted_dens_df) +
  stat_function(fun = dnorm, size = 0.5, color = 'red') +
  xlab("Standardized log(FPKM)") +
  ylab("Density")
rm(L, melted_dens_df, pd, fd, data.count)

# save.image("./data/monocle.RData")

###

load("./data/monocle.RData")

# alpha cell, beta cell, gamma cell, delta cell, acinar cell, ductal cell
cell_type <- "ductal cell"
cells <- cells[cells$Characteristics.cell.type. == cell_type, ]

# directory to save results
fil <- paste("monocle/", cell_type, "", sep = "/")

HSMM <- HSMM[, as.vector(cells$Source.Name)]

diff_test_res <- differentialGeneTest(HSMM[expressed_genes, ],
                fullModelFormulaStr = "~Characteristics.individual.", cores = 3)

write.csv(diff_test_res, paste(fil, "deg.csv", sep = ""),
          row.names = FALSE, quote = FALSE)

# alpha beta gamma ductal
ordering_genes <- diff_test_res[order(diff_test_res$qval, decreasing = FALSE)[1:150], ]

# delta acinar
ordering_genes <- diff_test_res[order(diff_test_res$qval, decreasing = FALSE)[1:200], ]

ordering_genes <- row.names(ordering_genes[ordering_genes$qval <= 1e-02, ])

HSMM <- setOrderingFilter(HSMM, ordering_genes)
plot_ordering_genes(HSMM)

#reduce data dimensionality
HSMM <- reduceDimension(HSMM, max_components = 2, method = 'DDRTree')

# order cells along the trajectory
# alpha beta gamma delta
HSMM <- orderCells(HSMM)

# acinar
HSMM <- orderCells(HSMM)
HSMM <- orderCells(HSMM, reverse = TRUE)

# ductal
HSMM <- orderCells(HSMM)
HSMM <- orderCells(HSMM, root_state = 3)

png(filename= paste(fil, "state.png", sep = ""), width = 800, height = 470)
plot_cell_trajectory(HSMM, color_by = "State")
dev.off()

png(filename= paste(fil, "pseudotime.png", sep = ""), width = 800, height = 470)
plot_cell_trajectory(HSMM, color_by = 'Pseudotime');
dev.off()

png(filename= paste(fil, "individual.png", sep = ""), width = 800, height = 470)
plot_cell_trajectory(HSMM, color_by = "Characteristics.individual.")
dev.off()

png(filename= paste(fil, "bmi.png", sep = ""), width = 800, height = 470)
plot_cell_trajectory(HSMM, color_by = "Characteristics.body.mass.index.")
dev.off()

png(filename= paste(fil, "age.png", sep = ""), width = 800, height = 470)
plot_cell_trajectory(HSMM, color_by = "Characteristics.age.")
dev.off()

### Beam
set.seed(123)

# branch point
bp = 1

BEAM_res <- BEAM(HSMM, branch_point = bp, cores = 3)
BEAM_res <- BEAM_res[order(BEAM_res$qval), ]

beam <- BEAM_res[order(BEAM_res$qval, decreasing = FALSE)[1:250], ]
beam <- as.vector(beam[beam$qval <= 1e-03, ]$G)

write.csv(BEAM_res, paste(fil, "deg_beam", bp, ".csv", sep = ""), row.names = FALSE, quote = FALSE)

t <- 5

Cairo::CairoPDF(paste(fil, "psuedotime heatmap_", bp, ".pdf", sep = ""), width=12, height=20)
my_branched_heatmap <- plot_genes_branched_heatmap(HSMM[beam, ],
       branch_point = bp, num_clusters = t, cores = 3,
       use_gene_short_name = T, show_rownames = T, 
       return_heatmap = T)
dev.off()

# plot important genes
g <- c("INS", "JUNB", "PTEN", "EGR1", "FBXO2", "TFF3")
g <- as.vector(diff_test_res[diff_test_res$gene_short_name %in% g, ]$G)
i = 1
b = 1
cl <- 3
while(i < length(g)){
  cds <- HSMM[g[i:min(length(g), i+5)],]
  p <- plot_genes_branched_pseudotime(cds, branch_point = bp, ncol = cl, color_by = "State", min_expr = 10)
  plot(p)
  i = i+6
  b = b+1
}
