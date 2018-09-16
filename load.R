library(biomaRt)

# load information file
data.info <- read.csv("./data/E-MTAB-5061.sdrf.txt", sep = '\t', stringsAsFactors = FALSE)
data.info <- data.info[, c(1, 6:12, 15)]

table(data.info$Characteristics.cell.type.)
x <- barplot(table(data.info$Characteristics.cell.type.), xaxt="n", main = "Histogram of celltypes")
labs <- paste(names(table(data.info$Characteristics.cell.type.)))
text(cex=0.7, x=x, y=-1.25, labs, xpd=TRUE, srt=90)

cells <- data.info[data.info$Characteristics.single.cell.well.quality. == "OK", ]
tab <- table(cells$Characteristics.individual., cells$Characteristics.cell.type.)
tab <- as.data.frame.matrix(tab)
t1 <- rowSums(tab)
t2 <- colSums(tab)
tab <- cbind(tab, total = t1)
tab <- rbind(tab, total = t2)
tab <- t(tab)
write.csv(tab, "./results/celltype_vs_individuals.csv")

rm(x, t1, t2, tab, labs)

data <- readLines("./data/pancreas_refseq_rpkms_counts_3514sc.txt")
cols <- data[1]
cols <- strsplit(cols, split = "\t")[[1]]
cols <- cols[-1]
cols <- c(cols, cols)
cols <- c("gene_name", "gene_id", cols)

data <- read.csv("./data/pancreas_refseq_rpkms_counts_3514sc.txt", 
                 sep = "\t", header = FALSE, stringsAsFactors = FALSE)
colnames(data) <- cols
data <- data[-c(1), ]
rownames(data) <- (1:nrow(data))

data.rpkm <- data[, c(1:2, 3:3516)]
data.count <- data[, c(1:2, 3517:7030)]

rm(data)

G <- c()
for(i in 1:dim(data.rpkm)[1]){
  G <- c(G, paste("G", i, sep = ""))
}

data.genes <- cbind(G, name = as.vector(data.rpkm$gene_name), 
                    id = as.vector(data.rpkm$gene_id))
data.genes <- data.frame(data.genes)

# genes annotations
ensembl  <- useMart(host="www.ensembl.org", "ENSEMBL_MART_ENSEMBL",
                    dataset="hsapiens_gene_ensembl")
genemap <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol", "entrezgene"),
                 filters = "hgnc_symbol",
                 values = as.vector(data.genes$name),
                 mart = ensembl)

genemap <- genemap[match(as.vector(data.genes$name), genemap$hgnc_symbol), ]
data.genes <- cbind(data.genes, genemap)
rownames(data.genes) <- 1:dim(data.genes)[1]

data.rpkm$gene_name = NULL
data.rpkm$gene_id = NULL
data.rpkm <- data.frame(data.rpkm)
rownames(data.rpkm) <- G

data.count$gene_name = NULL
data.count$gene_id = NULL
data.count <- data.frame(data.count)
rownames(data.count) <- G

rm(i, G, cols, ensembl, genemap)

# save.image("data/data.RData")
