library(AnnotationDbi)
library(biomaRt)
library(GenomicFeatures) 
library(rmarkdown)
library(data.table)
library(RMariaDB)
library(org.Hs.eg.db)
library("TxDb.Hsapiens.UCSC.hg38.knownGene")
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
library(pheatmap)

# set your working directory
setwd("WORKDIRECTORY")

# get all the ratiofiles from the folder using the fast read function
filenames <- list.files('Ratio_files/')
ratio_files_list <- lapply(filenames, function(x){
  fread(paste0("Ratio_files/",x))
})
# set the correct names so it is the same as in the rest of the paper .
names(ratio_files_list) <-
c("T2 Organoid", "T1 Normal tissue", "T1 Organoid", "T1 Tumor tissue", "T16 Normal", "T16 Organoid", "T16 Tumor tissue", "T11 Organoid", "T12 Normal",
  "T12 Organoid", "T12 Tumor tissue", "T13 Normal", "T13 Organoid", "T14 Normal","T14 Organoid", "T14 Tumor tissue", "T17 Normal", "T17 Organoid", "T18 Normal",
  "T18 Organoid", "T18 Tumor tissue", "T19 Normal", "T19 Organoid", "T20 Normal", "T20 Organoid", "T20 Tumor tissue", "T21 Normal", "T21 Organoid",
  "T22 Normal", "T22 Organoid", "T35 Normal", "T35 Organoid", "T34 Normal", "T34 Organoid", "T23 Normal", "T23 Organoid", "T24 Normal",
  "T24 Organoid", "T26 Normal", "T26 Organoid", "T32 Normal", "T32 Organoid", "T36 Normal", "T36 Organoid", "T37 Organoid", "T40 Normal",
  "T40 Organoid", "T41 Organoid", "T42 Normal", "T42 Organoid", "T43 Normal", "T43 Organoid", "T45 Organoid", "T46 Organoid", "T49 Normal",
  "T49 Organoid", "T25 Normal", "T25 Organoid", "TXX Organoid", "T27 Organoid", "T55 Organoid", "T55 Tumor tissue")

# devide the copynumbers to genomic ranges. 
gr <- lapply(names(ratio_files_list), function(x){
  GRanges(
    seqnames = paste0("chr",ratio_files_list[[x]]$Chromosome),
    ranges = unlist(lapply(ratio_files_list[[x]]$Gene, function(x){strsplit(x, ":")[[1]][2]})),
    strand = "+",
    ratio = ratio_files_list[[x]]$Ratio,
    medianratio = ratio_files_list[[x]]$MedianRatio,
    copynumber = ratio_files_list[[x]]$CopyNumber
  )})
# give the correct names
names(gr) <- names(ratio_files_list)

# link the genomic regions to geneIDs from a reference genome
# Making a reference with both genomic location and gene symbol
txid <- keys(txdb, "TXID")
df <- select(txdb, txid, "GENEID", "TXID")
genes <- as.data.frame(genes(txdb))
genes <- genes[order(genes$seqnames),]
gene_symbol <- select(org.Hs.eg.db, df$GENEID, c("SYMBOL", "GENENAME"))
gene_symbol <- gene_symbol[-which(duplicated(gene_symbol$ENTREZID)),]
gene_symbol <- gene_symbol[-which(is.na(gene_symbol$ENTREZID)),]
rownames(gene_symbol) <- gene_symbol$ENTREZID
genetable <- genes
genetable$gene_symbol <- gene_symbol[genetable$gene_id,2]


# Calculate the average copynumber for specific genes.
# Only scores a gene if the entire bin is locaten within the same gene.
ratio_tables_genes <- lapply(names(gr), function(x){
  ratio_counts <- as.data.frame(gr[[x]])
  ratio_counts$genesymbol <- NA
  for (row in c(1:dim(ratio_counts)[[1]])) {
    chr <- ratio_counts[row, 1]
    start <- ratio_counts[row,2]
    end <- ratio_counts[row,3]
    temp <- genetable[which(genetable$seqnames == as.character(chr)),]
    rowselect <- intersect(which(temp$start < start),which(temp$end > end))
    if (!is.null(rowselect)) {
      if(length(rowselect) > 1){
        symbol <- temp$gene_symbol[rowselect[1]]
        ratio_counts[row, 9] <- symbol
      }
      if(length(rowselect) == 1){
        symbol <- temp$gene_symbol[rowselect]
        ratio_counts[row, 9] <- symbol
      }}
  }
  return(ratio_counts)
})
# give the correct names
names(ratio_tables_genes) <- names(ratio_files_list)


# select for the most variable genes in literature (Cbio portal head and neck dataset), based on gene symbol 
CNA_genes <- read.table("CNA_Genes.txt", header = F, sep = "\t", skip = 1)
CNA_genes$V7 <- unlist(lapply(CNA_genes$V7, function(x){strsplit(x, "%")[[1]][1]}))
CNA_genes$V7 <- as.numeric(CNA_genes$V7)
CNA_genes <- CNA_genes[order(CNA_genes[,7], decreasing = T),]


Genes_of_interest <- c(CNA_genes$V1[-grep("^RN", CNA_genes$V1)][c(1:70)], "EGFR", "PIK3CA")

# calculate the average copynmeber per gene in the genelist. 
genecopytable <- lapply(Genes_of_interest, function(x){
  genevalues <- unlist(lapply(names(ratio_tables_genes), function(y){
    temp <- ratio_tables_genes[[y]][which(ratio_tables_genes[[y]]$genesymbol == x),]
    value <- mean(temp$ratio) * 2
    return(value)
  }))
  return(genevalues)
})

# make the table ready for use in a heatmap
names(genecopytable) <- Genes_of_interest
genecopytable <- as.data.frame(genecopytable)
rownames(genecopytable) <- names(ratio_tables_genes)
genecopytable <- t(genecopytable)
genecopytable <- genecopytable[-unlist(which(is.na(genecopytable[,1]))),]

# adding chromosome positions to gene names to show stretches of genes
annotationtable <- genetable[which(genetable$gene_symbol %in% list_to_change),]
rownames(annotationtable) <- annotationtable$gene_symbol
annotationtable$newname <- unlist(lapply(rownames(annotationtable), function(x){
  gene_symbol <- x
  chromosome <- annotationtable[x,"seqnames"]
  start <- annotationtable[x, "start"]
  returnvalue <- paste0(gene_symbol, "_", chromosome, "_", start)
  return(returnvalue)
}))
rownames(genecopytable) <- annotationtable[rownames(genecopytable), "newname"]

# Plot the copynumber changes in a heatmap
pheatmap(genecopytable[,grep("Tumor", colnames(genecopytable))],
         color = c("purple","violet", "#F8F8FF","#F8F8FF", "#FFF5F0","#FEE0D2",
                   "#FCBBA1", "#FC9272", "#FB6A4A" ,"#EF3B2C", "#CB181D",
                   "#A50F15", "#67000D", "#67000D", "#67000D", "#67000D", "#67000D", "#67000D", "#67000D", "#67000D",
                   "#67000D", "#67000D", "#67000D", "#67000D", "#67000D", "#67000D", "#67000D", "#67000D", "#67000D", "#67000D", "#67000D",
                   "#67000D", "#67000D", "#67000D", "#67000D", "#67000D", "#67000D", "#67000D", "#67000D", "#67000D", "#67000D", "#67000D"),
         main = "mean coverage of genes known to have copy number changes",
         cluster_rows = F,
         cluster_cols = T,
         border_color = "white",
         treeheight_col = 0,
         treeheight_row = 0, 
         gaps_col = 23)

