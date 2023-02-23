library(ggplot2)
library(pheatmap)
library(data.table)
library(maftools)
library(readxl)
library(NMF)
library(BiocManager)
setwd("/Users/gijsvanson/OneDrive - Prinses Maxima Centrum/Else_Rosie/")

# set the variants and their colors that are (not) used in the maf-tools analysis
syn <- c("synonymous_variant", "intron_variant")
df <- fread("maf_rightnames_variants_classified.maf")
vc <- names(table(df$Variant_Classification))
nonSyn <- setdiff(vc,syn)
colors <- rainbow(length(nonSyn))
names(colors) <- nonSyn

# read the maf-file
maf <- read.maf("maf_rightnames_variants_classified.maf", vc_nonSyn = nonSyn)

# filter out germline or normal samples.
maf <- maftools::filterMaf(maf, tsb = grep("Germ", maf@variants.per.sample$Tumor_Sample_Barcode, value = T, ignore.case = T))
maf <- maftools::filterMaf(maf, tsb = grep("N1", maf@variants.per.sample$Tumor_Sample_Barcode, value = T, ignore.case = T))

#read the metadata and set it in the clinical.data slot
metadata <- read_xlsx("WES list for Gijs.xlsx")
metadata <- as.data.frame(metadata)
rownames(metadata) <- metadata$Samplenames
metadata_filtered <- metadata[which(rownames(metadata) %in% unique(maf@data$Tumor_Sample_Barcode)),]
maf@clinical.data$HPV <- metadata_filtered[maf@clinical.data$Tumor_Sample_Barcode,"HPV"]
maf@clinical.data$SEX <- metadata_filtered[maf@clinical.data$Tumor_Sample_Barcode,"Sex"]
maf@clinical.data$`Fanconi Anemia` <- metadata_filtered[maf@clinical.data$Tumor_Sample_Barcode,"Fanconi Anemia"]
maf@clinical.data$Age <- unlist(lapply(metadata_filtered[maf@clinical.data$Tumor_Sample_Barcode,"Year of birth"], function(x){2022-x}))
maf@clinical.data$Treatment <- metadata_filtered[maf@clinical.data$Tumor_Sample_Barcode,"Previous treatment"]
maf@clinical.data$Type <- metadata_filtered[maf@clinical.data$Tumor_Sample_Barcode,"tumor type"]
maf@clinical.data$Location <- metadata_filtered[maf@clinical.data$Tumor_Sample_Barcode,"Anatomical Location"]
maf@clinical.data$Metestasis_status <- metadata_filtered[maf@clinical.data$Tumor_Sample_Barcode,"Primary/recurrence/met status"]
maf@clinical.data$`Fanconi Anemia`[which(maf@clinical.data$`Fanconi Anemia` == T)] <- "Yes"
maf_filtered@clinical.data$Location[which(maf_filtered@clinical.data$Location == "Cavum nasi")] <- "Nasal cavity"
maf@clinical.data$Fanconi <- maf@clinical.data$`Fanconi Anemia`
# Set the colors for the different metadata columns
HPVcols <- c("lightgray", "red")
Fanconi_cols <- c("lightgray", "blue")
Typecols <- c("aquamarine4", "coral3")
Locationcols <-c("green", "blue", "red", "yellow", 
                 "black", "magenta","green")
Treatmentcols <- c( "lightgray", "yellow","red")
names(HPVcols) <-unique(maf@clinical.data$HPV)
names(Fanconi_cols)<- unique(maf@clinical.data$Fanconi)
names(Typecols) <- unique(maf@clinical.data$Type)
names(Locationcols) <- unique(maf_filtered@clinical.data$Location)
names(Treatmentcols) <- unique(maf@clinical.data$Treatment)

# Filter out variants that do not map to a gene
maf_filtered <-filterMaf(maf, genes = '""')

# filter out samples that are not confirmed to be tumor and select for organoid only.
maf_filtered <- filterMaf(maf_filtered, tsb = grep("Tumor tissue", maf@variants.per.sample$Tumor_Sample_Barcode, value = T, ignore.case = T))
maf_filtered <- filterMaf(maf_filtered, tsb = grep("T18 ", maf@variants.per.sample$Tumor_Sample_Barcode, value = T, ignore.case = T))
maf_filtered <- filterMaf(maf_filtered, tsb = grep("T22 ", maf@variants.per.sample$Tumor_Sample_Barcode, value = T, ignore.case = T))
maf_filtered <- filterMaf(maf_filtered, tsb = grep("T25", maf@variants.per.sample$Tumor_Sample_Barcode, value = T, ignore.case = T))
maf_filtered <- filterMaf(maf_filtered, tsb = grep("T26", maf@variants.per.sample$Tumor_Sample_Barcode, value = T, ignore.case = T))
maf_filtered <- filterMaf(maf_filtered, tsb = grep("T29", maf@variants.per.sample$Tumor_Sample_Barcode, value = T, ignore.case = T))
maf_filtered <- filterMaf(maf_filtered, tsb = grep("T34", maf@variants.per.sample$Tumor_Sample_Barcode, value = T, ignore.case = T))
maf_filtered <- filterMaf(maf_filtered, tsb = grep("T40", maf@variants.per.sample$Tumor_Sample_Barcode, value = T, ignore.case = T))
maf_filtered <- filterMaf(maf_filtered, tsb = grep("T43 ", maf@variants.per.sample$Tumor_Sample_Barcode, value = T, ignore.case = T))
maf_filtered <- filterMaf(maf_filtered, tsb = grep("Normal", maf@variants.per.sample$Tumor_Sample_Barcode, value = T, ignore.case = T))

# select significant genes
maf_filtered <- subsetMaf(maf_filtered, genes=c("TP53", "PIK3CA", "FAT1", "HERC1", "HRAS", "CASP8", "CDKN2A",
                                                "AJUBA", "NOTCH1", "KMT2D", "NSD1", "HLA-A", "TGFBR2", "APOB", "FBXW7",
                                                "RB1", "PIK3R1", "TRAF3", "NFE2L2", "CLU3", "PTEN", "LAMA2"))
# make the plots and save them to a pdf-file
pdf("oncoplots.pdf", width = 20, height = 8)
plotmafSummary(maf_filtered, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE, color = colors)
oncoplot(maf_filtered, top = 30, showTumorSampleBarcodes = T, SampleNamefontSize = 1, legend_height = 4, legendFontSize = 1, barcode_mar = 8, removeNonMutated = F,
         clinicalFeatures = c( "Type","HPV", "Fanconi", "Treatment", "Location"),annotationColor = list("HPV"= HPVcols,
                                                                                                               "Fanconi" = Fanconi_cols,
                                                                                                               "Type" = Typecols,
                                                                                                               "Treatment" = Treatmentcols,
                                                                                                               "Location" = Locationcols),sortByAnnotation = T)
dev.off()

# Create a matrix to look at mutational signatures
tes_mat <- trinucleotideMatrix(maf_filtered, ref_genome = 'BSgenome.Hsapiens.UCSC.hg38', prefix = "chr", add = T, useSyn = T)
tes_sign <- estimateSignatures(tes_mat)
test_sign3 <- extractSignatures(tes_mat, n = 5, pConstant = 0.01)


# plot the signatures
p1 <- plotSignatures(test_sign3, contributions = T, show_barcodes = T, sig_db = T, show_title = T)
p2 <- plotSignatures(test_sign3)
dev.off()


#lollipop for fanconi figure
lollipopPlot(maf_filtered, gene = "FANCC", refSeqID = "NM_000136", labelPos = 23, showLegend = F, printCount = F, showMutationRate = F, axisTextSize = c(1,1),
             defaultYaxis = F, labPosSize = 3, repel = F, pointSize = 3, domainLabelSize = 3, labPosAngle = 45)



# make a maf-file for the oncopanel data
oncopanel <- data.frame("Hugo_Symbol" = c("PIK3CA", "TP53", "CCND1", "CDKN2A", 
                                           "PIK3CA", "TP53", "EGFR", 
                                          "TP53", "CCND1", "TP53", "CCND1",
                                          "NRAS", "PIK3CA"), 
                        "Entrez_Gene_Id" = NA, "Center" = NA, "NCBI_Build" = 38, "Chromosome" = c(1,2,3,4,5,6,7,8,9,10,11,12,13),
                        "Start_Position" = c(1,2,3,4,5,6,7,8,9,10,11,12,13), "End_Position" = c(13,12,11,10,9,8,7,6,5,4,3,2,1),
                        "Strand" = NA, "Variant_Classification" = c("Missense mutation","Missense mutation", "Amplification", "Homozygous deletion", 
                                                                    "Missense mutation", "Missense mutation","Amplification", 
                                                                    "Missense mutation", "Amplification", "Missense mutation", "Amplification",
                                                                    "Missense mutation", "Missense mutation"),
                        "Variant_Type" = "SNP", "Reference_Allele" = "G", "Tumor_Seq_Allele1" = "A", "Tumor_Seq_Allele2" = "C", 
                        "Tumor_Sample_Barcode" = c("T15 Organoid", "T15 Organoid", "T15 Organoid", "T15 Organoid",
                                                  "T68 Organoid", "T68 Organoid", "T68 Organoid",
                                                   "T69 Organoid", "T69 Organoid", "T73 Organoid", "T73 Organoid",
                                                  "T66 Organoid", "T66 Organoid"),
                        "Protein_Change" = NA, "i_TumorVAF_WU" = 1, "i_transcript_name"=NA)

# write the file to disk
fwrite(oncopanel, "oncopanel.maf", row.names = F, sep = "\t")

# maftools on the oncopanel data. 
syn <- c("synonymous_variant", "intron_variant")
df <- fread("oncopanel.maf")
head(df)
vc <- names(table(df$Variant_Classification))
nonSyn <- setdiff(vc,syn)
colors <- rainbow(length(nonSyn))
names(colors) <- nonSyn
maf <- read.maf("oncopanel.maf", vc_nonSyn = nonSyn)
maf@data
oncoplot(maf, top = 30, showTumorSampleBarcodes = T, SampleNamefontSize = 1, legend_height = 4, legendFontSize = 1, barcode_mar = 8, removeNonMutated = F)

# add metadata 
metadata_onco <- data.frame("Sample_name" = c("T15 Organoid", "T69 Organoid", "T73 Organoid", "T68 Organoid", "T66 Organoid"),
                            "Sex" = c("F", "F", "M", "F", "M"), 
                            "Year of birth" = c(1956, 1970, 1939, 1952, 1956),
                            "Location" = c("Lymph node", "Oral cavity", "salivary gland", "Larynx", "Larynx"),
                            "Type" = c("SCC", "SCC", "AC", "SCC", "SCC"),
                            "Sample type" =c("Organoid", "Organoid", "Organoid", "Organoid", "Organoid"),
                            "Fanconi" = c("No","No","No","No","No"),      
                            "HPV" = c("negative","negative","negative","negative","negative"),
                            "Treatment" = c("untreated","untreated","untreated","untreated","untreated"))
rownames(metadata_onco) <- metadata_onco$Sample_name
metadata_onco <- as.data.table(metadata_onco)
maf@clinical.data <- metadata_onco
maf@clinical.data
maf@clinical.data$Tumor_Sample_Barcode <- maf@clinical.data$Sample_name

oncoplot(maf, top = 30, showTumorSampleBarcodes = T, SampleNamefontSize = 1, legend_height = 4, legendFontSize = 1, barcode_mar = 8, removeNonMutated = F,
clinicalFeatures = c( "Type","HPV", "Fanconi", "Treatment", "Location"),annotationColor = list("HPV"= HPVcols2,
                                                                                               "Fanconi" = Fanconi_cols2,
                                                                                               "Type" = Typecols,
                                                                                               "Treatment" = Treatmentcols2,
                                                                                               "Location" = Locationcols2))

Locationcols2 <- c("red", "black", "magenta", "green")
names(Locationcols2) <- c("Larynx",  "Oral cavity", "salivary gland", "Lymph node")
HPVcols2 <- c("lightgray")
names(HPVcols2) <- "negative"
Fanconi_cols2 <- "lightgray"
names(Fanconi_cols2) <- "No"
Treatmentcols2 <- "lightgray"
names(Treatmentcols2) <- "untreated"
