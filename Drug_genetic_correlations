# clinical significance correlations

genecopynumbers_all <- ratio_tables_genes
genecopytable_all <- lapply(unique(genecopynumbers_all$`T2 Organoid`$genesymbol), function(x){
  genevalues <- unlist(lapply(names(genecopynumbers_all), function(y){
    temp <- ratio_tables_genes[[y]][which(ratio_tables_genes[[y]]$genesymbol == x),]
    value <- mean(temp$ratio) * 2
    return(value)
  }))
  return(genevalues)
})
names(genecopytable_all) <- unique(genecopynumbers_all$`T2 Organoid`$genesymbol)
genecopytable_all <- as.data.frame(genecopytable_all)
rownames(genecopytable_all) <- names(genecopynumbers_all)
genecopytable_all <- t(genecopytable_all)
maf_4_correlations <- maf_filtered
samples <- unique(maf_4_correlations@data$Tumor_Sample_Barcode)

zscores <- read_xlsx("TT_Zscores_Forgijs_131220222.xlsx", sheet = 2)
zscores <- as.data.frame(zscores)
rownames(zscores) <- zscores$Compound
zscores <- zscores[,-1]
zscores <- zscores[,which(colnames(zscores) %in% as.character(samples))]
x <- rownames(zscores)[1]
responders <- lapply(rownames(zscores), function(x){
  return(colnames(zscores)[which(zscores[grep(x, rownames(zscores)),]<0)])
})
names(responders) <- rownames(zscores)



maf_filtered@clinical.data$Nutlin <- F
maf_filtered@clinical.data$Nutlin[which(maf_filtered@clinical.data$Tumor_Sample_Barcode %in% responders$`Nutlin-3`)] <- T
maf_filtered@clinical.data$Everolimus <- F
maf_filtered@clinical.data$Everolimus[which(maf_filtered@clinical.data$Tumor_Sample_Barcode %in% responders$Everolimus)] <- T
maf_filtered@clinical.data$Niraparib <- F
maf_filtered@clinical.data$Niraparib[which(maf_filtered@clinical.data$Tumor_Sample_Barcode %in% responders$Niraparib)] <- T
maf_filtered@clinical.data$Alpelisib <- F
maf_filtered@clinical.data$Alpelisib[which(maf_filtered@clinical.data$Tumor_Sample_Barcode %in% responders$Alpelisib)] <- T
maf_filtered@clinical.data$AZD4547 <- F
maf_filtered@clinical.data$AZD4547[which(maf_filtered@clinical.data$Tumor_Sample_Barcode %in% responders$AZD4547)] <- T
maf_filtered@clinical.data$EPZ015666 <- F
maf_filtered@clinical.data$EPZ015666[which(maf_filtered@clinical.data$Tumor_Sample_Barcode %in% responders$EPZ015666)] <- T
maf_filtered@clinical.data$Alpelisib <- F
maf_filtered@clinical.data$Tipifarnib[which(maf_filtered@clinical.data$Tumor_Sample_Barcode %in% responders$Tipifarnib)] <- T

nut_clinenrichment <- maftools::clinicalEnrichment(maf = maf_filtered, clinicalFeature = "Nutlin",pathways = F)
evo_clinenrichment <- maftools::clinicalEnrichment(maf = maf_filtered, clinicalFeature = "Everolimus",pathways = F)
nira_clinenrichment <- maftools::clinicalEnrichment(maf = maf_filtered, clinicalFeature = "Niraparib",pathways = F)
alpe_clinenrichment <- maftools::clinicalEnrichment(maf = maf_filtered, clinicalFeature = "Alpelisib",pathways = F)
azd_clinenrichment <- maftools::clinicalEnrichment(maf = maf_filtered, clinicalFeature = "AZD4547",pathways = F)
epz_clinenrichment <- maftools::clinicalEnrichment(maf = maf_filtered, clinicalFeature = "EPZ015666",pathways = F)
tipi_clinenrichment <- maftools::clinicalEnrichment(maf = maf_filtered, clinicalFeature = "Tipifarnib",pathways = F)
write.table(nut_clinenrichment$groupwise_comparision, "nutlin_clinical_enrichment.csv")
write.table(evo_clinenrichment$groupwise_comparision, "everolimus_clinical_enrichment.csv")
write.table(nira_clinenrichment$groupwise_comparision, "niraparib_clinical_enrichment.csv")
write.table(alpe_clinenrichment$groupwise_comparision, "alpelisib_clinical_enrichment.csv")
write.table(azd_clinenrichment$groupwise_comparision, "azd4547_clinical_enrichment.csv")
write.table(epz_clinenrichment$groupwise_comparision, "epz015666_clinical_enrichment.csv")
write.table(tipi_clinenrichment$groupwise_comparision, "tipifarnib_clinical_enrichment.csv")



genecopy_filtered <- genecopytable_all[,colnames(zscores)]


Nutlin_genecopy <- data.frame("average_responders" = rowMeans(genecopy_filtered[,responders$`Nutlin-3`]),
                              "average_non_responders" = rowMeans(genecopy_filtered[,-which(colnames(genecopy_filtered) %in% responders$`Nutlin-3`)]),
                              "t-test" = 1)
for(x in c(1:length(rownames(genecopy_filtered)))) {
  if(!any(is.na(genecopy_filtered[x,]))){
    Nutlin_genecopy[x,3] <- t.test(x = genecopy_filtered[x,which(colnames(genecopy_filtered) %in% responders$`Nutlin-3`)],
           y = genecopy_filtered[x,-which(colnames(genecopy_filtered) %in% responders$`Nutlin-3`)])$p.value}
}

Everolimus_genecopy <- data.frame("average_responders" = rowMeans(genecopy_filtered[,responders$Everolimus]),
                              "average_non_responders" = rowMeans(genecopy_filtered[,-which(colnames(genecopy_filtered) %in% responders$Everolimus)]),
                              "t-test" = 1)
for(x in c(1:length(rownames(genecopy_filtered)))) {
  if(!any(is.na(genecopy_filtered[x,]))){
    Everolimus_genecopy[x,3] <- t.test(x = genecopy_filtered[x,which(colnames(genecopy_filtered) %in% responders$Everolimus)],
                                   y = genecopy_filtered[x,-which(colnames(genecopy_filtered) %in% responders$Everolimus)])$p.value}
}

Niraparib_genecopy <- data.frame("average_responders" = rowMeans(genecopy_filtered[,responders$Niraparib]),
                                  "average_non_responders" = rowMeans(genecopy_filtered[,-which(colnames(genecopy_filtered) %in% responders$Niraparib)]),
                                  "t-test" = 1)
for(x in c(1:length(rownames(genecopy_filtered)))) {
  if(!any(is.na(genecopy_filtered[x,]))){
    Niraparib_genecopy[x,3] <- t.test(x = genecopy_filtered[x,which(colnames(genecopy_filtered) %in% responders$Niraparib)],
                                       y = genecopy_filtered[x,-which(colnames(genecopy_filtered) %in% responders$Niraparib)])$p.value}
}

Alpelisib_genecopy <- data.frame("average_responders" = rowMeans(genecopy_filtered[,responders$Alpelisib]),
                                 "average_non_responders" = rowMeans(genecopy_filtered[,-which(colnames(genecopy_filtered) %in% responders$Alpelisib)]),
                                 "t-test" = 1)
for(x in c(1:length(rownames(genecopy_filtered)))) {
  if(!any(is.na(genecopy_filtered[x,]))){
    Alpelisib_genecopy[x,3] <- t.test(x = genecopy_filtered[x,which(colnames(genecopy_filtered) %in% responders$Alpelisib)],
                                      y = genecopy_filtered[x,-which(colnames(genecopy_filtered) %in% responders$Alpelisib)])$p.value}
}

AZD4547_genecopy <- data.frame("average_responders" = rowMeans(genecopy_filtered[,responders$AZD4547]),
                                 "average_non_responders" = rowMeans(genecopy_filtered[,-which(colnames(genecopy_filtered) %in% responders$AZD4547)]),
                                 "t-test" = 1)
for(x in c(1:length(rownames(genecopy_filtered)))) {
  if(!any(is.na(genecopy_filtered[x,]))){
    AZD4547_genecopy[x,3] <- t.test(x = genecopy_filtered[x,which(colnames(genecopy_filtered) %in% responders$AZD4547)],
                                      y = genecopy_filtered[x,-which(colnames(genecopy_filtered) %in% responders$AZD4547)])$p.value}
}

EPZ015666_genecopy <- data.frame("average_responders" = rowMeans(genecopy_filtered[,responders$EPZ015666]),
                               "average_non_responders" = rowMeans(genecopy_filtered[,-which(colnames(genecopy_filtered) %in% responders$EPZ015666)]),
                               "t-test" = 1)
for(x in c(1:length(rownames(genecopy_filtered)))) {
  if(!any(is.na(genecopy_filtered[x,]))){
    EPZ015666_genecopy[x,3] <- t.test(x = genecopy_filtered[x,which(colnames(genecopy_filtered) %in% responders$EPZ015666)],
                                    y = genecopy_filtered[x,-which(colnames(genecopy_filtered) %in% responders$EPZ015666)])$p.value}
}

Tipifarnib_genecopy <- data.frame("average_responders" = rowMeans(genecopy_filtered[,responders$Tipifarnib]),
                                 "average_non_responders" = rowMeans(genecopy_filtered[,-which(colnames(genecopy_filtered) %in% responders$Tipifarnib)]),
                                 "t-test" = 1)
for(x in c(1:length(rownames(genecopy_filtered)))) {
  if(!any(is.na(genecopy_filtered[x,]))){
    Tipifarnib_genecopy[x,3] <- t.test(x = genecopy_filtered[x,which(colnames(genecopy_filtered) %in% responders$Tipifarnib)],
                                      y = genecopy_filtered[x,-which(colnames(genecopy_filtered) %in% responders$Tipifarnib)])$p.value}
}

pdf("drug_response_copynumbers_top50.pdf")
pheatmap(Nutlin_genecopy[order(Nutlin_genecopy$t.test, decreasing = F)[c(1:50)],c(1,2)], cluster_rows = T, cluster_cols = F, main = "Nutlin gene copies")
pheatmap(Everolimus_genecopy[order(Everolimus_genecopy$t.test, decreasing = F)[c(1:50)],c(1,2)], cluster_rows = T, cluster_cols = F, main = "Everolimus gene copies")
pheatmap(Niraparib_genecopy[order(Niraparib_genecopy$t.test, decreasing = F)[c(1:50)],c(1,2)], cluster_rows = T, cluster_cols = F, main = "Niraparib gene copies")
pheatmap(Alpelisib_genecopy[order(Alpelisib_genecopy$t.test, decreasing = F)[c(1:50)],c(1,2)], cluster_rows = T, cluster_cols = F, main = "Alpelisib gene copies")
pheatmap(AZD4547_genecopy[order(AZD4547_genecopy$t.test, decreasing = F)[c(1:50)],c(1,2)], cluster_rows = T, cluster_cols = F, main = "AZD4547 gene copies")
pheatmap(EPZ015666_genecopy[order(EPZ015666_genecopy$t.test, decreasing = F)[c(1:50)],c(1,2)], cluster_rows = T, cluster_cols = F, main = "EPZ015666 gene copies")
pheatmap(Tipifarnib_genecopy[order(Tipifarnib_genecopy$t.test, decreasing = F)[c(1:50)],c(1,2)], cluster_rows = T, cluster_cols = F, main = "Tipifarnib gene copies")
pheatmap(AZD4547_genecopy[grep("FGF", rownames(AZD4547_genecopy)),c(1,2)], cluster_rows = T, cluster_cols = F, main = "AZD4547 FGF gene copies (not significant)")
dev.off()

