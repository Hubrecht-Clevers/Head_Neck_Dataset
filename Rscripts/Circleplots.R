library(circlize)

# function to make a circleplot with genomic info, such as structual variants, SNPS and copynumber variations.
# this function can use a vcf_file as input, as well as copynumber info from the FREEC-package.
# Structual variants can be inputted straight from the MANTA, structual variant caller. 
# Genomes can be selected, so far only hg38 will give you an output.
# Please specify whether your sample was sequenced WGS or WES and for correct copy number calculations, sex info is needed.
make_circos_plot <- function(snv_in, ratiofile_in, sv_in = NULL, genome = "hg38", samplename, outputfile, WES = F, SEX="M"){
  # get all the SNVs from the file and make them into a table
  circos.clear()
  snv_file_info <- file.info(snv_in)
  # checking to see if there are snps in your vcf-file.
  test1 <- system(paste0("cat ", snv_in, " | wc -l"), intern = T)
  test2 <- system(paste0("cat ", snv_in, " | grep ^# | wc -l"), intern = T)
  if (test1 != test2) {
    # convert your vcf file to a workable table
    snv_tmp <- read.table(snv_in,comment.char="#")
    snv <- as.data.frame(cbind(paste0("chr",as.character(snv_tmp[,1])), snv_tmp[,2], snv_tmp[,2]+1))
    colnames(snv) <- c("CHROM", "START", "END")
    if(any(snv$CHROM %in% c("chrMT", grep("chrUn", snv$CHROM, value = T)))){
      snv <- snv[-which(snv$CHROM %in% c("chrMT", grep("chrUn", snv$CHROM, value = T))),]
    }
    snv$CHROM <- as.character(snv$CHROM)
    snv$START <- as.numeric(snv$START)
    snv$END <- as.numeric(snv$END)
  }
  if (test1 == test2) {
    message("No SNPs found in coding regions")
  }
  
  # get the copynumbers from the ratiofiles
  cnv <- read.table(ratiofile_in, header = T)
  if(WES){
    cnv2 <- data.frame("Chromosome" = cnv$Chromosome, "Start" = cnv$Start, "End" = unlist(lapply(cnv$Gene, function(x){strsplit(x, "-")[[1]][2]})), 
                       "Ratio" = cnv$Ratio, "CopyNumber" = cnv$CopyNumber)
    # correct for sex differences in chromosome copy numbers
    if (SEX == "F") {
      cnv2$CopyNumber[which(cnv2$Chromosome == "Y")] <- 2
    }
    if (SEX == "M") {
      cnv2$CopyNumber[which(cnv2$Chromosome == "X")] <- cnv2$CopyNumber[which(cnv2$Chromosome == "X")]*2
      cnv2$CopyNumber[which(cnv2$Chromosome == "Y")] <- cnv2$CopyNumber[which(cnv2$Chromosome == "Y")]*2
    }
  }
  # if the data is from WGS, it will use bins.
  if (!WES) {
    cnv2 <- data.frame("Chromosome" = cnv$Chromosome, "Start" = cnv$Start, "End" = unlist(lapply(cnv$Start, function(x){x+999})), 
                       "Ratio" = cnv$Ratio, "CopyNumber" = cnv$CopyNumber)
    test <- lapply(seq(1,dim(cnv2)[[1]], 50)[-length(seq(1,dim(cnv2)[[1]], 50))], function(x){
      if(all(cnv2$Chromosome[c(x,x+49)] == cnv$Chromosome[x])){
        chr <- cnv2$Chromosome[x]
        start <- cnv2$Start[x]
        end <- cnv2$End[x+49]
        ratio <- median(cnv2$Ratio[c(x:x+49)])
        # correct for sex differences in chromosome copy numbers
        if (chr == "X") {
          if (SEX == "M") {
            copynumber <- median(cnv2$CopyNumber[c(x:x+49)])*2
          }
        }
        if (chr == "Y"){
          if(SEX == "M"){
            copynumber <- median(cnv2$CopyNumber[c(x:x+49)])*2
          }
          if(SEX == "F"){
            copynumber <- 2
          }}
        if (chr %in% c(1:22)) {
          copynumber <- median(cnv2$CopyNumber[c(x:x+49)])
        }
        tmp <- data.frame("Chromosome" = chr, "Start" = start, "End" = end, 
                          "Ratio" = ratio, "CopyNumber" = copynumber)
        return(tmp)
      }
      
    })
    test2 <- test[-which(unlist(lapply(test, function(x){length(x)})) != 5)]
    test2 <- do.call("rbind",test2)
    cnv2 <- test2
  }
  # gather all duplications
  dup <- NULL
  if (any(cnv2$CopyNumber>2)) {
    
    list()
    dup=cnv2[cnv2$CopyNumber>2,]
    dup[,1]=paste("chr",as.character(dup[,1]),sep="")
    dup$Chromosome <- as.character(dup$Chromosome)
    dup$Start <- as.numeric(dup$Start)
    dup$End <- as.numeric(dup$End)
  }
  if (any(cnv2$CopyNumber < 2)) {
    
    del <- NULL
    # gather all deletions
    del=cnv2[cnv2$CopyNumber<2,]
    del[,1]=paste("chr",as.character(del[,1]),sep="")
    del$Chromosome <- as.character(del$Chromosome)
    del$Start <- as.numeric(del$Start)
    del$End <- as.numeric(del$End)
  }
  # make the structual variants table ready for use
  if (!is.null(sv_in)) {
    sv <- read.table(sv_in,comment.char="#")
    colnames(sv) <- c("CHROM", "START", "TYPE", "REF", "ALT", "V6", "FILTER", "INFO", "GT_legend", "GT")
    sv <- sv[sv$FILTER == "PASS",]
    sv <- sv[grep("1/1", sv$GT),]
    all_chrom <- unique(unlist(c(lapply(sv$CHROM, function(x){paste0("chr",x)}), "chrUnxx")))
    if(any(unlist(lapply(sv$TYPE, function(x){strsplit(x, ":")[[1]][1]})) == "MantaBND")){
      sv_bnd <- sv[grep("BND", sv$TYPE),]
      svTable=data.frame(paste0("chr", sv_bnd[,1]), as.numeric(sv_bnd[,2]), 
                         paste0("chr",as.character(unlist(lapply(sv_bnd$ALT, function(x){strsplit(strsplit(x,"]")[[1]][2],":")[[1]][1]})))),
                         as.numeric(unlist(lapply(sv_bnd$ALT, function(x){strsplit(strsplit(x,"]")[[1]][2],":")[[1]][2]}))))
      svTable <- svTable[-which(is.na(svTable[,4])),]
    }
    # Filter out all variants that map partly or completely to unknown chromosomes
    if(any(svTable[,1] %in% grep("chrUn", all_chrom, value = T))){
      if(!all(svTable[,1] %in% grep("chrUn", all_chrom, value = T))){
        if(length(grep("Un", svTable[,1])) > 0){
          svTable <- svTable[-grep("Un", svTable[,1]),]}
      }}
    if(any(svTable[,1] == "chrMT")){
      if(!all(svTable[,1] == "chrMT")){
        if(length(grep("MT", svTable[,1])) > 0){
          svTable <- svTable[-grep("MT", svTable[,1]),]}
      }}
    if (any(svTable[,1] %in% c("chrMT", grep("chrUn", all_chrom, value = T)))) {
      sv_in <- NULL
    }
    if(any(svTable[,3] %in% grep("chrUn", all_chrom, value = T))){
      if(!all(svTable[,3] %in% grep("chrUn", all_chrom, value = T))){
        if(length(grep("Un", svTable[,3])) > 0){
          svTable <- svTable[-grep("Un", svTable[,3]),]}
      }}
    if(any(svTable[,3] == "chrMT")){
      if(!all(svTable[,3] == "chrMT")){
        if(length(grep("MT", svTable[,3])) > 0){
          svTable <- svTable[-grep("MT", svTable[,3]),]}
      }}
    if (any(svTable[,3] %in% c("chrMT", grep("chrUn", all_chrom, value = T)))) {
      sv_in <- NULL
    }}
  # make the actual plots
  circos.clear()
  circos.par("start.degree" = 90, cell.padding = c(0, 0, 0, 0))
  
  circos.initializeWithIdeogram(species = "hg38")
  title(samplename)
  if (snv_file_info$size > 170000) {
    circos.genomicTrackPlotRegion(snv,stack = T, panel.fun = function(region, value, ...) {
      circos.genomicPoints(region, value, cex = 0.05, pch = 9,col='orange' , ...)
    })
  }
  if (snv_file_info$size < 170000) {
    snv <- data.frame("col1" = NA, "col2" = NA, "col3" = NA)
    circos.genomicTrackPlotRegion(snv,stack = T, panel.fun = function(region, value, ...) {
      circos.genomicPoints(region, value, cex = 0.05, pch = 9,col='orange' , ...)
    })
  }
  if (!is.null(dup)) {
    circos.genomicTrackPlotRegion(dup, stack = TRUE,panel.fun = function(region, value, ...) {
      circos.genomicRect(region, value, col = "red",bg.border = NA, cex=1 , border = "red", ...)
    })}
  if(!is.null(del)){
    circos.genomicTrackPlotRegion(del, stack = TRUE,panel.fun = function(region, value, ...) {
      circos.genomicRect(region, value, col = "blue",border = "blue", bg.border = NA, cex=1 , ...)
    })}
  
  if (is.null(sv_in)){message("no structual variants found")}
  if (!is.null(sv_in)) {
    if(!is.null(svTable)){
      svTable <- svTable[-which(svTable$paste0..chr...sv_bnd...1.. == svTable$paste0..chr...as.character.unlist.lapply.sv_bnd.ALT..function.x...),]
      for (i in 1:dim(svTable)[1]) {
        circos.link(as.character(svTable[i,1]), as.numeric(svTable[i,2]), as.character(svTable[i,3]), as.numeric(svTable[i,4]))
      }
      legend(0.6,0.85,legend="SV-TRANSLOCATION",col="black",lty=1,cex=0.75,lwd=1.2,bty='n')
    }}
  legend(0.7,1.02,legend=c("SNV", "CNV-DUPLICATION","CNV-DELETION"),col=c("orange","red","blue"),pch=c(16,15,15),cex=0.75,title="Tracks:",bty='n')
  
  dev.copy2pdf(file = outputfile)
  
  #clear up mess
  circos.clear()
  snv_tmp <- NULL
  snv <- NULL
  cnv <- NULL
  cnv2 <- NULL
  del <- NULL
  dup <- NULL
  sv <- NULL
  sv_bnd <- NULL
  svTable <- NULL
}

#usage:
make_circos_plot("VCF_file_input", "RATIO_file_input", "SV_MANTA_input", "GENOME_VERSION", "NAME", "OUTPUT_PDF_FILENAME", "WES-Y/N", "SEX")