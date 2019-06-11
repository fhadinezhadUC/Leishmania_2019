annotate.final.geneset <- function(genefilepath){
  # genefilepath is the path to the file integrated_tse_ara.txt
  library(gsubfn)
  library(ggplot2)
  library(seqinr)
  library(gdata)
  #source("https://bioconductor.org/biocLite.R")
  library(Biostrings)
  genefilepath <-
    "/home/fatemeh/Leishmania_2019/Leishmania_2019/Results/Integrated_Genes/integrated_tse_ara.txt"
  genefile <-
    read.table(genefilepath, header = TRUE, colClasses = "character")
  
  # cuttoff score for TSE 50, for ARA 107
  genefile[genefile$arascore == "notfound", ]$arascore <- -1
  genefile$arascore <- as.integer(genefile$arascore)
  genefile[genefile$tsescore == "notfound", ]$tsescore <- -1
  genefile$tsescore <- as.integer(genefile$tsescore)
  
  dismiss0 <-
    (genefile$foundby == "both") &
    (genefile$arascore < 107 | genefile$tsescore < 50)
  dismiss1 <-
    (genefile$foundby == "ara") &
    (genefile$arascore < 107) # most of these genes were pseudo or truncated or both
  dismiss2 <- (genefile$foundby == "tse") & (genefile$tsescore < 50)
  genefile <- genefile[(!dismiss0 & !dismiss1 & !dismiss2), ]
  
  # genefunc column is added which shows the presented gene identity in table
  # Genes marked as ?? include: pseudo|truncated genes, genes with unmatched identity, genes with unassigned identity|anticodon by any of genefinders.
  
  ambiguty1 <-
    genefile$tsefunc == "" &  genefile$arafunc == "" # 2 genes
  ambiguty2 <-
    genefile$tsenote != "notfound" # 2 genes not the same 2 genes in ambiguty1
  ambiguty3 <-
    genefile$foundby == "both" &
    genefile$tsefunc != genefile$arafunc # 22 genes
  ambiguty4 <- logical(length = nrow(genefile))
  for (i in 1:nrow(genefile)) {
    if (genefile$foundby[i] != "ara")
      ambiguty4[i] <- length(grep("n|N", genefile$tsegeneseq[i])) == 1
    else
      ambiguty4[i] <- length(grep("n|N", genefile$arageneseq[i])) == 1
  }
  
  ambiguties <- ambiguty1 | ambiguty2 | ambiguty3 | ambiguty4
  genefile$genefunc <- ""
  genefile[ambiguties, ]$genefunc <- "??"
  
  #table(genefile[!ambiguties, ]$foundby)
  #ara both 
  #36 3525 
  
  genefile[!ambiguties, ]$genefunc <-
    genefile[!ambiguties, ]$arafunc
  
  outputpath <-
    "/home/fatemeh/Leishmania_2019/Leishmania_2019/Document/"
  create.summary.table(genefile)
  prepare.tsfm.input(genefile)
  
}
prepare.tsfm.input <- function(genefile){
  library(gsubfn)
  library(ggplot2)
  library(seqinr)
  tsfmInput <- genefile[genefile$genefunc!="??",]
  # Write the table in a file and the rest will be done in TriTrypAlignment.R script
  # Write one file to have all the data to read everytime we need the genefile
  resultpath <- "/home/fatemeh/Leishmania_2019/Leishmania_2019/Results/Integrated_Genes/"
  write.table(tsfmInput,col.names = TRUE,file = paste(resultpath,"tsfm_input_geneset.txt",sep = ""))
  # read in the HomoC tRNA genes downloaded from http://gtrnadb.ucsc.edu/GtRNAdb2/genomes/eukaryota/Hsapi38/Hsapi38-seq.html
  #Homo_geneset <- read.fasta("/home/fatemeh/Leishmania_2019/Leishmania_2019/Results/tsfminput_final/hg38-tRNAs.fa",seqtype="DNA",as.string = TRUE)
  # create a dataframe with columns
  #merged both fasta files 
  #cat hg38-tRNAs.fa tsfm_input_geneset_EditedCovea.fasta > tsfm_finalinput_HomoC.fasta
  #covea TRNA2-euk.cm tsfm_finalinput_HomoC.fasta > tsfm_finalinput_HomoC.covea
  #WARNING: Warning: unrecognized character - in sequence
  # edit the file tsfm_finalinput_HomoC.covea using the functions in TriTrypAlignment.R and save the result in tsfm_input_geneset_EditedCovea.fasta
  # Split genes based on genome names or clusters (clusters of Leishmania and Trypanosoma) using the script GenomeClustering.R
  # Use script splitFuncClass.sh to split the genes based on their function generate the input tsfm files in clustal format
}
create.summary.table <- function(genefile) {
  # create a table with columns: Annotation, Intersection, ARAonly, Union
  
  firstcol <- c("#tRNA","#N/#G","Min Gene Length","Max Gene Length","%intron","%G" ,"%C" ,"%T" ,"%A","A",  "C",  "D",  "E",  "F",  "G",  "H",  "I",  "K",  "L",  "M",  "N",  "P",  "Q",  "R",  "S",  "T",  "V",  "W",  "X",  "Y",  "Z", "??")
  resulttable <- data.frame(firstcol,rep(0,length(firstcol)),rep(0,length(firstcol)),rep(0,length(firstcol)))
  names(resulttable) <- c("Annotation", "Intersection", "ARAonly", "Union")
  
  isboth <- genefile$foundby=="both"
  isaraonly <- genefile$foundby=="ara"
  intersectDF <- genefile[isboth,]
  ARAonlyDF <- genefile[isaraonly,]
  genefile$geneseq <- ""
  genefile[isaraonly,]$geneseq <- genefile[isaraonly,]$arageneseq
  genefile[!isaraonly,]$geneseq <- genefile[!isaraonly,]$tsegeneseq
  
  resulttable[1,2:4] <- c(nrow(intersectDF),nrow(ARAonlyDF),nrow(genefile))
  resulttable[2,2:4] <- c(sum(nchar(intersectDF$tsegeneseq))/nrow(intersectDF),sum(nchar(ARAonlyDF$arageneseq))/nrow(ARAonlyDF),sum(nchar(genefile$geneseq))/nrow(genefile)) # sum of gene length / # genes from the first row
  resulttable[3,2:4] <- c(min(nchar(intersectDF$tsegeneseq)),min(nchar(ARAonlyDF$arageneseq)),min(nchar(genefile$geneseq)))
  resulttable[4,2:4] <- c(max(nchar(intersectDF$tsegeneseq)),max(nchar(ARAonlyDF$arageneseq)),max(nchar(genefile$geneseq)))
  resulttable[5, 2:4] <-
    c(
      100 * nrow(bothdf[bothdf$tseintronbegin != 0, ]) / nrow(bothdf),
      100 * nrow(ARAonlyDF[ARAonlyDF$araintronbegin != "nointron", ]) / nrow(ARAonlyDF),
      100 * (nrow(bothdf[bothdf$tseintronbegin != 0, ]) + nrow(ARAonlyDF[ARAonlyDF$araintronbegin !=
                                                                           "nointron", ])) / nrow(genefile)
    )
  resulttable[6,2:4] <- c(Gpercentage(paste(intersectDF$tsegeneseq,collapse = "")),Gpercentage(paste(ARAonlyDF$arageneseq,collapse = "")),Gpercentage(paste(genefile$geneseq,collapse = "")))
  resulttable[7,2:4] <- c(Cpercentage(paste(intersectDF$tsegeneseq,collapse = "")),Cpercentage(paste(ARAonlyDF$arageneseq,collapse = "")),Cpercentage(paste(genefile$geneseq,collapse = "")))
  resulttable[8,2:4] <- c(Tpercentage(paste(intersectDF$tsegeneseq,collapse = "")),Tpercentage(paste(ARAonlyDF$arageneseq,collapse = "")),Tpercentage(paste(genefile$geneseq,collapse = "")))
  resulttable[9,2:4] <- c(Apercentage(paste(intersectDF$tsegeneseq,collapse = "")),Apercentage(paste(ARAonlyDF$arageneseq,collapse = "")),Apercentage(paste(genefile$geneseq,collapse = "")))
  
  # make a table of Class frequencies for each set:
  intersect_ClassFreq <- as.data.frame(table(intersectDF$genefunc))
  names(intersect_ClassFreq) <- c("class","freq")
  AraOnly_ClassFreq <- as.data.frame(table(ARAonlyDF$genefunc))
  names(AraOnly_ClassFreq) <- c("class","freq")
  Union_ClassFreq <- as.data.frame(table(genefile$genefunc))
  names(Union_ClassFreq) <- c("class","freq")
  intersect_ClassFreq$class <- as.character(intersect_ClassFreq$class)
  AraOnly_ClassFreq$class <- as.character(AraOnly_ClassFreq$class)
  Union_ClassFreq$class <- as.character(Union_ClassFreq$class)
  
  for (i in 1:nrow(intersect_ClassFreq)) {
    curr_class <- intersect_ClassFreq$class[i]
    resulttable[resulttable$Annotation==curr_class,][,2] <- intersect_ClassFreq$freq[i]
  }
  for (i in 1:nrow(AraOnly_ClassFreq)) {
    curr_class <- AraOnly_ClassFreq$class[i]
    resulttable[resulttable$Annotation==curr_class,][,3] <- AraOnly_ClassFreq$freq[i]
  }
  for (i in 1:nrow(Union_ClassFreq)) {
    curr_class <- Union_ClassFreq$class[i]
    resulttable[resulttable$Annotation==curr_class,][,4] <- Union_ClassFreq$freq[i]
  }
  
  resulttable$Intersection <- round(resulttable$Intersection,digits = 0)
  resulttable$ARAonly <- round(resulttable$ARAonly,digits = 0)
  resulttable$Union <- round(resulttable$Union,digits = 0)
  resulttable$Annotation <- as.character(resulttable$Annotation)
  Latex.file.prep(resulttable,"resulttable")
  
}
gcContent <- function(x) {
  x <- DNAString(toupper(x))
  alf <- alphabetFrequency(x, as.prob = TRUE)
  sum(alf[c("G", "C")]) * 100
}
atContent <- function(x) {
  x <- DNAString(toupper(x))
  alf <- alphabetFrequency(x, as.prob = TRUE)
  sum(alf[c("A", "T")]) * 100
}
Apercentage <- function(x) {
  x <- DNAString(toupper(x))
  alf <- alphabetFrequency(x, as.prob = TRUE)
  sum(alf[c("A")]) * 100
}
Cpercentage <- function(x) {
  x <- DNAString(toupper(x))
  alf <- alphabetFrequency(x, as.prob = TRUE)
  sum(alf[c("C")]) * 100
}
Gpercentage <- function(x) {
  x <- DNAString(toupper(x))
  alf <- alphabetFrequency(x, as.prob = TRUE)
  sum(alf[c("G")]) * 100
}
Tpercentage <- function(x) {
  x <- DNAString(toupper(x))
  alf <- alphabetFrequency(x, as.prob = TRUE)
  sum(alf[c("T")]) * 100
}

Latex.file.prep <- function(df,filename){
  firstline <- ""
  Llines <- paste(firstline, collapse = "&")
  Llines <- paste(Llines,"\\\\",sep="")
  #df$clusters <- as.character(df$clusters)
  for (i in 1:nrow(df)) {
    currline <- paste(df[i,],collapse = "&")
    currline <- paste(currline,"\\\\",sep="")
    Llines <- c(Llines,currline)
  }
  filename <- paste("/home/fatemeh/Leishmania_2019/Leishmania_2019/Results/",clade,".txt",sep = "")
  writeLines(Llines, filename,sep = "\n")
  
}
