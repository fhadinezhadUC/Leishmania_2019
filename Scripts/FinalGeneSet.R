library(gsubfn)
library(ggplot2)
library(seqinr)

mergeData <- function(){
  genefile <- prepare.tsfmGeneset()
  genecountdf <- as.data.frame(table(genefile$sourceOrg))
  #genomefile <- genome.nuc.composition
  filepath <- "/home/fatemeh/Leishmania_2019/Leishmania_2019/Results/GenomeNucComposition.txt"
  genomefile <- read.table(filepath, header = TRUE, colClasses = "character")
  genomefile  <- genomefile[order(seqcounts),]
  genomefile$genecounts <- 0
  for (i in 1:nrow(genomefile)) {
    currname <- genomefile$organism_shortname[i]
    if(currname=="LmajorSD75.1")
      currname <- "LmajorSD75"
    genomefile$genecounts[i] <- genecountdf[genecountdf$Var1==currname,]$Freq
  }
  # add the column genesMissing
  # latex
  firstline <- c("organism","Aperc","Tperc","Cperc","Gperc","ATp","GCp","seqcounts","genecounts")
  Llines <- paste(firstline, collapse = "&")
  Llines <- paste(Llines,"\\\\",sep="")
  for (i in 1:nrow(genomefile)) {
    currline <- paste(genomefile[,2:ncol(genomefile)][i,],collapse = "&")
    currline <- paste(currline,"\\\\",sep="")
    Llines <- c(Llines,currline)
  }
  writeLines(Llines, "/home/fatemeh/Leishmania_2019/Leishmania_2019/Results/GenomeNucComposition_Latex.txt",sep = "\n")
}
mismatchedGenes.annotation <- function(genefile){
  resultpath <-
    "/home/fatemeh/Leishmania_2019/Leishmania_2019/Results/Integrated_Genes/"
  genefilepath <-
    "/home/fatemeh/Leishmania_2019/Leishmania_2019/Results/Integrated_Genes/integrated_tse_ara.txt"
  genefile <-
    read.table(genefilepath, header = TRUE, colClasses = "character")
  genefile[genefile$arascore=="notfound",]$arascore <- -1 
  genefile$arascore <- as.integer(genefile$arascore)
  genefile[genefile$tsescore=="notfound",]$tsescore <- -1 
  genefile$tsescore <- as.integer(genefile$tsescore)
  dismiss <- genefile$tsefunc == "" &  genefile$arafunc == ""
  genefile <- genefile[!dismiss, ]
  #ara both  tse 
  #741 3595   30 
  bothdf <- genefile[genefile$foundby=="both",]
  UndetDF <- bothdf[bothdf$tsefunc!=bothdf$arafunc,]
  #nrow(UndetDF)
  #[1] 35
  
}
pseudogenes.annotation <- function(){
  # genefile before filtering but after adding the extra columns
  
  #genefile[genefile$foundby == "both", ]$genefunc <- genefile[genefile$foundby == "both", ]$tsefunc
  #genefile[genefile$tsenote!="notfound",]$genefunc <- toupper(genefile[genefile$tsenote!="notfound",]$genefunc)
  #isarraonly <- genefile$foundby =="ara"
  #genefile$genefunc <- toupper(genefile$genefunc)
  
  #genefile[isarraonly & genefile$tsenote=="notfound",]$genefunc <- tolower(genefile[isarraonly& genefile$tsenote=="notfound",]$genefunc)
  genefile[genefile$foundby == "both", ]$genefunc <- (genefile[genefile$foundby == "both", ]$arafunc)
  genefile <- genefile[genefile$foundby=="both",]
  genefile[genefile$tsefunc!= genefile$arafunc, ]$genefunc <- tolower(genefile[genefile$tsefunc!= genefile$arafunc, ]$genefunc)
  # keep the aragorn function in a column while makeing clusterDF

  
  clusterdf <- create.clusterDF(genefile)
  bigclusters <- clusterdf[nchar(clusterdf$clusterSeq2) > 0,]
  temdf <- as.data.frame(table(bigclusters$clusterSeq2))
  temdf$Var1 <- as.character(temdf$Var1)
  islowerarr <- logical(length = nrow(tempdf))
  for (i in 1:nrow(temdf)) {
    islowerarr[i] <- all(grepl("[[:upper:]]", strsplit(temdf$Var1[i], "")[[1]]))
  }
  df2 <- temdf[!islowerarr,]
  df2C <- temdf[islowerarr,]
  # df2$foundby <- ""
  # for (i in 1:nrow(df2)) {
  #   df2$foundby[i] <- paste(clusterdf[clusterdf$clusterSeq2==df2$Var1[i],]$foundby,collapse = ",") 
  # }
  df2$arafunc <- ""
  df2$arascore <- ""
  df2$tsescore <- ""
  for (i in 1:nrow(df2)) {
    df2$arafunc[i] <- paste(unique(clusterdf[clusterdf$clusterSeq2==df2$Var1[i],]$ara),collapse = ",")
    df2$arascore[i] <- paste(unique(clusterdf[clusterdf$clusterSeq2==df2$Var1[i],]$arascore),collapse = ",")
    df2$tsescore[i] <- paste(unique(clusterdf[clusterdf$clusterSeq2==df2$Var1[i],]$tsescore),collapse = ",")
    
  }
  
  df2similar <- df2[,-c(4,5)]
  df2similar$aracluster <- gsub(",","",df2similar$arafunc)
  df2 <- df2similar
  Latex.file.prep(df2,"mismatchIdentity")
  # frequency of each cluster within other clusters
  df2$similarARAcounts <- 0
  df2$similarTSEcounts <- 0
  for (i in 1:nrow(df2)) {
    currclus<-df2$Var1[i]
    currclusara<-df2$aracluster[i]
    for (j in 1:nrow(df2C)) {
      if(length(grep(tolower(currclus),tolower(df2C$Var1[j]),invert = FALSE))==1)
        df2$similarTSEcounts[i] = df2$similarTSEcounts[i] + df2C$Freq[j]
      if(length(grep(tolower(currclusara),tolower(df2C$Var1[j]),invert = FALSE))==1)
        df2$similarARAcounts[i] = df2$similarARAcounts[i] + df2C$Freq[j]
      
    }
  }
  #create clusters
  #clusterdf <- create.clusterDF(genefile)
  #clusterfreqT <- cluster.feq.table(clusterdf,"T","L",20,22)
  # 
  #Tclus <- clusterfreqT[nchar(clusterfreqT$cluster) > 2, ]
  #Lclus <- clusterfreqL[nchar(clusterfreqL$cluster) > 2, ]
  cluster_setDF <- Cluster.variations(genefile)
  # for each cluster_set keep those that have at least one lower case letter in the clusters 
  cluster_setDF2 <- cluster_setDF
  cluster_setnumers <- unique(cluster_setDF2$cluster_set)
  for (i in 1:length(cluster_setnumers)) {
    j <- cluster_setnumers[i]
    currclusters <- as.character(cluster_setDF2[cluster_setDF2$cluster_set==j,]$clusters)
    cluStr <- paste(as.character(currclusters),collapse = "")
    if(all(grepl("[[:upper:]]", strsplit(cluStr, "")[[1]])))
    {cluster_setDF2 <- cluster_setDF2[cluster_setDF2$cluster_set!=j,]}
  }
  cluster_setDF2$clusters <- as.character(cluster_setDF2$clusters)
  cluster_setDF2 <- cluster_setDF2[order(cluster_setDF2$cluster_set,cluster_setDF2$clusters),]
  Latex.file.prep(cluster_setDF2,"cluster_araonly_sets")
  
}
create.cluster.tables.TLO <- function(genefile){
  ################## final genes set for tsfm #######################
  genefile$class <- toupper(substr(genefile$sourceOrg,1,1))
  clusterdf <- create.clusterDF(genefile)
  clusterfreqT <- cluster.feq.table(clusterdf,"T","L",20,22)
  clusterfreqL <- cluster.feq.table(clusterdf,"L","T",22,20)
  
  Tclus <- clusterfreqT[nchar(clusterfreqT$cluster) > 2, ]
  Lclus <- clusterfreqL[nchar(clusterfreqL$cluster) > 2, ]
  
  clusterdf_gclustered <- cluster.genome(clusterdf)
  clusterdf_gclustered[clusterdf_gclustered$GenomeClust=="AfricanTrypanosome",]$class <- "T_AF"
  clusterdf_gclustered[clusterdf_gclustered$GenomeClust=="AmericanTrypanosome",]$class <- "T_Am"
  
  clusterfreqTAF <- cluster.feq.table(clusterdf_gclustered,"T_AF","T_Am",6,13)
  clusterfreqTAM <- cluster.feq.table(clusterdf_gclustered,"T_Am","T_AF",13,6)
  TAFclus <- clusterfreqTAF[nchar(clusterfreqTAF$cluster) > 2, ]
  TAMclus <- clusterfreqTAM[nchar(clusterfreqTAM$cluster) > 2, ]
  
  Latex.file.prep(Tclus[,c(1,2,4,6)],"T")
  Latex.file.prep(Lclus[,c(1,2,4,6)],"L")
  Latex.file.prep(TAFclus[,c(1,2,4,6)],"TAF")
  Latex.file.prep(TAMclus[,c(1,2,4,6)],"TAM")
  
  Tclus0 <- clusterfreqT[nchar(clusterfreqT$cluster) > 3 & clusterfreqT$Perc1 > 5, ]
  Lclus0 <- clusterfreqL[nchar(clusterfreqL$cluster) > 3 & clusterfreqL$Perc1 > 5, ]
  Latex.file.prep(Tclus0[,c(1,2,4,6)],"T0")
  Latex.file.prep(Lclus0[,c(1,2,4,6)],"L0")
  ###################################################################
  
  #A   C   D   E   F   G   H   I   K   L   M   N   P   Q   R   S   T   V   W   X   Y   Z 
  #111  24  68  74  41 120  45  84  96 180  47  66 112  83 182 116 129 119  24  44  56  22 
  
  #A   C   D   E   F   G   H   I   K   L   M   N   P   Q   R   S   T   V   W   X   Y   Z 
  #80  33  29  70  58  89  31  73  70 125  42  48  69  67 134  96  71  94  22  22  23  47 
  
  #clusterdf <- create.clusterDF(genefile)
  # Tclusterfreq <- cluster.feq.table(clusterdf,"T")
  # Tclusterfreq <- Tclusterfreq[order(-Tclusterfreq$TgenomePerc),]
  # Lclusterfreq <- cluster.feq.table(clusterdf,"L")
  # Lclusterfreq <- Lclusterfreq[order(-Lclusterfreq$LgenomePerc),]
  # 
  # Tdf <- Tclusterfreq[,c(1,2,4,6)]
  # Ldf <- Lclusterfreq[,c(1,2,4,6)]
  # # subset clusters of length > 2 that have occurred more than once 
  # Tdf2 <- Tdf[nchar(Tdf$cluster)>3 & Tdf$TgenomePerc > 5, ]
  # Ldf2 <- Ldf[nchar(Ldf$cluster)>3 & Ldf$LgenomePerc > 5, ]
  # Latex.file.prep(Tdf2,"T")
  # Latex.file.prep(Ldf2,"L")
  # latex files
  
  
  # clusterdf$clusterSeq2 <- paste(clusterdf$clusterSeq2,clusterdf$clusterDir,sep = ",")
  # isbig <- nchar(clusterdf$clusterSeq2) > 3
  # bigclusterDF <- clusterdf[isbig,]
  # 
  # clustercounts <- as.data.frame(table(LbigclusterDF$clusterSeq2))
  # clustercounts$cluster <- ""
  # clustercounts$direction <- ""
  # clustercounts$Var1 <- as.character(clustercounts$Var1)
  # for (i in 1:nrow(clustercounts)) {
  #   clustercounts$cluster[i] <- unlist(strsplit(clustercounts$Var1[i],split = ","))[1]
  #   clustercounts$direction[i] <- unlist(strsplit(clustercounts$Var1[i],split = ","))[2]
  # }
  # 
  # clustercounts <- clustercounts[clustercounts$Freq > 9,]
  # topclusters <- clustercounts[order(clustercounts$cluster),]
  # tempfile <- genefile
  # tem <- read.table("/home/fatemeh/Leishmania_2019/Leishmania_2019/Results/Integrated_Genes/undetDF2.txt",header = FALSE,skip = 1,fill = TRUE)
  # IDto <- tem$V2
  # IDfrom <- tem$V4
  # shifts <- genefile$geneid %in% tem$V1
  # valid <- tem$V1 %in% genefile$geneid
  # genefile[shifts,]$genefunc <- paste("(",tolower(IDto[valid]),"<",tolower(IDfrom[valid]),")",sep = "")
}
Latex.file.prep <- function(df,clade){
  firstline <- c("cluster","dir","Perc1","Perc2")
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
cluster.feq.table <- function(clusterdf, clade1,clade2,genomecounts1,genomecounts2){
  clusterdf$clusterSeq3 <-
    paste(clusterdf$clusterSeq2, clusterdf$clusterDir, sep = ",")
  clusterdf1 <-
    clusterdf[clusterdf$class == clade1 & nchar(clusterdf$clusterSeq2) > 1, ]
  clusterdf2 <-
    clusterdf[clusterdf$class == clade2 &
                nchar(clusterdf$clusterSeq2) > 1, ]
  clusterFreq <- as.data.frame(table(clusterdf1$clusterSeq3))
  
  names(clusterFreq) <- c("clusterdir", "freq")
  clusterFreq$cluster <- ""
  clusterFreq$dir <- ""
  clusterFreq$freq1 <- 0
  clusterFreq$Perc1 <- 0
  clusterFreq$freq2 <- 0
  clusterFreq$Perc2 <- 0
  clusterFreq$clusterdir <- as.character(clusterFreq$clusterdir)
  for (i in 1:nrow(clusterFreq)) {
    clusterFreq$cluster[i] <-
      unlist(strsplit(clusterFreq$clusterdir[i], split = ","))[1]
    clusterFreq$dir[i] <-
      unlist(strsplit(clusterFreq$clusterdir[i], split = ","))[2]
    clusterFreq$freq1[i] <-
      nrow(clusterdf[clusterdf$clusterSeq2 == clusterFreq$cluster[i] &
                       clusterdf$clusterDir == clusterFreq$dir[i] &
                       clusterdf$class == clade1,])
    clusterFreq$freq2[i] <-
      nrow(clusterdf[clusterdf$clusterSeq2 == clusterFreq$cluster[i] &
                       clusterdf$clusterDir == clusterFreq$dir[i] &
                       clusterdf$class == clade2,])
    clusterFreq$Perc1[i] <-
      length(unique(clusterdf[clusterdf$clusterSeq2 == clusterFreq$cluster[i] & clusterdf$clusterDir == clusterFreq$dir[i] &
                                clusterdf$class == clade1, ]$sourceOrg))/genomecounts1
    clusterFreq$Perc2[i] <-
      length(unique(clusterdf[clusterdf$clusterSeq2 == clusterFreq$cluster[i] & clusterdf$clusterDir == clusterFreq$dir[i] &
                                clusterdf$class == clade2, ]$sourceOrg))/genomecounts2
  }
  clusterFreq$Perc1 <- round(clusterFreq$Perc1*100,digits = 0)
  clusterFreq$Perc2 <- round(clusterFreq$Perc2*100,digits = 0)
  outputdf <- clusterFreq[,3:8][order(-clusterFreq[,3:8]$Perc1),]
  outputdf
}
genome.nuc.composition <- function(){
  library(gsubfn)
  library(ggplot2)
  library(seqinr)
  filenames <-
    list.files(
      "/home/fatemeh/Leishmania_2019/GenomeData/",
      pattern = "*.fasta",
      full.names = TRUE
    )
  seqcounts <- integer(length = length(filenames))
  genomelength <- integer(length = length(filenames))
  Cperc <- integer(length = length(filenames))
  Gperc <- integer(length = length(filenames))
  Aperc <- integer(length = length(filenames))
  Tperc <- integer(length = length(filenames))
  GCp <- integer(length = length(filenames))
  ATp <- integer(length = length(filenames))
  geneCounts <- integer(length = length(filenames))
  aveGeneLen <- integer(length = length(filenames))
  organism_longname <- character(length = length(filenames))
  organism_shortname <- character(length = length(filenames))
  for (i in 1:length(filenames)) {
    seq = read.fasta(filenames[i],seqtype="DNA",as.string = TRUE)
    allseq <- paste(paste(seq,collapse = ""),collapse = "")
    allseq <- DNAString(allseq)
    GCp[i] <- gcContent(allseq)
    ATp[i] <- atContent(allseq)
    Tperc[i] <- Tpercentage(allseq)
    Aperc[i] <- Apercentage(allseq)
    Cperc[i] <- Cpercentage(allseq)
    Gperc[i] <- Gpercentage(allseq)
    genomelength[i] <- nchar(allseq)
    seqcounts[i] <- length(seq)
    temparr <- unlist(strsplit(as.character(getAnnot(seq[[1]])), " "))[3]
    organism_longname[i] <- substring(temparr,10)
    arr <- unlist(strsplit(filenames[i], "/"))
    shortname <- arr[length(arr)]
    shortname2 <- gsub("TriTrypDB-41_", "", shortname)
    shortname3 <- gsub("_Genome.fasta$", "", shortname2)
    organism_shortname[i] <- shortname3
  }
  Aperc <- round(Aperc,digits = 2)*100
  Tperc <- round(Tperc,digits = 2)*100
  Cperc <- round(Cperc,digits = 2)*100
  Gperc <- round(Gperc,digits = 2)*100
  ATp <- round(ATp,digits = 2)*100
  GCp <- round(GCp,digits = 2)*100
  genomeDF <- data.frame(organism_longname,organism_shortname,Aperc,Tperc,Cperc,Gperc,ATp,GCp,seqcounts)
  resultpath <- "/home/fatemeh/Leishmania_2019/Leishmania_2019/Results/"
  write.table(genomeDF,col.names = TRUE,file = paste(resultpath,"GenomeNucComposition.txt",sep = ""))
  
  filepath <- "/home/fatemeh/Leishmania_2019/Leishmania_2019/Results/GenomeNucComposition.txt"
  genomefile <- read.table(filepath, header = TRUE, colClasses = "character")
  # write with fixed width
  n <- data.frame("organism_longname","organism_shortname","Aperc","Tperc","Cperc","Gperc","ATp","GCp","seqcounts")
  names(n) <- c("organism_longname","organism_shortname","Aperc","Tperc","Cperc","Gperc","ATp","GCp","seqcounts")
  genomefile  <- genomefile[order(seqcounts),]
  genomefile2 <- genomefile
  genomefile[, ] <- lapply(genomefile[, ], as.character)
  genomefile$Aperc <- paste(genomefile$Aperc,"%",sep = "")
  genomefile$Cperc <- paste(genomefile$Cperc,"%",sep = "")
  genomefile$Gperc <- paste(genomefile$Gperc,"%",sep = "")
  genomefile$Tperc <- paste(genomefile$Tperc,"%",sep = "")
  genomefile$ATp <- paste(genomefile$ATp,"%",sep = "")
  genomefile$GCp <- paste(genomefile$GCp,"%",sep = "")
  library(gdata)
  write.fwf(
    rbind(n, genomefile),
    colnames = FALSE,
    width = c(57, 32, 10, 10, 10, 10, 10, 10, 10),
    file = paste(resultpath,"GenomeNucComposition2.txt",sep = "")
  )
  genomefile2
  # read like this:
 
}
cluster.genome <- function(genefile){
  genefile$GenomeClust <- ""
  for (i in 1:nrow(genefile)) {
    currGenome <- genefile$sourceOrg[i]
    if (currGenome == "LspMARLEM2494" |
        currGenome == "LenriettiiLEM3045")
      genefile$GenomeClust[i] <- "LenriettiComplex"
    if (currGenome == "TbruceigambienseDAL972" |
        currGenome == "TbruceiLister427" |
        currGenome == "TbruceiTREU927" |
        currGenome == "TevansiSTIB805" |
        currGenome == "TcongolenseIL3000" |
        currGenome == "TvivaxY486")
      genefile$GenomeClust[i] <- "AfricanTrypanosome"
    if (currGenome == "TgrayiANR4" |
        currGenome == "TrangeliSC58" |
        currGenome == "TcruziCLBrener" |
        currGenome == "TcruziCLBrenerEsmeraldo-like" |
        currGenome == "TcruziCLBrenerNon-Esmeraldo-like" |
        currGenome == "TcruzicruziDm28c" |
        currGenome == "TcruziDm28c" |
        currGenome == "TcruziEsmeraldo" |
        currGenome == "TcruziJRcl4" |
        currGenome == "TcruzimarinkelleiB7" |
        currGenome == "TcruziSylvioX10-1" |
        currGenome == "TcruziSylvioX10-1-2012" |
        currGenome == "TcruziTulacl2" )#|
        #currGenome == "TtheileriEdinburgh")
    genefile$GenomeClust[i] <- "AmericanTrypanosome"
    if (currGenome == "CfasciculataCfCl" |
        currGenome == "LseymouriATCC30220" |
        currGenome == "LpyrrhocorisH10")
      genefile$GenomeClust[i] <- "Leishmania1"
    if (currGenome == "LmajorFriedlin" |
        currGenome == "LmajorLV39c5" |
        currGenome == "LmajorSD75" |
        currGenome == "LturanicaLEM423" |
        currGenome == "LarabicaLEM1108" |
        currGenome == "LtropicaL590" |
        currGenome == "LaethiopicaL147" |
        currGenome == "LgerbilliLEM452")
      genefile$GenomeClust[i] <- "Leishmania2"
    if (currGenome == "LdonovaniBHU1220" |
        currGenome == "LdonovaniBPK282A1" |
        currGenome == "LinfantumJPCM5")
      genefile$GenomeClust[i] <- "LDonovaniComplex"
    if (currGenome == "LamazonensisMHOMBR71973M2269" |
        currGenome == "LmexicanaMHOMGT2001U1103")
      genefile$GenomeClust[i] <- "LMexicanaComplex"
    if (currGenome == "LbraziliensisMHOMBR75M2904" |
        currGenome == "LbraziliensisMHOMBR75M2903" |
        currGenome == "LpanamensisMHOMPA94PSC1" |
        currGenome == "LpanamensisMHOMCOL81L13")
      genefile$GenomeClust[i] <- "Lvianna"
  }
  #table(paste(genefile$GenomeClust,genefile$class,sep = "|"))
  #AfricanTrypanosome|T AmericanTrypanosome|T                    |B                    |E                    |L    LDonovaniComplex|L         Leishmania1|C         Leishmania1|L 
  #413                   980                    69                   103                    78                   253                   105                   198 
  #Leishmania2|L    LenriettiComplex|L    LMexicanaComplex|L             Lvianna|L                    |P 
  #672                   162                   149                   331                    61 
  genefile
}
# clustersize distribution 
clustersize.dist.visualize <- function(genefile){
  clusterdf <- create.clusterDF(genefile)
  
  # dataframe c("clusterLen","Tfreq","Lfreq","Ofreq","percOfgenes")
  maxcluslen <- max(nchar(clusterdf$clusterSeq2))
  sizeDF <- data.frame(seq(1,maxcluslen,1),rep(0,maxcluslen))
  names(sizeDF) <- c("ClusterLen","Tfreq")
  sizeDF$Lfreq <- 0
  sizeDF$Ofreq <- 0
  sizeDF$percOfgenes <- 0
  for (i in 1: nrow(sizeDF)) {
    sizeDF$Tfreq[i] <- nrow(clusterdf[clusterdf$class=="T" & nchar(clusterdf$clusterSeq2) == i,])
    sizeDF$Lfreq[i] <- nrow(clusterdf[clusterdf$class=="L" & nchar(clusterdf$clusterSeq2) == i,])
    sizeDF$Ofreq[i] <- nrow(clusterdf[clusterdf$class=="O" & nchar(clusterdf$clusterSeq2) == i,])
    size_i_freq <- sizeDF$Tfreq[i] + sizeDF$Lfreq[i] + sizeDF$Ofreq[i]
    gene_num <- size_i_freq * i 
    sizeDF$percOfgenes[i] <- round((gene_num/nrow(genefile))*100,digits = 0)
  }
  
  library(ggplot2)
  library(reshape2)
  #library(reshape)
  library(plyr)
  sizeDF2 <- melt(sizeDF, id = c("ClusterLen","percOfgenes"))
  # Sort the data by ClusterLen and variable
  # Calculate the cumulative sum of the variable value for each cluster
  # Create the plot
  #sizeDF2_sorted <- arrange(sizeDF2, ClusterLen, variable,value)
  sizeDF2$ClusterLen <- factor(sizeDF2$ClusterLen)
  sizeDF2$OLTorder <- 0
  sizeDF2[sizeDF2$variable == "Ofreq", ]$OLTorder <- 0
  sizeDF2[sizeDF2$variable == "Lfreq", ]$OLTorder <- 1
  sizeDF2[sizeDF2$variable == "Tfreq", ]$OLTorder <- 2
  sizeDF2_sorted <-
    sizeDF2[order(sizeDF2$ClusterLen, sizeDF2$OLTorder), ]
  #sizeDF2$variable <- as.factor(sizeDF2$variable)
  df_cumsum <- ddply(sizeDF2_sorted, "ClusterLen",
                     transform, label_ypos = cumsum(value))
  df_cumsum$label_ypos2 <- ""
  df_cumsum$label2 <- ""
  for (i in 1:nrow(df_cumsum)) {
    if (i %% 3 == 0)
    {
      df_cumsum$label_ypos2[i] <- toString(df_cumsum$label_ypos[i])
      df_cumsum$label2[i] <-
        paste(toString(df_cumsum$percOfgenes[i]), "%", sep = "")
    }
    else
    {
      df_cumsum$label_ypos2[i] = ""
      df_cumsum$label2[i] = ""
    }
  }
  df_cumsum$label1 <- ""
  df_cumsum$label1 <- df_cumsum$value
  df_cumsum[df_cumsum$value < 15,]$label1 <- ""
  ggplot(data = df_cumsum, aes(x = ClusterLen, y = value, fill = variable)) +
    geom_bar(stat = "identity") +
    geom_text(
      aes(y = label_ypos, label = label1),
      vjust = 1.3,
      color = "gray20",
      size = 4
    ) +
    scale_fill_brewer(palette = "Paired") + geom_text(aes(y = as.integer(label_ypos2), label =
                                                            label2),
                                                      vjust = -0.3,
                                                      size = 5,color = "darkgreen") + 
    scale_fill_manual(values = c("#00AFBB", "#E7B800", "lightgray"),name = "Genomes",labels = c("Trypanosoma", "Leishmania", "Others")) +
    theme_minimal() + xlab("Cluster Length") + ylab("Frequency")
  
}
Cluster.variations <- function(genefile){
  clusterdf <- create.clusterDF(genefile)
  bigclusters <- clusterdf[nchar(clusterdf$clusterSeq2) > 2,]
  clusters <- names(table(bigclusters$clusterSeq2))
  # percentage of pseudos that have occured as singleton
  
  cluster_set <-integer(length = length(clusters))
  setnumber <- 1
  for (i in 1:length(clusters)) {
    for (j in 1:length(clusters)) {
      S1 <-
        substring(tolower(clusters[i]), seq(1, nchar(clusters[i]), 1), seq(1, nchar(clusters[i]), 1))
      S2 <-
        substring(tolower(clusters[j]), seq(1, nchar(clusters[j]), 1), seq(1, nchar(clusters[j]), 1))
      C1 <- setdiff(S1, S2)
      C2 <- setdiff(S2, S1)
      C3 <- intersect(S1,S2)
      minL<- min(length(S1),length(S2))
      issubstr <- grep(tolower(clusters[i]),tolower(clusters[j]),invert = TRUE)
      if ((length(C1) < 1 | length(C2) < 1) ) {
        if (cluster_set[i] != 0 & cluster_set[j] == 0)
        {
          cluster_set[j] <-  cluster_set[i]
          #setnumber = setnumber + 1
        }
        
        if (cluster_set[j] != 0 & cluster_set[i] == 0)
        {
          cluster_set[i] <-  cluster_set[j]
          #setnumber = setnumber + 1
        }
        if (cluster_set[j] != 0 & cluster_set[i] != 0)
        {
          if(length(unique(S1)) > length(unique(S2)))
            cluster_set[j] <-  cluster_set[i]
          else
            cluster_set[i] <-  cluster_set[j]
          #setnumber = setnumber + 1
        }
        if (cluster_set[j] == 0 & cluster_set[i] == 0)
        {
          cluster_set[i] = setnumber
          cluster_set[j] = setnumber
          setnumber = setnumber + 1
        }
      }
    }
  }
  
  cluster_setDF <- data.frame(clusters,cluster_set)
  cluster_setDF <- cluster_setDF[order(cluster_setDF$cluster_set),]
  # for each cluster add the genomeclasses that have it 
  # and the direction of each classes
  cluster_setDF$Genomeclass <- ""
  cluster_setDF$Tdirs <- ""
  cluster_setDF$Ldirs <- ""
  cluster_setDF$Odirs <- ""
  
  for (i in 1:nrow(cluster_setDF)) {
    currclus <- cluster_setDF$clusters[i]
    cluster_setDF$Genomeclass[i] <- paste(unique(bigclusters[bigclusters$clusterSeq2==currclus,]$class),collapse = "" )
    cluster_setDF$Tdirs[i] <- paste(unique(bigclusters[bigclusters$clusterSeq2==currclus & bigclusters$class=="T",]$clusterDir),collapse = "|" )
    cluster_setDF$Ldirs[i] <- paste(unique(bigclusters[bigclusters$clusterSeq2==currclus & bigclusters$class=="L",]$clusterDir),collapse = "|" )
    cluster_setDF$Odirs[i] <- paste(unique(bigclusters[bigclusters$clusterSeq2==currclus & bigclusters$class=="O",]$clusterDir),collapse = "|" )
  }
  cluster_setDF <- cluster_setDF[order(cluster_setDF$cluster_set,cluster_setDF$Genomeclass),]
  Latex.file.prep(cluster_setDF,"cluster_araonly")
  # make a nother distance method based on diff == 0 
  # merge two arrays with the second array having the higher priority
  cluster_setDF
}
tree.visualization <- function(){
  # source("https://bioconductor.org/biocLite.R")
  # biocLite("BiocUpgrade") # you may need this
  biocLite("ggtree")
  library("ggtree")
  tree <- read.tree("/home/fatemeh/Leishmania_2019/Leishmania_2019/Document/tree2.nwk")
  #ggplot(tree) + geom_tree() + theme_tree()
  p <- ggtree(tree)
  #p + geom_nodepoint()
  #p + geom_tippoint()
  p + geom_tiplab(offset = 1,size=3,hjust =1)
  p+ geom_cladelabel(node=1, label="Some random clade", color="red")
  }
#source("https://bioconductor.org/biocLite.R")
library(Biostrings)

gcContent <-
  function(x)
  {
    x <- DNAString(toupper(x))
    alf <- alphabetFrequency(x, as.prob = TRUE)
    sum(alf[c("G", "C")])*100
  }
atContent <-
  function(x)
  {
    x <- DNAString(toupper(x))
    alf <- alphabetFrequency(x, as.prob = TRUE)
    sum(alf[c("A", "T")])*100
  }
Apercentage <-
  function(x)
  {
    x <- DNAString(toupper(x))
    alf <- alphabetFrequency(x, as.prob = TRUE)
    sum(alf[c("A")])*100
  }
Cpercentage <-
  function(x)
  {
    x <- DNAString(toupper(x))
    alf <- alphabetFrequency(x, as.prob = TRUE)
    sum(alf[c("C")])*100
  }
Gpercentage <-
  function(x)
  {
    x <- DNAString(toupper(x))
    alf <- alphabetFrequency(x, as.prob = TRUE)
    sum(alf[c("G")])*100
  }
Tpercentage <-
  function(x)
  {
    x <- DNAString(toupper(x))
    alf <- alphabetFrequency(x, as.prob = TRUE)
    sum(alf[c("T")])*100
  }
create.clusterDF <- function(tsfmGeneDF){
  geneset <- assignCluster(tsfmGeneDF,1000)
  m <- matrix(nrow = max(geneset$cluster),ncol = 13)
  clusterDF<- as.data.frame(m)
  names(clusterDF) <- c("sourceOrg","soureSeq","clusterID","clusterSeq1","clusterSeq2", "clusterBegin","clusterEnd","clusterDir","clusterLen","foundby","ara","arascore","tsescore")
  for (i in 1:max(geneset$cluster)) {
    clusterDF$sourceOrg[i] <- geneset[geneset$cluster==i,]$sourceOrg[1]
    clusterDF$soureSeq[i] <- geneset[geneset$cluster==i,]$sourceseq[1]
    clusterDF$clusterID[i] <- i
    clusterDF$clusterSeq1[i] <- paste(geneset[geneset$cluster==i,]$genecode,collapse = ",")
    clusterDF$clusterSeq2[i] <- paste(geneset[geneset$cluster==i,]$genefunc,collapse = "")
    clusterDF$clusterBegin[i] <- sort(geneset[geneset$cluster==i,]$begin)[1]
    clusterDF$clusterEnd[i] <- sort(geneset[geneset$cluster==i,]$end)[nrow(geneset[geneset$cluster==i,])]
    clusterDF$clusterDir[i] <- paste(geneset[geneset$cluster==i,]$direction,collapse = "")
    clusterDF$clusterLen[i] <- nrow(geneset[geneset$cluster==i,])
    clusterDF$foundby[i] <- paste(geneset[geneset$cluster==i,]$foundby,collapse = "|")
    clusterDF$ara[i] <- paste(geneset[geneset$cluster==i,]$tsefunc,collapse = ",")
    clusterDF$arascore[i] <- paste(geneset[geneset$cluster==i,]$arascore,collapse = ",")
    clusterDF$tsescore[i] <- paste(geneset[geneset$cluster==i,]$tsescore,collapse = ",")
  }
  clusterDF$class <- toupper(substr(clusterDF$sourceOrg,1,1))
  clusterDF[clusterDF$class!="T" & clusterDF$class!="L",]$class <- "O"
  clusterDF
}
assignCluster <- function(geneDF, clusterdistance) {
  geneDF <-
    geneDF[order(geneDF$sourceOrg, geneDF$sourceseq, geneDF$begin), ]
  geneDF$cluster <- 0
  geneDF$cluster[1] <- 1
  setnumber <- 1
  
  for (i in 1:(nrow(geneDF) - 1)) {
    geneDist <-
      abs(geneDF$begin[i + 1] - geneDF$begin[i])
    if ((geneDist < clusterdistance) &&
        (geneDF$sourceOrg[i] == geneDF$sourceOrg[i + 1]) &&
        (geneDF$sourceseq[i + 1] == geneDF$sourceseq[i]))
      geneDF$cluster[i + 1] <- setnumber
    else
    {
      setnumber <- setnumber + 1
      geneDF$cluster[i + 1] <- setnumber
    }
  }
  geneDF
}

prepare.tsfmGeneset <- function(){
  resultpath <-
    "/home/fatemeh/Leishmania_2019/Leishmania_2019/Results/Integrated_Genes/"
  genefilepath <-
    "/home/fatemeh/Leishmania_2019/Leishmania_2019/Results/Integrated_Genes/integrated_tse_ara.txt"
  genefile <-
    read.table(genefilepath, header = TRUE, colClasses = "character")
  
  # cuttoff score for TSE 50, for ARA 107
  genefile[genefile$arascore=="notfound",]$arascore <- -1 
  genefile$arascore <- as.integer(genefile$arascore)
  genefile[genefile$tsescore=="notfound",]$tsescore <- -1 
  genefile$tsescore <- as.integer(genefile$tsescore)
  
  dismiss0 <- (genefile$foundby=="both") & (genefile$arascore < 107 | genefile$tsescore < 50)
  dismiss1 <- (genefile$foundby=="ara") & (genefile$arascore < 107) # most of these genes were pseudo or truncated or both
  dismiss2 <- (genefile$foundby=="tse") & (genefile$tsescore < 50)
  genefile <- genefile[(!dismiss0 & !dismiss1 & !dismiss2),]
  #ara both  tse 
  #36 3554    1 
  # report these genes
  #dismiss4 <- genefile1$tsefunc == "" &  genefile1$arafunc == ""
  #genefile2 <- genefile1[!dismiss4, ]
  #nonpseudo_nontrunc <- genefile2$tsenote == "notfound"
  #genefile2 <- genefile2[nonpseudo_nontrunc, ]
  # Add a column begin, end and genecode to the dataframe to be able to assign clusters later
  # genefile$begin <- ""
  # genefile$end <- ""
  # genefile$genecode <- ""
  # genefile[genefile$foundby == "both" |
  #            genefile$foundby == "tse", ]$begin <-
  #   genefile[genefile$foundby == "both" |
  #              genefile$foundby == "tse", ]$tsebegin
  # genefile[genefile$foundby == "ara", ]$begin <-
  #   genefile[genefile$foundby == "ara", ]$arabegin
  # genefile[genefile$foundby == "both" |
  #            genefile$foundby == "tse", ]$end <-
  #   genefile[genefile$foundby == "both" |
  #              genefile$foundby == "tse", ]$tseend
  # genefile[genefile$foundby == "ara", ]$end <-
  #   genefile[genefile$foundby == "ara", ]$araend
  # 
  # genefile[genefile$foundby == "both" |
  #            genefile$foundby == "tse", ]$genecode <-
  #   paste("(",tolower(genefile[genefile$foundby == "both" |
  #                                genefile$foundby == "tse", ]$tseac),")",genefile[genefile$foundby == "both" |
  #                                                                                   genefile$foundby == "tse", ]$tsefunc,sep="")
  # genefile[genefile$foundby == "ara", ]$genecode <-
  #   paste(genefile[genefile$foundby == "ara", ]$araac,genefile[genefile$foundby == "ara", ]$arafunc,sep="")
  # 
  
  # genefunc column is added which shows the presented gene identity in table
  # Genes marked as ?? include: pseudo|truncated genes, genes with unmatched identity, genes with unassigned identity|anticodon by any of genefinders.
  
  ambiguty1 <- genefile$tsefunc == "" &  genefile$arafunc == "" # 2 genes
  ambiguty2 <- genefile$tsenote != "notfound" # 2 genes not the same 2 genes in ambiguty1
  ambiguty3 <- genefile$foundby=="both" &  genefile$tsefunc != genefile$arafunc # 22 genes
  ambiguty4 <- logical(length = nrow(genefile))
  for (i in 1:nrow(genefile)) {
    if(genefile$foundby[i]!="ara")
        ambiguty4[i] <- length(grep("n|N",genefile$tsegeneseq[i]))==1
    else
        ambiguty4[i] <- length(grep("n|N",genefile$arageneseq[i]))==1
  }
  
  ambiguties <- ambiguty1 | ambiguty2 | ambiguty3 | ambiguty4
  genefile$genefunc <- ""
  genefile[ambiguties,]$genefunc <- "??"
  # 30 ??
  #table(genefile[!ambiguties,]$foundby)
  #ara both 
  #36 3525 
  genefile[!ambiguties,]$genefunc <-  genefile[!ambiguties,]$arafunc
  
  outputpath <- "/home/fatemeh/Leishmania_2019/Leishmania_2019/Document/"
  create.summary.table(genefile,outputpath)

  genefile
}
create.summary.table <- function(genefile, outputpath) {
  # create a table with columns: Annotation, Intersection, ARAonly, Union
  # with rows
  # 23 functional class

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
