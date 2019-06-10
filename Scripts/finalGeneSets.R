library(gdata)
library(readr)
library(Hmisc)
library(stringr)
library(ggplot2)
library(ggrepel)
library(gridExtra)

finalgenesets <- function(){
  
  # preparing final gene set for Leishmania paper
  # LeishGeneDF <- prepare.LeishGeneset()
  tsfmGeneDF <- prepare.tsfmGeneset()
  #tsfmGeneDF[tsfmGeneDF$tsenote!="notfound",]$genefunc <- paste(tsfmGeneDF[tsfmGeneDF$tsenote!="notfound",]$genefunc,"*",sep = "")
  #tsfmGeneDF[tsfmGeneDF$foundby=="tse",]$genefunc <- paste(tsfmGeneDF[tsfmGeneDF$foundby=="tse",]$genefunc,"(t)",sep = "")
  tsfmGeneDF[tsfmGeneDF$foundby=="ara",]$genefunc <- paste(tsfmGeneDF[tsfmGeneDF$foundby=="ara",]$genefunc,"(a)",sep = "")
  vague <- (tsfmGeneDF$foundby=="both") & (tsfmGeneDF$arafunc != tsfmGeneDF$tsefunc)
  tsfmGeneDF[vague,]$genefunc <- paste("(",tolower(tsfmGeneDF[vague,]$tsefunc),"|",tolower(tsfmGeneDF[vague,]$arafunc),")",sep = "")
  tem <- read.table("/home/fatemeh/Leishmania_2019/Leishmania_2019/Results/Integrated_Genes/undetDF2.txt",header = FALSE,skip = 1,fill = TRUE)
  IDto <- tem$V2
  IDfrom <- tem$V4
  shifts <- tsfmGeneDF$geneid %in% tem$V1
  valid <- tem$V1 %in% tsfmGeneDF$geneid
  tsfmGeneDF[shifts,]$genefunc <- paste("(",tolower(IDto[valid]),"<",tolower(IDfrom[valid]),")",sep = "")
  clusterDF <- create.clusterDF(tsfmGeneDF)
  
  
  clusterDF$clusterSeq2 <- paste(clusterDF$clusterSeq2,clusterDF$clusterDir,sep = ",")
  
  # TrypDF <- tsfmGeneDF[tolower(substr(tsfmGeneDF$sourceOrg,1,1))!="l",]
  # LeishDF <- tsfmGeneDF[tolower(substr(tsfmGeneDF$sourceOrg,1,1))!="t",]
  # clusterDF_Tryp <- create.clusterDF(TrypDF)
  # clusterDF_Leish <- create.clusterDF(LeishDF)

  # 
  # 
  # group1 <- c("TbruceigambienseDAL972","TbruceiLister427","TbruceiTREU927", "TevansiSTIB805")
  # group2 <- c("TcruziCLBrenerNon-Esmeraldo-like", "TcruziCLBrenerEsmeraldo-like")
  # group3 <- c("LaethiopicaL147","LarabicaLEM1108","LbraziliensisMHOMBR75M2903","LbraziliensisMHOMBR75M2904","LdonovaniBHU1220",
  #             "LdonovaniBPK282A1","LenriettiiLEM3045","LgerbilliLEM452","LinfantumJPCM5","LmajorFriedlin","LmajorLV39c5","LmajorSD75","LmexicanaMHOMGT2001U1103",
  #             "LpanamensisMHOMCOL81L13","LpanamensisMHOMPA94PSC1","LspMARLEM2494","LtropicaL590","LturanicaLEM423")
  # 
  # group4 <- c("LpanamensisMHOMCOL81L13","LpanamensisMHOMPA94PSC1","LbraziliensisMHOMBR75M2904","LbraziliensisMHOMBR75M2903")#?????????
  # group5 <- c("LamazonensisMHOMBR71973M2269","LmexicanaMHOMGT2001U1103")
  # group6 <- c("LinfantumJPCM5", "LdonovaniBPK282A1","LdonovaniBHU1220","LaethiopicaL147","LtropicaL590","LgerbilliLEM452","LarabicaLEM1108","LturanicaLEM423","LmajorFriedlin","LmajorLV39c5")
   group7 <- c("TcruzicruziDm28c","TcruziDm28c","TcruziJRcl4","TcruziSylvioX10-1","TcruziSylvioX10-1-2012")
  # group8 <- c("TcruziCLBrenerNon-Esmeraldo-like", "TcruziCLBrenerEsmeraldo-like","TcruziCLBrener")
  # group8 <- c("TcruzimarinkelleiB7","TcruziEsmeraldo","TcruziTulacl2")
  #group9 <- c("LpyrrhocorisH10","LseymouriATCC30220","TvivaxY486","TgrayiANR4","TcongolenseIL3000","PconfusumCUL13","LtarentolaeParrotTarII","LspMARLEM2494","LenriettiiLEM3045" ,"EmonterogeiiLV88","BayalaiB08-376", "CfasciculataCfCl")
  # 
  G1 <- c("LgerbilliLEM452","LarabicaLEM1108","LmajorFriedlin","LmajorLV39c5","LturanicaLEM423")
  G2 <- c("LaethiopicaL147","LtropicaL590")
  G3 <- c("LinfantumJPCM5", "LdonovaniBPK282A1","LdonovaniBHU1220")
  #G4 <- c("LmexicanaMHOMGT2001U1103","LamazonensisMHOMBR71973M2269")
  G5 <- c("LpanamensisMHOMCOL81L13","LpanamensisMHOMPA94PSC1","LbraziliensisMHOMBR75M2904","LbraziliensisMHOMBR75M2903")
  #G6 <- c("LtarentolaeParrotTarII","LenriettiiLEM3045")
  #G7 <- c("LpyrrhocorisH10","LseymouriATCC30220","CfasciculataCfCl")
  G8 <- c("TcruziCLBrenerNon-Esmeraldo-like", "TcruziCLBrenerEsmeraldo-like")
  G9 <- c("TbruceigambienseDAL972","TbruceiLister427","TbruceiTREU927", "TevansiSTIB805")
  summerizedVisualization(clusterDF,G1,"G1")
  summerizedVisualization(clusterDF,G2,"G2")
  summerizedVisualization(clusterDF,G3,"G3")
  #summerizedVisualization(clusterDF,G4,"G4")
  summerizedVisualization(clusterDF,G5,"G5")
  #summerizedVisualization(clusterDF,G6,"G6")
  #summerizedVisualization(clusterDF,G7,"G7")
  summerizedVisualization(clusterDF,G8,"G8")
  summerizedVisualization(clusterDF,G9,"G9")
  
  # summerizedVisualization(clusterDF,group1,"Group1")
  # summerizedVisualization(clusterDF,group2,"Group2")
  # summerizedVisualization(clusterDF,group3,"Group3")
  # summerizedVisualization(clusterDF,group4,"Group4")
  # summerizedVisualization(clusterDF,group5,"Group5_")
  # summerizedVisualization(clusterDF,group6,"Group6")
  # summerizedVisualization(clusterDF,group7,"Group7")
  # summerizedVisualization(clusterDF,group8,"Group8")
  # summerizedVisualization(clusterDF,group9,"Group9")
  
  visializeAllClusters(clusterDF_Tryp)
  visializeAllClusters(clusterDF_Leish)
  
  
  # check which one of the genes finders has highest score and what is the function for the highest score 
  # make a df with columns: "tsefunc", "tseScore", "arafunc","araScore","maxFunc","maxScore"
  # maxfun <- character(length = nrow(scoresdf_train))
  # maxarr <- numeric(length = nrow(scoresdf_train))
  # tsescores <- numeric(length = nrow(scoresdf_train))
  # tsefunc <- character(length = nrow(scoresdf_train))
  # arascore <- numeric(length = nrow(scoresdf_train))
  # arafunc <- character(length = nrow(scoresdf_train))
  # for (i in 1:nrow(scoresdf_train)) {
  #   maxfun[i] <- names(scoresdf_train)[scoresdf_train[i, ] == max(scoresdf_train[i, ])]
  #   maxarr[i] <- max(scoresdf_train[i, ])
  #   tsefunc[i] <- vagueDF[vagueDF$geneid == undetDF$geneid[i],]$tsefunc
  #   arafunc[i] <- vagueDF[vagueDF$geneid == undetDF$geneid[i],]$arafunc
  #   #tsescore[i] <- scoresdf_train[i,]$tsefunc[i]
  #   #arascore[i] <- scoresdf_train[i,]$arafunc[i]
  # }
  # undetinfo <- data.frame(maxfun,maxarr,tsefunc,arafunc)
  

  
  # plot for each big cluster
  visualizeBigclusters(clusterDF)
  visializeAllClusters(clusterDF)
  # isbig <- nchar(clusterDF$clusterSeq2) > 4
  # bigclusterDF <- clusterDF[isbig,]
  # clusSeqs <- names(table(bigclusterDF$clusterSeq2))
  # ClusCountsDF <- as.data.frame(table(bigclusterDF$clusterSeq2))
  # names(ClusCountsDF) <- c("cluster","repeats")
  # ClusCountsDF$Organisms <- ""
  # for (i in 1:length(clusSeqs)) {
  #   currentClus <- bigclusterDF$clusterSeq2 == clusSeqs[i]
  #   ClusCountsDF$Organisms[i] <- paste(bigclusterDF[currentClus,]$sourceOrg,collapse = ",")
  # }
  # write.fwf(
  #   ClusCountsDF,
  #   file = "/home/fatemeh/Leishmania_2019/Leishmania_2019/Document/Figures/HomoClusters2/ClusCountsDF.txt",
  #   #sep = "\n",
  #   colnames = FALSE
  # )
  group1 <- c("TbruceigambienseDAL972","TbruceiLister427","TbruceiTREU927", "TevansiSTIB805")
  group2 <- c("TcruziCLBrenerNon-Esmeraldo-like", "TcruziCLBrenerEsmeraldo-like")
  group3 <- c("LaethiopicaL147","LarabicaLEM1108","LbraziliensisMHOMBR75M2903","LbraziliensisMHOMBR75M2904","LdonovaniBHU1220",
              "LdonovaniBPK282A1","LenriettiiLEM3045","LgerbilliLEM452","LinfantumJPCM5","LmajorFriedlin","LmajorLV39c5","LmajorSD75","LmexicanaMHOMGT2001U1103",
              "LpanamensisMHOMCOL81L13","LpanamensisMHOMPA94PSC1","LspMARLEM2494","LtropicaL590","LturanicaLEM423")
  
  summerizedVisualization(clusterDF,group1,"Group1")
  summerizedVisualization(clusterDF,group2,"Group2")
  summerizedVisualization(clusterDF,group3,"Group3")
  
  summerizedVisualization(clusterDF,group1,"Group11")
  summerizedVisualization(clusterDF,group2,"Group22")
  summerizedVisualization(clusterDF,group3,"Group33")
  
}

summerizedVisualization <- function(clusterDF,group1,filen){
  isbig <- nchar(clusterDF$clusterSeq2) > 0
  bigclusterDF <- clusterDF[isbig,]
  ### ??????????????????????????? order does not work if I have the name of the dataframe ?????
  bigclusterDF <- bigclusterDF[order(bigclusterDF$sourceOrg,bigclusterDF$soureSeq,bigclusterDF$clusterBegin),]
  # assign IDs to each cluster of same sourceOrg
  bigclusterDF$LocalClusID <- 0
  bigclusterDF$LocalClusID[1] <- 1
  
  bigclusterDF$seqnum <- 0
  bigclusterDF$seqnum[1] <- 1
  seqnumcount <- 1
  setnumber <- 1
  for (i in 1:(nrow(bigclusterDF)-1)) {
    
   
    
    if(bigclusterDF$sourceOrg[i+1]==bigclusterDF$sourceOrg[i])
    {
      setnumber <- setnumber + 1
      bigclusterDF$LocalClusID[i+1] <- setnumber
    }
    else
    {
      seqnumcount <- 1
      setnumber <- 1
      bigclusterDF$LocalClusID[i+1] <- setnumber
    }
    
    
    if(bigclusterDF$soureSeq[i + 1] != bigclusterDF$soureSeq[i])
    {
      seqnumcount <- seqnumcount + 1
      bigclusterDF$seqnum[i+1] <-  seqnumcount
    }
    else
    {
      bigclusterDF$seqnum[i+1] <- seqnumcount
    }
  }
  
  
  #############
  isgroup1 <- bigclusterDF$sourceOrg %in% group1
  bigclusterDF <- bigclusterDF[isgroup1,]
  ###########
  sources <- nrow(table(bigclusterDF$sourceOrg))
  bigclusterDF$sourceOrg <- factor(bigclusterDF$sourceOrg)
  bigclusterDF$isbig <- FALSE
  # for (i in 1:nrow(bigclusterDF)) {
  #   bigclusterDF$isbig[i] <- nchar(bigclusterDF$clusterSeq2[i]) > 4
  # }
  for (i in 1:nrow(bigclusterDF)) {
    bigclusterDF$isbig[i] <- all(grepl("[[:upper:]]", strsplit(bigclusterDF$clusterSeq2[i], "")[[1]]))
  }
  
  bigclusterDF$col <- "0"
  bigclusterDF[(bigclusterDF$seqnum %% 2) == 0,]$col <- "1"
  #bigclusterDF[!bigclusterDF$isbig,]$col <- "2"
  g <- ggplot(
    data = bigclusterDF,
    aes(
      x = bigclusterDF$LocalClusID,
      y = bigclusterDF$sourceOrg,
      label = bigclusterDF$clusterSeq2,
      color = bigclusterDF$col
    )
  ) + geom_point() + theme(
    legend.position = "none") + geom_label_repel(
      aes(label = bigclusterDF$clusterSeq2),
      box.padding   = 0.1,
      point.padding = 0.1,
      segment.color = 'grey50',force = 1,nudge_x = 0.1
    )
 
  myheight <- sources * 2
  ggsave(file = paste("/home/fatemeh/Leishmania_2019/Leishmania_2019/Document/Figures/clusters_nofilter/",filen,".jpg",sep = ""), g,width = 120, height = myheight, units = "cm",limitsize = FALSE)
  
  
}
visializeAllClusters <- function(clusterDF){
  isbig <- nchar(clusterDF$clusterSeq2) > 0
  bigclusterDF <- clusterDF[isbig,]
  bigclusterDF <- bigclusterDF[order(bigclusterDF$sourceOrg,bigclusterDF$soureSeq,bigclusterDF$clusterBegin),]
  # assign IDs to each cluster of same sourceOrg
  bigclusterDF$LocalClusID <- 0
  bigclusterDF$LocalClusID[1] <- 1
  setnumber <- 1
  for (i in 1:(nrow(bigclusterDF)-1)) {
    if(bigclusterDF$sourceOrg[i+1]==bigclusterDF$sourceOrg[i])
    {
      setnumber <- setnumber + 1
      bigclusterDF$LocalClusID[i+1] <- setnumber
    }
    else
    {
      setnumber <- 1
      bigclusterDF$LocalClusID[i+1] <- setnumber
    }
  }
  #bigclusterDF <- bigclusterDF[order(bigclusterDF$clusterSeq2),]
  bigclusterDF$sourceOrg <- factor(bigclusterDF$sourceOrg)
  # Y axis are sourceOrg, X is the clusters, points are labeled with the clusterseq
  
  for (i in 1:nrow(bigclusterDF)) {
  bigclusterDF$islower[i] <- all(grepl("[[:upper:]]", strsplit(bigclusterDF$clusterSeq2[i], "")[[1]]))
  }
  # before plotting sort the sourceOrg based on the 
  # bigclusterDF$isbig <- FALSE
  # for (i in 1:nrow(bigclusterDF)) {
  #   bigclusterDF$isbig[i] <- nchar(bigclusterDF$clusterSeq2[i]) > 5 
  # }
  g <- ggplot(
    data = bigclusterDF,
    aes(
      x = bigclusterDF$LocalClusID,
      y = bigclusterDF$sourceOrg,
      label = bigclusterDF$clusterSeq2,
      color = bigclusterDF$islower
    )
  ) + geom_point() + theme(
    legend.position = "none") + geom_label_repel(
      aes(label = bigclusterDF$clusterSeq2),
      box.padding   = 0.1,
      point.padding = 0.1,
      segment.color = 'grey50',force = 1,nudge_x = 0.1
    )
  
  ggsave(file = "/home/fatemeh/Leishmania_2019/Leishmania_2019/Document/Figures/vagues.jpg", g,width = 120, height = 70, units = "cm",limitsize = FALSE)
  
}
visualizeBigclusters <- function(clusterDF){
  isclus <- nchar(clusterDF$clusterSeq2) > 0
  allclusterDF <- clusterDF[isclus,]
  
  isbig <- nchar(clusterDF$clusterSeq2) > 4
  bigclusterDF <- clusterDF[isbig,]
  clusSeqs <- names(table(bigclusterDF$clusterSeq2))
  m <- matrix(nrow = length(clusSeqs), ncol =  length(clusSeqs))
  distanceDF <- as.data.frame(m)
  for (i in 1:length(clusSeqs)) {
    for (j in 1:length(clusSeqs)) {
      distanceDF[i, j] <- adist(clusSeqs[i], clusSeqs[j])
    }
  }
  hc <- hclust(as.dist(distanceDF), method='ward.D2')
  plot(hc, cex = 0.6, hang = -1)
  sub_groups <- cutree(hc, h = 1)
  sub_groups[9] <- 8
  sub_groups[11] <- 10
  sub_groups[15] <- 13
  for (i in names(table(sub_groups))){
    i <- as.integer(i)
    clusterSeq_arr <- clusSeqs[sub_groups==sub_groups[i]]
    # plot those sourceOrg that have this cluster
    sources<- names(table(allclusterDF[allclusterDF$clusterSeq2 %in% clusterSeq_arr,]$sourceOrg))
    if(length(sources) == 1)
      next
    mydf <- allclusterDF[allclusterDF$sourceOrg %in% sources,]
    mydf <- mydf[order(mydf$sourceOrg,mydf$soureSeq,mydf$clusterBegin),]
    # assign IDs to each cluster of same sourceOrg
    mydf$LocalClusID <- 0
    mydf$LocalClusID[1] <- 1
    setnumber <- 1
    for (j in 1:(nrow(mydf)-1)) {
      if(mydf$sourceOrg[j+1]==mydf$sourceOrg[j])
      {
        setnumber <- setnumber + 1
        mydf$LocalClusID[j+1] <- setnumber
      }
      else
      {
        setnumber <- 1
        mydf$LocalClusID[j+1] <- setnumber
      }
    }
    mydf$sourceOrg <- factor(mydf$sourceOrg)

    mydf$mainclus <- FALSE
    mydf$mainclus <- mydf$clusterSeq2 %in% clusterSeq_arr
    
    mydf$
    g <- ggplot(
      data = mydf,
      aes(
        x = mydf$LocalClusID,
        y = mydf$sourceOrg,
        label = mydf$clusterSeq2,
        color = mydf$mainclus
      )
    ) + geom_point() + theme(
      legend.position = "none",axis.text.y = element_text(face = "bold")) + geom_label_repel(
        aes(label = mydf$clusterSeq2),
        box.padding   = 0.1,
        point.padding = 0.1,
        segment.color = 'grey50',force = 1,nudge_x = 0.1
      ) + xlab("tRNA gene clusters in order")+ ylab("Organism")
    filename <- paste("/home/fatemeh/Leishmania_2019/Leishmania_2019/Document/Figures/HomoClusters2/",i,".jpg",sep = "")
    myheight <- length(sources) * 2
    ggsave(filename,g,width = 50, height = myheight, units = "cm")
    
  }
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
  dismiss1 <- (genefile$foundby!="ara") & (genefile$arascore < 107 | genefile$tsescore < 50) # most of these genes were pseudo or truncated or both
  dismiss2 <- (genefile$foundby=="ara") & (genefile$arascore < 107)
  genefile1 <- genefile[(!dismiss1 & !dismiss2),]
  
  dismiss3 <- genefile1$tsefunc == "" &  genefile1$arafunc == ""
  genefile2 <- genefile1[!dismiss3, ]
  #table( genefile2$tsenote)
  #notfound truncated 
  #3586         2
  # remove 6 truncated|pseudo genes for tsfm input
  nonpseudo_nontrunc <- genefile2$tsenote == "notfound"
  genefile2 <- genefile2[nonpseudo_nontrunc, ]
  
  #genefile <- genefile3
  
  #vaguedf <- data.frame(genefile[vage,]$tsegeneseq,genefile[vage,]$geneid)
  #names(vaguedf) <- c("sequences","geneid")
  # bothdf <- genefile[genefile$foundby == "both",]
  # vage <- bothdf$tsefunc != bothdf$arafunc
  #   
  # vaguegeneID <- trainDF_intersection$geneid %in% bothdf[vage,]$geneid 
  # 
  # trainDF_intersection[vaguegeneID,]
  # vaguescoreDF <- calculate.score(pos_neg_list, trainDF_intersection[vaguegeneID,], 0)
  # maxscore <- numeric(length = ncol(vaguescoreDF))
  # minfunc <- character(length = ncol(vaguescoreDF))
  # for (i in 1:nrow(vaguescoreDF)) {
  #   maxscore[i] <- max(vaguescoreDF[i,])
  #   minfunc[i] <- names(vaguescoreDF)[vaguescoreDF[i,]== maxscore[i]]
  # }
  # temp <- bothdf[vage,]$geneid %in% trainDF_intersection$geneid
  # table( minfunc == bothdf[vage,][temp,]$tsefunc)
  # bothdf[vage,][temp,][minfunc != bothdf[vage,][temp,]$tsefunc,]
  # _________________________________________________________________________________________
  # Add a column begin, end and genecode to the dataframe to be able to assign clusters later
  genefile <- genefile2
  genefile$begin <- ""
  genefile$end <- ""
  genefile$genecode <- ""
  genefile$genefunc <- ""
  
  genefile[genefile$foundby == "both" |
             genefile$foundby == "tse", ]$begin <-
    genefile[genefile$foundby == "both" |
               genefile$foundby == "tse", ]$tsebegin
  genefile[genefile$foundby == "ara", ]$begin <-
    genefile[genefile$foundby == "ara", ]$arabegin
  genefile[genefile$foundby == "both" |
             genefile$foundby == "tse", ]$end <-
    genefile[genefile$foundby == "both" |
               genefile$foundby == "tse", ]$tseend
  genefile[genefile$foundby == "ara", ]$end <-
    genefile[genefile$foundby == "ara", ]$araend
  
  genefile[genefile$foundby == "both" |
             genefile$foundby == "tse", ]$genecode <-
    paste("(",tolower(genefile[genefile$foundby == "both" |
               genefile$foundby == "tse", ]$tseac),")",genefile[genefile$foundby == "both" |
                                                               genefile$foundby == "tse", ]$tsefunc,sep="")
  genefile[genefile$foundby == "ara", ]$genecode <-
    paste(genefile[genefile$foundby == "ara", ]$araac,genefile[genefile$foundby == "ara", ]$arafunc,sep="")
  
  genefile[genefile$foundby == "both" |
             genefile$foundby == "tse", ]$genefunc <-
    genefile[genefile$foundby == "both" |
               genefile$foundby == "tse", ]$tsefunc
  genefile[genefile$foundby == "ara", ]$genefunc <-
    genefile[genefile$foundby == "ara", ]$arafunc
  
  genefile$begin <- as.integer(genefile$begin)
  genefile$end <- as.integer(genefile$end)
  # ______________________________________________________________________________________
  
 
  #table(genefile3$foundby)
  #ara both  tse 
  #36 3571    1 
  #table(genefile3[genefile3$foundby=="both",]$arafunc!=genefile3[genefile3$foundby=="both",]$tsefunc)
  #FALSE  TRUE 
  #3538    33 
  # lets implemet the the method for finding similar clusters first
  
  genefile
}
create.clusterDF <- function(tsfmGeneDF){
  geneset <- assignCluster(tsfmGeneDF,1000)
  m <- matrix(nrow = max(geneset$cluster),ncol = 9)
  clusterDF<- as.data.frame(m)
  names(clusterDF) <- c("sourceOrg","soureSeq","clusterID","clusterSeq1","clusterSeq2", "clusterBegin","clusterEnd","clusterDir","clusterLen")
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
  }
  clusterDF
}
geneClusterVisualization <- function() {
  # read gene file
  geneDF <- read.input()
  geneDF$begin <- as.integer(geneDF$begin)
  geneDF$end <- as.integer(geneDF$end)
  clusterdistance <- 5000
  geneDF <- assignCluster(geneDF, clusterdistance)
  geneOrganization(geneDF)
  geneScore(geneDF)
  
  ara_geneDF <- geneDF[geneDF$foundby == "ara",]
}
read.input <- function() {
  # this function will read integrated_Tse_Ara gene file (output of Integrate_Tse_Ara script)
  # Dismiss genes with no identity assigned to them (15 genes) the detected anticodon for these were either 2bp or NNN
  # Dismiss genes that are marked as pseudo or truncated by tse. (39 genes 17 of them found only by tse and 8 found by both ara and tse)
  # for genes in intersection of two genefinders with same predicted identity, add the model at the end of gene id. (3560 genes)
  # for those with different predicted identity it will add _undet at the end of geneid (33  genes)
  # for genes found by ony tse, it will add _tse<genemodel> at the end of gene id. (6)
  # for genes found by ony tse, it will add _ara<genemodel> at the end of gene id. (741)
  # keep all the information needed to make a fasta file for doing the alignment
  resultpath <-
    "/home/fatemeh/Leishmania_2019/Leishmania_2019/Results/Integrated_Genes/"
  genefilepath <-
    "/home/fatemeh/Leishmania_2019/Leishmania_2019/Results/Integrated_Genes/integrated_tse_ara.txt"
  genefile <-
    read.table(genefilepath, header = TRUE, colClasses = "character")
  dismiss <- genefile$tsefunc == "" &  genefile$arafunc == ""
  genefile <- genefile[!dismiss, ]
  nonpseudo_nontrunc <- genefile$tsenote == "notfound"
  genefile <- genefile[nonpseudo_nontrunc, ]
  #_________________________________________________________________________________________________
  both <-
    (genefile$foundby == "both") &
    (genefile$tsefunc == genefile$arafunc)
  undet <-
    (genefile$foundby == "both") &
    (genefile$tsefunc != genefile$arafunc)
  genefile$geneid[both] <-
    paste(genefile$geneid[both], "_", genefile$tsefunc[both], sep = "")
  genefile$geneid[undet] <-
    paste(genefile$geneid[undet], "_undet", sep = "")
  #_________________________________________________________________________________________________
  tseonly <- genefile$foundby == "tse"
  genefile$geneid[tseonly] <-
    paste(genefile$geneid[tseonly], "_tse", genefile$tsefunc[tseonly], sep = "")
  #_________________________________________________________________________________________________
  araonly <- genefile$foundby == "ara"
  genefile$geneid[araonly] <-
    paste(genefile$geneid[araonly], "_ara", genefile$arafunc[araonly], sep = "")
  #_________________________________________________________________________________________________
  genefile$begin <- ""
  genefile$end <- ""
  genefile[genefile$foundby == "both" |
             genefile$foundby == "tse", ]$begin <-
    genefile[genefile$foundby == "both" |
               genefile$foundby == "tse", ]$tsebegin
  genefile[genefile$foundby == "ara", ]$begin <-
    genefile[genefile$foundby == "ara", ]$arabegin
  genefile[genefile$foundby == "both" |
             genefile$foundby == "tse", ]$end <-
    genefile[genefile$foundby == "both" |
               genefile$foundby == "tse", ]$tseend
  genefile[genefile$foundby == "ara", ]$end <-
    genefile[genefile$foundby == "ara", ]$araend
  drops <- c("tsebegin", "tseend", "arabegin", "araend")
  genefile <- genefile[, !(names(genefile) %in% drops)]
  genefile
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
#####################################################################
geneOrganization <- function(geneDF) {
  araCluster <- geneDF[geneDF$foundby == "ara",]$cluster
  # visualize these clusters
  # how many members these clusters have?
  cluster_countDF <- data.frame(table(geneDF$cluster))
  names(cluster_countDF) <- c("cluster", "genecount")
  table(cluster_countDF$genecount > 1)
  #FALSE  TRUE
  #1382   859
  # number of clusters with more than one gene
  table(cluster_countDF[cluster_countDF$cluster %in% araCluster,]$genecount > 1)
  # FALSE  TRUE
  # 680    30
  ara_cluster_countDF <-
    cluster_countDF[cluster_countDF$cluster %in% araCluster,]
  isbig <- ara_cluster_countDF$genecount > 1
  bigclusters <- ara_cluster_countDF[isbig,]
  # visualizing the big clusters
  geneDF_bigclusters <-
    geneDF[geneDF$cluster %in% bigclusters$cluster,]
  # number of genes in big clusters to be visualized
  nrow(geneDF_bigclusters)
  # 82 out of 741
  # make genecodes with Anticodon+
  geneDF_bigclusters$genecode  <-
    paste(geneDF_bigclusters$araac, geneDF_bigclusters$arafunc, sep = "")
  # visualize each sequence: y axis= sourceseq. x axis= clusters
  # sort based on the sequence and location
  geneDF_bigclusters <-
    geneDF_bigclusters[order(geneDF_bigclusters$sourceseq, geneDF_bigclusters$begin),]
  # Subtract all the locations by the mean location of each sourceseq
  for (i in 1:nrow(table(geneDF_bigclusters$sourceseq))) {
    sourceseq <- names(table(geneDF_bigclusters$sourceseq)[i])
    minloc <-
      min(geneDF_bigclusters[geneDF_bigclusters$sourceseq == sourceseq,]$begin)
    geneDF_bigclusters[geneDF_bigclusters$sourceseq == sourceseq,]$begin <-
      geneDF_bigclusters[geneDF_bigclusters$sourceseq == sourceseq,]$begin - minloc
    geneDF_bigclusters[geneDF_bigclusters$sourceseq == sourceseq,]$end <-
      geneDF_bigclusters[geneDF_bigclusters$sourceseq == sourceseq,]$end - minloc
  }
  
  
  
  myclusters <- levels(factor(geneDF_bigclusters$cluster))
  listofplots <-
    lapply(myclusters, function(i) {
      print(i)
      mydata = geneDF_bigclusters[geneDF_bigclusters$cluster == i,]
      ggplot(
        data = mydata,
        aes(
          x = mydata$begin,
          y = mydata$sourceseq,
          label = mydata$genecode,
          color = mydata$foundby
        )
      ) + geom_point() + theme(
        legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_text(angle = 90),
        plot.margin = unit(c(0, 0, 0, 0), "cm"),
        title = mydata$sourceOrg[1]
      ) + geom_label_repel(
        aes(label = mydata$genecode),
        box.padding   = 0.35,
        point.padding = 0.5,
        segment.color = 'grey50'
      )
    })
  
  print(length(listofplots))
  g <- do.call(arrangeGrob, c(listofplots))
  ggsave(file = "/home/fatemeh/big_ara_cluster.jpg", g)
  
}
geneScore <- function(geneDF) {
  resultpath <-
    "/home/fatemeh/Leishmania_2019/Leishmania_2019/Results/Integrated_Genes/"
  genefilepath <-
    "/home/fatemeh/Leishmania_2019/Leishmania_2019/Results/Integrated_Genes/integrated_tse_ara.txt"
  geneDF <-
    read.table(genefilepath, header = TRUE, colClasses = "character")
  geneDF2 <-
    geneDF[geneDF$foundby == "ara" | geneDF$foundby == "both",]
  geneDF2$arascore <- as.integer(geneDF2$arascore)
  geneDF2$Foundby <- "tse"
  geneDF2[geneDF2$foundby!="ara",]$Foundby <- "both"
  geneDF2[geneDF2$tsenote!="notfound",]$Foundby <- "both: \n pseudo or \n truncated"
  #"Aragorn score of genes found by both genefinders TSE and ARA, and genes found by only ARA"
  ggplot(geneDF2,
         aes(
           x = geneDF2$geneid,
           y = geneDF2$arascore,
           color = Foundby
         )) + theme_bw() +
    geom_point() + xlab("Genes") + ylab("Aragorn Score") + ggtitle(
      ""
    )  +
    theme(
      text = element_text(size = 25),
      axis.text.x = element_blank(),
      plot.title = element_text(
        color = "black",
        size = 12,
        hjust = 0.5
      ),
      axis.ticks = element_blank()
    ) + scale_y_continuous(breaks =seq(min(geneDF2$arascore),max(geneDF2$arascore),1))
  
  # plot for tse scores
  geneDF3 <-
    geneDF[geneDF$foundby == "tse" | geneDF$foundby == "both",]
  geneDF3$tsescore <- as.integer(geneDF3$tsescore)
  #tRNAscan-SE score of genes found by both genefinders TSE and ARA, and genes found by only TSE
  ggplot(geneDF3,
         aes(
           x = geneDF3$geneid,
           y = geneDF3$tsescore,
           color = foundby
         )) + theme_bw() +
    geom_point() + xlab("Genes") + ylab("TSE Score") + ggtitle(
      ""
    )  +
    theme(
      text = element_text(size = 25),
      axis.text.x = element_blank(),
      plot.title = element_text(
        color = "black",
        size = 12,
        hjust = 0.5
      ),
      axis.ticks = element_blank()
    ) 
}