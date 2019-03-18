# This function will read the outlier.fasta file and subtract the from the intersection gene file. 
# Make a model for each functional class and saves the results in a list of dataframes called profile_list
# At the end will also look at the distribution of the scores for each model and make sure that they are normally distributed 
classifier <- function() {
  genedirpath <-
    "/home/fatemeh/Leishmania_2019/Leishmania_2019/Results/Integrated_Genes/TriTryp_EditedCovea.fasta"
  outlierpath <-
    "/home/fatemeh/Leishmania_2019/Leishmania_2019/Results/Integrated_Genes/TriTrypgenemodel_Intersection/TriTryp_EditedCovea/outlier/outliers.txt"
  outputdir <-
    "/home/fatemeh/Leishmania_2019/Leishmania_2019/Results/Integrated_Genes/"
  trainDF <- read.genefile(genedirpath) # 3543 genes
  # gaps are removed from outier sequences by OD-seq package
  outlierDF <- read.genefile(outlierpath) # 39 outier
  trainDF <- remove.outliers(trainDF, outlierDF) # 3504 genes left
  trainDF <- make.clusters(trainDF)
  # make the 22 models
  #trainDF <- trainDF[trainDF$cluster == "T", ]
  profile_list <- profile.cm(trainDF)
  # print(paste("we have", length(names(profile_list)), "functional classes:"))
  # print(names(profile_list))
  # use the trainDF as the undet class for now
  scoresdf <- find.seqScore(profile_list, trainDF)
  score.dist(scoresdf, trainDF$funclass)
  
  zscores <- find.Zscore(scoresdf, trainDF$funclass)
  zscores <- round(zscores)
  trainDF
}

make.clusters <- function(trainDF){
  trainDF$cluster <- ""
  for (i in 1:nrow(trainDF)) {
    firstletter <- substr(trainDF$geneid[i],1,1)
    if(firstletter=="L")
      trainDF$cluster[i] <- "L"
    if(firstletter=="T")
      trainDF$cluster[i] <- "T"
  }
  trainDF
}
#__________________________________________________________________________________________
remove.outliers <- function(trainDF,outlierDF){
  isoutlier <- trainDF$geneid %in% outlierDF$geneid
  trainDF <- trainDF[!isoutlier,] 
}
#__________________________________________________________________________________________
read.genefile <- function(genedirpath){
  # this function read the gene file into a dataframe as teh training set with columns: "header" "sequences" "funclass" "geneid"
  # at the end it will remove sequences of undet model
  coveadf <-
    read.table(genedirpath)
  coveaArr <- as.character(coveadf$V1)
  headers <- coveaArr[seq(1, length(coveaArr), 2)]
  sequences <- coveaArr[seq(2, length(coveaArr), 2)]
  funclass <- character(length = length(headers))
  geneid <- character(length = length(headers))
  for (i in 1:length(headers)) {
    temp <- unlist(strsplit(headers[i], split = "_"))
    geneid[i] <-
      gsub(">", "", paste(temp[1:length(temp) - 1], collapse = "_"))
    funclass[i] <- temp[length(temp)]
  }
  GeneDF <- data.frame(headers, sequences, funclass, geneid)
  for (i in 1:ncol(GeneDF)) {
    GeneDF[, i] <- as.character(GeneDF[, i])
  }
  GeneDF <-
    GeneDF[(GeneDF$funclass != "undet"), ]# & (GeneDF$cluster == "T"), ]
  GeneDF
}
#__________________________________________________________________________________________
score.dist <- function(scoresdf,trainDF_Fun){

  par(mfrow = c(4, 6))
  for (i in 1:ncol(scoresdf)) {
    samemodel <- names(scoresdf)[i] == trainDF_Fun
    currentmodel_score <- scoresdf[samemodel, i]
    hist(
      currentmodel_score, 
      prob=TRUE, 
      main = "",
      xlab = paste(
        "score against model",
        names(scoresdf)[i],
        sep = " "),
        ylab = paste("freq (total=", length(currentmodel_score), ")", sep = "")
      )
    lines(density(currentmodel_score), col="blue", lwd=2)
    #plot(density(currentmodel))
  }
  #par(mar=c(0.5, 0.5, 0.5, 0.5))
  # plotlist <- list()
  # count = 1
  # filenm <- paste("/home/fatemeh/Leishmania_2019/Leishmania_2019/Document/","modeldistribution", ".pdf", sep = "")
  # for (i in 1:ncol(scoresdf)) {
  #   samemodel <- names(scoresdf)[i] == trainDF_Fun
  #   currentmodel <- scoresdf[samemodel, i]
  #   currentmdf <- data.frame(currentmodel,seq(1,length(currentmodel),1))
  #   names(currentmdf) <- c("score","index")
  #   p <- ggplot(currentmdf, aes(x=currentmdf$score)) + 
  #     geom_histogram(aes(y=..density..),  
  #                    binwidth=.5,
  #                    colour="black", fill="white") +
  #     geom_density(alpha=.2, fill="#FF6666") 
  #   plotlist[[count]] <-  p
  #   count <- count + 1
  # }
  # ml <- do.call("grid.arrange", c(plotlist, ncol = 4))
  # ggsave(ml,
  #        width = 20,
  #        height = 30,
  #        filename = "")
}
# ___________________________________________________________________________________________
find.anticodon <- function(undetDF){
  genefilepath <-
    "/home/fatemeh/Leishmania_2019/Leishmania_2019/Results/Integrated_Genes/integrated_tse_ara.txt"
  geneDF <-
    read.table(genefilepath, header = TRUE, colClasses = "character")
  undetDF$tseac <- ""
  undetDF$araac <- ""
  for (i in 1:nrow(undetDF)) {
    undetDF$tseac[i] <- geneDF[undetDF$geneid[i]==geneDF$geneid,]$tseac
    undetDF$araac[i] <- geneDF[undetDF$geneid[i]==geneDF$geneid,]$araac
  }
  undetDF
}
# ___________________________________________________________________________________________
find.Zscore <- function(scoresdf,trainDF_Fun){
  stats <-
    data.frame(names(scoresdf), rep(0, ncol(scoresdf)), rep(0, ncol(scoresdf)))
  names(stats) <- c("FUNC", "mean", "std")
  for (i in 1:ncol(scoresdf)) {
    samemodel <- names(scoresdf)[i] == trainDF_Fun
    stats[i, 2] <-  mean(scoresdf[samemodel, i])
    stats[i, 3] <- sd(scoresdf[samemodel, i])
  }
  
  # make a copy of scoresdf and instead of each score calculate the z-score
  zscores <- scoresdf
  for (i in 1:ncol(scoresdf)) {
    zscores[, i] <- (scoresdf[, i] - stats[i, 2]) /  stats[i, 3]
  }
  zscores
}
#___________________________________________________________________________________________
update.genefiles <- function() {
  resultpath <-
    "/home/fatemeh/Leishmania_2019/Leishmania_2019/Results/Integrated_Genes/"
  genefilepath <-
    "/home/fatemeh/Leishmania_2019/Leishmania_2019/Results/Integrated_Genes/integrated_tse_ara.txt"
  geneDF <-
    read.table(genefilepath, header = TRUE, colClasses = "character")
  geneDF$Func <- ""
  geneDF$Func[geneDF$tsefunc == geneDF$arafunc] <-
    geneDF$tsefunc[geneDF$tsefunc == geneDF$arafunc]
  
  undetClassScorePath <-
    paste(resultpath, "UndetClassScore.txt", sep = "")
  undetClassScoreDF <-
    read.table(undetClassScorePath, header = TRUE)
  
  for (i in 1:ncol(undetClassScoreDF)) {
    undetClassScoreDF[, i] <- as.character(undetClassScoreDF[, i])
  }
  
  # Update the EditCovea.fasta file
  coveaDF <-
    read.table(paste(dirpath, "EditedCovea.fasta", sep = ""))
  coveaArr <- as.character(coveaDF$V1)
  headers <- coveaArr[seq(1, length(coveaArr), 2)]
  sequences <- coveaArr[seq(2, length(coveaArr), 2)]
  coveageneid <- headers
  for (i in 1:length(sequences)) {
    temp <- unlist(strsplit(headers[i], split = "_"))
    if (substring(headers[i], nchar(headers[i]), nchar(headers[i])) == "_")
    {
      headers[i] <- paste(headers[i], "undet", sep = "")
    }
    temp <- unlist(strsplit(headers[i], split = "_"))
    coveageneid[i] <-
      gsub(">", "", paste(temp[1:length(temp) - 1], collapse = "_"))
  }
  
  coveaDF <- data.frame(coveageneid, headers, sequences)
  
  for (i in 1:ncol(coveaDF)) {
    coveaDF[, i] <- as.character(coveaDF[, i])
  }
  
  for (i in 1:nrow(undetClassScoreDF)) {
    FUNC <- "undet"
    if (undetClassScoreDF$TseFuncScore[i] >  undetClassScoreDF$AraFuncScore[i])
      FUNC = undetClassScoreDF$TseFunc[i]
    else if (undetClassScoreDF$TseFuncScore[i] <  undetClassScoreDF$AraFuncScore[i])
      FUNC = undetClassScoreDF$AraFunc[i]
    
    coveaDF[coveaDF$coveageneid == undetClassScoreDF$geneid[i],]$headers <-
      paste(">", coveaDF[coveaDF$coveageneid == undetClassScoreDF$geneid[i],]$coveageneid, "_", FUNC, sep = "")
    geneDF[geneDF$geneid == undetClassScoreDF$geneid[i],]$Func <-
      FUNC
    # make undet sequences to remove out of this loop
    if (FUNC == "undet")
      coveaDF <-
      coveaDF[!(coveaDF$coveageneid == undetClassScoreDF$geneid[i]),]
    
  }
  
  write.fwf(
    data.frame(coveaDF$headers,
               coveaDF$sequences),
    file = paste(dirpath, "EditedCovea2.fasta", sep = ""),
    sep = "\n",
    colnames = FALSE
  )
  
  # edit the gene file
  write.table(
    geneDF,
    col.names = TRUE,
    file = paste(resultpath, "integrated_tse_ara2.txt", sep = "")
  )
  
}
#___________________________________________________________________________________________
find.seqScore <- function(profile_list, undetClass) {
  # this function will score each sequence against all the models using leave one out method
  # the final scores are saved in scoresdf
  m <- matrix(ncol = length(profile_list),
              nrow = nrow(undetClass),
              0)
  scoresdf <- as.data.frame(m)
  names(scoresdf) <- names(profile_list)
  for (x in 1:length(undetClass$sequences)) {
    myseq <- undetClass$sequences[x]
    types <- names(profile_list)
    scores <- integer(length = length(profile_list))
    for (j in 1:length(profile_list)) {
      currentDF <- profile_list[[j]]
      # total is sum of any column in our profile
      total <- sum(currentDF[, 1])
      # if this sequence was annotated as a sequence of model names(profile_list[j]) then total = total - 1
      if (undetClass$funclass[x] == names(profile_list[j]))
        total = total - 1
      for (i in 1:ncol(currentDF)) {
        currChar <- substring(myseq, i, i)
        if (currChar == "-")
          currChar <- "gap"
        model_freq <- currentDF[currChar, i]
        if (undetClass$funclass[x] == names(profile_list[j]))
          model_freq <- currentDF[currChar, i] - 1
        scores[j] <- scores[j] + log(model_freq / total)
        if (types[j] != names(profile_list[j]))
          print("shit!")
      }
    }
    
    #undetClass$funclass[x] <- types[scores == max(scores)]
    scoresdf[x,] <- scores
    
  }
  scoresdf
}
#___________________________________________________________________________________________
profile.cm <- function(trainDF) {
  # used psedo counts to avoid zero for cells
  # read in the aligned sequences from file EditedCovea.fasta
  func_list <- split(trainDF, f = trainDF$funclass)
  # for each Df make a dataframe to keep the profile cm
  seqlength <- nchar(trainDF$sequences[1])
  m <- matrix(nrow = 5, ncol = seqlength)
  prf <- data.frame(m)
  rownames(prf) <- c("A", "C", "G", "T", "gap")
  profile_list <- rep(list(prf), length(func_list))
  names(profile_list) <- names(func_list)
  for (i in 1:length(func_list)) {
    currentfuncdf <- func_list[[i]]
    currentfuncdf$sequences <- as.character(currentfuncdf$sequences)
    m <- matrix(nrow = 5, ncol = seqlength)
    prf <- data.frame(m)
    rownames(prf) <- c("A", "C", "G", "T", "gap")
    for (j in 1:seqlength) {
      # pseudo counts of one for each cell
      Afunc = 1
      Cfunc = 1
      Tfunc = 1
      Gfunc = 1
      gapfunc = 1
      # go over all the sequences and see how many of them are A,C,G,T in their jth position
      for (t in 1:nrow(currentfuncdf)) {
        if (substring(currentfuncdf$sequences[t], j, j) == "A")
          Afunc = Afunc + 1
        if (substring(currentfuncdf$sequences[t], j, j) == "C")
          Cfunc = Cfunc + 1
        if (substring(currentfuncdf$sequences[t], j, j) == "G")
          Gfunc = Gfunc + 1
        if (substring(currentfuncdf$sequences[t], j, j) == "T")
          Tfunc = Tfunc + 1
        if (substring(currentfuncdf$sequences[t], j, j) == "-")
          gapfunc = gapfunc + 1
      }
      # we do not devide them by total here to be able to do the leave one out cross validation
      # we have an array of number number of sequences in each model saved
      prf[1, j] <- Afunc # / nrow(currentfuncdf) # A
      prf[2, j] <- Cfunc # / nrow(currentfuncdf) # C
      prf[3, j] <- Gfunc # / nrow(currentfuncdf) # G
      prf[4, j] <- Tfunc # / nrow(currentfuncdf) # T
      prf[5, j] <- gapfunc # / nrow(currentfuncdf) # gap
    }
    # read the tsfm table generated for cluster Leishmania and update the model profiles
    
    # devide each column in prf by the number of sequences in the current funcclass
    # make sure sum of each column is asame
    profile_list[[i]] <- prf
    names(profile_list[i]) == names(func_list[i])
  }
  
  
  profile_list
  
}
#____________________________________________________________________________________________
