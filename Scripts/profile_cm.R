# This function will make a model for each functional class and saves the results in a list of dataframes called profile_list
# Later, we align the new undet identity genes (last part of their geneID has _undet) with the aligned file EditedCovea.fasta
# initiators are not added yet!!!
# At the end will also look at the distribution of the scores for each model and make sure that they are normally distributed 
classifier <- function() {
  # seperate undet genes from det genes as training set and test
  dirpath <-
    "/home/fatemeh/Leishmania_2019/Leishmania_2019/Results/Integrated_Genes/"
  coveadf <-
    read.table(paste(dirpath, "EditedCovea.fasta", sep = ""))
  coveaArr <- as.character(coveadf$V1)
  headers <- coveaArr[seq(1, length(coveaArr), 2)]
  sequences <- coveaArr[seq(2, length(coveaArr), 2)]
  funclass <- character(length = length(headers))
  cluster <- character(length = length(headers))
  # assign clusters 
  for (i in 1:length(funclass)) {
    #print(i)
    temp <- unlist(strsplit(headers[i], split = "_"))
    firstchar <- substr(temp[1], 2, 2)
    if (firstchar == "L")
      cluster[i] = "L"
    if (firstchar == "T")
      cluster[i] = "T"
    funclass[i] <- temp[length(temp)]
    if (substring(headers[i], nchar(headers[i]), nchar(headers[i])) == "_")
    {
      funclass[i] <- "undet"
      headers[i] <- paste(headers[i], "undet", sep = "")
    }
  }
  coveaDF <- data.frame(headers, sequences, funclass, cluster)
  for (i in 1:ncol(coveaDF)) {
    coveaDF[, i] <- as.character(coveaDF[, i])
  }
  trainDF <-
    coveaDF[(coveaDF$funclass != "undet"),]# & (coveaDF$cluster == "T"), ]
  
  # make the 22 models 
  seqlength <- nchar(sequences[1])
  profile_list <- profile.cm(trainDF, seqlength)
  
  print(paste("we have", length(names(profile_list)), "functional classes:"))
  print(names(profile_list))
  #______________________________________________________________________________________
  # use the trainDF as the undet class for now
  trainDF$geneid <- ""
  for (i in 1:nrow(trainDF)) {
    myhead <- trainDF$headers[i]
    temp <- unlist(strsplit(myhead, split = "_"))
    geneid <-
      gsub(">", "", paste(temp[1:length(temp) - 1], collapse = "_"))
    trainDF$geneid[i] <- geneid
  }
  scoresdf <- find.seqScore(profile_list, trainDF)
  
  zscores <- find.Zscore(scoresdf,trainDF$funclass)
  zscores <- round(zscores)
  scoresdf <- round(scoresdf)
  trainDF$Score1 <- -999
  trainDF$Model1 <- ""
  trainDF$Zscore <- -999
  trainDF$ModelZ <- ""
  trainDF$funcscore <- ""
  for (i in 1:nrow(scoresdf)) {
    trainDF$funcscore[i] <- scoresdf[i,names(scoresdf)==trainDF$funclass[i]]
    trainDF$Model1[i] <- names(which.max(scoresdf[i, ]))
    trainDF$Score1[i] <- scoresdf[i, which.max(scoresdf[i, ])[1]]
    trainDF$Zscore[i] <- zscores[i, which.max(zscores[i, ])[1]]
    trainDF$ModelZ[i] <- names(which.max(zscores[i, ]))
  }
  # # all the Leishmania genomes had the highest score with their prediceted model
  # table(trainDF$Model1 == trainDF$funclass)
  # table(trainDF$ModelZ == trainDF$funclass)
  # table(trainDF$ModelZ == trainDF$Model1)
  
  # combine those with model1 or modelZ different from the predicted model 
  
  diffrent <- (trainDF$Model1 != trainDF$funclass) #| (trainDF$ModelZ != trainDF$funclass)
  
  diffrentDF <- trainDF[diffrent, ]
  diffscores <- scoresdf[diffrent, ]
  diffzscores <- zscores[diffrent, ]
  
  # write two files : 1. the diffrentDF with first three columns
  # 2. the diffscores with the headers
  # 2. diffzscores with the headers
  # find the class function of those sequences that did not have same func between tse and ara
  undetDF<- data.frame(diffrentDF$geneid,diffrentDF$sequences,diffrentDF$funclass,diffrentDF$funcscore)
  names(undetDF) <- c("geneid","sequence","func","func_score")
  # add one column as potential shifts 
  # undetDF$pot_shifts_z <- ""
  # for (i in 1:nrow(diffzscores)) {
  #   undetDF$pot_shifts_z[i] <- paste(names(diffzscores)[ which(diffzscores[i, ] == max(diffzscores[i, ]))],collapse = " " , sep = " ")
  # }
  # 
  undetDF$pot_shifts <- ""
  undetDF$pot_shifts_s <- ""
  for (i in 1:nrow(diffzscores)) {
    undetDF$pot_shifts[i] <- paste(names(diffscores)[ which(diffscores[i, ] == max(diffscores[i, ]))],collapse = " " , sep = " ")
    undetDF$pot_shifts_s[i] <- paste(diffscores[i, which(diffscores[i, ] == max(diffscores[i, ]))],collapse = " " , sep = " ")
  }
  
  undetDF <- find.anticodon(undetDF)
  names <- data.frame("geneid","func","func_score","pot_shifts","pot_shifts_s","tseac","araac")
  names(names) <- c("geneid","func","func_score","pot_shifts","pot_shifts_s","tseac","araac")
  library(gdata)
  
  write.fwf(
    rbind(names,undetDF[,-(2)]),
    width = c(45, 7,10,10,12,7,7),
    file = paste(dirpath, "undetDF2.txt", sep = ""),
    colnames = FALSE
  )
  
  diffscores$geneid <- undetDF$geneid
  diffzscores$geneid <- undetDF$geneid
  
  n <- data.frame(t(names(diffscores)))
  names(n) <- names(diffscores)
  for (i in 1:ncol(n)) {
    n[,i] <- as.character(n[,i])
  }
  
  write.fwf(
    rbind(n,diffscores),
    width = c(rep(4,ncol(diffscores)-1),70),
    file = paste(dirpath, "undetDFscores.txt", sep = ""),
    colnames = FALSE
  )
  
  n <- data.frame(t(names(diffzscores)))
  names(n) <- names(diffzscores)
  for (i in 1:ncol(n)) {
    n[,i] <- as.character(n[,i])
  }
  write.fwf(
    rbind(n,diffzscores),
    width = c(rep(4,ncol(diffzscores)-1),70),
    file = paste(dirpath, "undetDFzscores.txt", sep = ""),
    colnames = FALSE
  )
  
  # find the anticodon for these genes
  
  trainDF
}

#__________________________________________________________________________________________
score.dist <- function(scoresdf,trainDF_Fun){

  par(mfrow=c(4,6))
  for (i in 1:ncol(scoresdf)) {
       samemodel <- names(scoresdf)[i] == trainDF_Fun
       currentmodel <- scoresdf[samemodel, i]
       hist(currentmodel, freq = FALSE)
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
  # also, calculate the z-score here
  # the final scores are saved in scoresdf
  # each column in scoresdf is for one model
  # those entries in that column that have the funcclass same as the name of the column
  # calculate the mean and std for that model and save it
  
  # resultfunc <- character(length = length(undetClass$sequences))
  # make a dataframe ncol = length(profile_list)), number of rows is nrow(undetClass)
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
      if (undetClass$funclass[x] == names(profile_list[j]))
        total = total - 1
      for (i in 1:ncol(currentDF)) {
        currChar <- substring(myseq, i, i)
        if (currChar == "-")
          currChar <- "gap"
        # if this sequence was annotated as a sequence of model names(profile_list[j]) then total = total - 1
        model_freq <- currentDF[currChar, i]
        if (undetClass$funclass[x] == names(profile_list[j]))
          model_freq <- currentDF[currChar, i] - 1
        scores[j] <- scores[j] + log(model_freq / total)
        if (types[j] != names(profile_list[j]))
          print("shit!")
      }
    }
    
    #undetClass$funclass[x] <- types[scores == max(scores)]
    scoresdf[x, ] <- scores
    
  }
  scoresdf
}

#___________________________________________________________________________________________
profile.cm <- function(trainDF, seqlength) {
  # used psedo counts to avoid zero for cells
  # read in the aligned sequences from file EditedCovea.fasta
  func_list <- split(trainDF, f = trainDF$funclass)
  # for each Df make a dataframe to keep the profile cm
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
