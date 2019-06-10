# This function will read the outlier.fasta file and subtract the from the intersection gene file.
# Make a model for each functional class and saves the results in a list of dataframes called profile_list
# At the end will also look at the distribution of the scores for each model and make sure that they are normally distributed
classifier <- function() {
  genedirpath <-
    "/home/fatemeh/Leishmania_2019/Leishmania_2019/Results/Integrated_Genes/TriTryp_EditedCovea_union.fasta"
  outlierpath <-
    "/home/fatemeh/Leishmania_2019/Leishmania_2019/Results/Integrated_Genes/TriTrypgenemodel_Intersection/TriTryp_EditedCovea/outlier/outliers.txt"
  outputdir <-
    "/home/fatemeh/Leishmania_2019/Leishmania_2019/Results/Integrated_Genes/"
  trainDF <- read.genefile(genedirpath) # 3543 genes
  undetDF <-
    trainDF[(trainDF$funclass == "undet"),]
  trainDF <-
    trainDF[(trainDF$funclass != "undet"),]# & (GeneDF$cluster == "T"), ]
  araonlyDF <- trainDF[nchar(trainDF$funclass) > 1, ]
  # gaps are removed from outier sequences by OD-seq package
  #outlierDF <- read.genefile(outlierpath) # 39 outier
  #trainDF <- remove.outliers(trainDF, outlierDF) # 3504 genes left
  trainDF_intersection <- trainDF[nchar(trainDF$funclass) == 1, ]
  pos_neg_list <- create.profile(trainDF_intersection)
  scoresdf_train <-
    calculate.score(pos_neg_list, undetDF, 1)
  araonlyDF$funclass <- gsub("ara", "", araonlyDF$funclass)
  #araonlyDF$funclass[araonlyDF$funclass=="M"] <- "X"
  scoresdf_araonly <- calculate.score(pos_neg_list, araonlyDF, 1)
  score.dist(
    scoresdf_train,
    trainDF_intersection$funclass,
    scoresdf_araonly,
    araonlyDF$funclass
  )
  
  # mark ara only genes whose score is  higher than the range of the scores of each model
  #find.true.aragenes(scoresdf_train, scoresdf_araonly)
  
  trainDF_Fun <- trainDF_intersection$funclass
  
  minscore_fore <- numeric(length = ncol(scoresdf_train))
  maxscore_back <- numeric(length = ncol(scoresdf_train))
  truegene <- logical(length = nrow(scoresdf_araonly))
  for (i in 1:ncol(scoresdf_train)) {
    foremodel <- names(scoresdf_train)[i] == trainDF_Fun
    foremodel_score <- scoresdf_train[foremodel, i]
    backmodel <- names(scoresdf_train)[i] != trainDF_Fun
    backmodel_score <- scoresdf_train[backmodel, i]
    # find the min of each model's foreground score
    minscore_fore[i] <- min(foremodel_score)
    maxscore_back[i] <- max(backmodel_score)
  }
  
  scoresdf_araonly_logical <-
    as.data.frame(matrix(
      data = FALSE,
      nrow = nrow(scoresdf_araonly),
      ncol = ncol(scoresdf_araonly)
    ))
  names(scoresdf_araonly_logical) <- names(scoresdf_araonly)
  for (i in 1:nrow(scoresdf_araonly)) {
    scoresdf_araonly_logical[i, ] <-
      (scoresdf_araonly[i, ] > minscore_fore) &
      (scoresdf_araonly[i, ] > maxscore_back)
    if (sum(scoresdf_araonly_logical[i, ] == TRUE) > 1)
    {
      print(sum(scoresdf_araonly_logical[i, ] == TRUE))
      print(names(scoresdf_araonly_logical)[scoresdf_araonly_logical[i, ] ==
                                              TRUE])
      print(araonlyDF$funclass[i])
      print(scoresdf_araonly[i, ])
      # make sure that it matched only one of the functions and that is the one the gene is marked with by aragorn
      
    }
    if ((sum(scoresdf_araonly_logical[i, ] == TRUE) == 1)) {
      if (araonlyDF$funclass[i] == names(scoresdf_araonly_logical)[scoresdf_araonly_logical[i, ] ==
                                                                   TRUE])
        truegene[i] <- TRUE
    }
  }
  
  maxfun <- character(length = nrow(scoresdf))
  maxarr <- numeric(length = nrow(scoresdf))
  for (i in 1:nrow(scoresdf)) {
    maxfun[i] <- names(scoresdf)[scoresdf[i, ] == max(scoresdf[i, ])]
    maxarr[i] <- max(scoresdf[i, ])
  }
  maxfunara <- paste("ara", maxfun, sep = "")
  maxarr[maxfunara == araonlyDF$funclass]
  araonlyDF[maxfunara == araonlyDF$funclass, ]
  
  outliersDF <- trainDF[maxfun != trainDF$funclass, ]
  trainDF <- remove.outliers(trainDF, outliersDF)
  table(maxfun == trainDF$funclass)
  types[scoresdf == maxarr]
  
  trainDF <- make.clusters(trainDF)
  # make the 22 models
  # trainDF <- trainDF[trainDF$cluster == "T", ]
  profile_list <- profile.cm(trainDF)
  # print(paste("we have", length(names(profile_list)), "functional classes:"))
  # print(names(profile_list))
  # use the trainDF as the undet class for now
  scoresdf <- find.seqScore(profile_list, trainDF)
  
  zscores <- find.Zscore(scoresdf, trainDF$funclass)
  zscores <- round(zscores)
  trainDF
}

# find.true.aragenes(scoresdf_train, trainDF_Fun, scoresdf_araonly, ara_Fun) {
#   minscore <- numeric(length = ncol(scoresdf_train))
#   truegene <- logical(length = ncol(scoresdf_train))
#   for (i in 1:ncol(scoresdf_train)) {
#     foremodel <- names(scoresdf_train)[i] == trainDF_Fun
#     foremodel_score <- scoresdf_train[foremodel, i]
#     # find the min of each model's foreground score
#     minscore[i] <- min(foremodel_score[, i])
#   }
#   for (i in 1:nrow(scoresdf_araonly)) {
#     scoresdf_araonly[i, ] > minscore
#   }
# }
#__________________________________________________________________________________________
calculate.score <- function(pos_neg_list, trainDF, araonly) {
  pos_profile_list <- pos_neg_list[[1]]
  neg_profile_list <- pos_neg_list[[2]]
  m <- matrix(ncol = length(pos_profile_list),
              nrow = nrow(trainDF),
              0)
  scoresdf <- as.data.frame(m)
  types <- names(pos_profile_list)
  names(scoresdf) <- types
  
  for (x in 1:length(trainDF$sequences)) {
    myseq <- trainDF$sequences[x]
    scores <- integer(length = length(types))
    for (j in 1:length(types)) {
      current_posDF <- pos_profile_list[[j]]
      current_negDF <- neg_profile_list[[j]]
      for (i in 1:ncol(current_posDF)) {
        currChar <- substring(myseq, i, i)
        if (currChar == "-")
          currChar <- "gap"
        if (araonly == 0)
          #substr(trainDF$funclass[x], 1, 3) != "ara")
        {
          if (trainDF$funclass[x] == names(pos_profile_list[j]))
            scores[j] <-
              scores[j] + log((current_posDF[currChar, i] - 1) / current_negDF[currChar, i])
          else
            scores[j] <-
              scores[j] + log(current_posDF[currChar, i] / (current_negDF[currChar, i] -
                                                              1))
        }
        else
          scores[j] <-
            scores[j] + log(current_posDF[currChar, i] / current_negDF[currChar, i])
      }
    }
    scoresdf[x, ] <- scores
  }
  scoresdf
}
#__________________________________________________________________________________________
make.clusters <- function(trainDF) {
  trainDF$cluster <- ""
  for (i in 1:nrow(trainDF)) {
    firstletter <- substr(trainDF$geneid[i], 1, 1)
    if (firstletter == "L")
      trainDF$cluster[i] <- "L"
    if (firstletter == "T")
      trainDF$cluster[i] <- "T"
  }
  trainDF
}
#__________________________________________________________________________________________
remove.outliers <- function(trainDF, outlierDF) {
  isoutlier <- trainDF$geneid %in% outlierDF$geneid
  trainDF <- trainDF[!isoutlier, ]
}
#__________________________________________________________________________________________
read.genefile <- function(genedirpath) {
  # this function read the gene file into a dataframe as teh training set with columns: "header" "sequences" "funclass" "geneid"
  # at the end it will remove sequences of undet model
  #genedirpath <-
  #  "/home/fatemeh/Leishmania_2019/Leishmania_2019/Results/Integrated_Genes/TriTryp_EditedCovea_Initiator.fasta"
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
  #GeneDF <-
  #  GeneDF[(GeneDF$funclass != "undet"),]# & (GeneDF$cluster == "T"), ]
  GeneDF
}
#__________________________________________________________________________________________
score.dist <- 
  function(scoresdf,
           trainDF_Fun ,
           scoresdf_araonly,
           ara_Fun) {
    par(mfrow = c(4, 6))
    for (i in 1:ncol(scoresdf)) {
      print(i)
      foremodel <- names(scoresdf)[i] == trainDF_Fun
      foremodel_score <- scoresdf[foremodel, i]
      #currentmodel_score <- scoresdf[, i]
      backmodel <- names(scoresdf)[i] != trainDF_Fun
      backmodel_score <- scoresdf[backmodel, i]
      
      foremodel_ara <- names(scoresdf)[i] == ara_Fun
      foremodel_score_ara <- scoresdf_araonly[foremodel_ara, i]
      #currentmodel_score <- scoresdf[, i]
      backmodel_ara <- names(scoresdf)[i] != ara_Fun
      backmodel_score_ara <- scoresdf_araonly[backmodel_ara, i]
      
      hist(
        foremodel_score,
        prob = TRUE,
        main = "",
        xlab = paste("score against model",
                     names(scoresdf)[i],
                     sep = " "),
        ylab = paste("freq (total=", length(foremodel_score), ")", sep = ""),
        col = "blue"
      )
      #lines(density(foremodel_score), col="blue", lwd=2)
      
      hist(backmodel_score,
           prob = TRUE,
           add = T,
           col = "red")
      #lines(density(backmodel_score), col="red", lwd=2)
      c1 = rgb(255,
               255,
               0,
               max = 255,
               alpha = 100,
               names = "yellow")
      c2 = rgb(0,
               100,
               0,
               max = 255,
               alpha = 80,
               names = "lt.yellow")
      if (length(foremodel_score_ara) == 0)
        next
      hist(
        foremodel_score_ara,
        prob = TRUE,
        add = T,
        col = c1
      )
      #lines(density(foremodel_score_ara), col="green", lwd=2)
      # if(length(backmodel_score_ara)==0)
      #   next
      # hist(
      #   backmodel_score_ara,
      #   prob=TRUE,
      #   add=T,col=c2
      # )
      #lines(density(backmodel_score_ara), col="black", lwd=2)
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
find.anticodon <- function(undetDF) {
  genefilepath <-
    "/home/fatemeh/Leishmania_2019/Leishmania_2019/Results/Integrated_Genes/integrated_tse_ara.txt"
  geneDF <-
    read.table(genefilepath, header = TRUE, colClasses = "character")
  undetDF$tseac <- ""
  undetDF$araac <- ""
  for (i in 1:nrow(undetDF)) {
    undetDF$tseac[i] <- geneDF[undetDF$geneid[i] == geneDF$geneid, ]$tseac
    undetDF$araac[i] <-
      geneDF[undetDF$geneid[i] == geneDF$geneid, ]$araac
  }
  undetDF
}
# ___________________________________________________________________________________________
find.Zscore <- function(scoresdf, trainDF_Fun) {
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
    
    coveaDF[coveaDF$coveageneid == undetClassScoreDF$geneid[i], ]$headers <-
      paste(">", coveaDF[coveaDF$coveageneid == undetClassScoreDF$geneid[i], ]$coveageneid, "_", FUNC, sep = "")
    geneDF[geneDF$geneid == undetClassScoreDF$geneid[i], ]$Func <-
      FUNC
    # make undet sequences to remove out of this loop
    if (FUNC == "undet")
      coveaDF <-
      coveaDF[!(coveaDF$coveageneid == undetClassScoreDF$geneid[i]), ]
    
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
    scoresdf[x, ] <- scores
    
  }
  scoresdf
}
#___________________________________________________________________________________________
create.profile <- function(trainDF) {
  seqlength <- nchar(trainDF$sequences[1])
  m <- matrix(nrow = 5, ncol = seqlength)
  prf <- data.frame(m)
  rownames(prf) <- c("A", "C", "G", "T", "gap")
  functions <- unique(trainDF$funclass)
  pos_profile_list <- rep(list(prf), length(functions))
  neg_profile_list <- rep(list(prf), length(functions))
  names(pos_profile_list) <- functions
  names(neg_profile_list) <- functions
  for (i in 1:length(functions)) {
    m <- matrix(nrow = 5, ncol = seqlength)
    pos_prf <- data.frame(m)
    rownames(pos_prf) <- c("A", "C", "G", "T", "gap")
    neg_prf <- data.frame(m)
    rownames(neg_prf) <- c("A", "C", "G", "T", "gap")
    for (j in 1:seqlength) {
      # pseudo counts of one for each cell
      # go over all the sequences and see how many of them are A,C,G,T in their jth position
      isA = substring(trainDF$sequences, j, j) == "A"
      poscounts <- sum(trainDF[isA, ]$funclass == functions[i])
      negcounts <- sum(trainDF[isA, ]$funclass != functions[i])
      pos_prf[1, j] <- (poscounts + 1)
      neg_prf[1, j] <- negcounts + 1
      isC = substring(trainDF$sequences, j, j) == "C"
      poscounts <- sum(trainDF[isC, ]$funclass == functions[i])
      negcounts <- sum(trainDF[isC, ]$funclass != functions[i])
      pos_prf[2, j] <- (poscounts + 1)
      neg_prf[2, j] <- negcounts + 1
      isG = substring(trainDF$sequences, j, j) == "G"
      poscounts <- sum(trainDF[isG, ]$funclass == functions[i])
      negcounts <- sum(trainDF[isG, ]$funclass != functions[i])
      pos_prf[3, j] <- (poscounts + 1)
      neg_prf[3, j] <- negcounts + 1
      isT = substring(trainDF$sequences, j, j) == "T"
      poscounts <- sum(trainDF[isT, ]$funclass == functions[i])
      negcounts <- sum(trainDF[isT, ]$funclass != functions[i])
      pos_prf[4, j] <- (poscounts + 1)
      neg_prf[4, j] <- negcounts + 1
      isGap = substring(trainDF$sequences, j, j) == "-"
      poscounts <- sum(trainDF[isGap, ]$funclass == functions[i])
      negcounts <- sum(trainDF[isGap, ]$funclass != functions[i])
      pos_prf[5, j] <- (poscounts + 1)
      neg_prf[5, j] <- negcounts + 1
    }
    pos_profile_list[[i]] <- pos_prf
    neg_profile_list[[i]] <- neg_prf
    names(pos_profile_list[i]) == functions[i]
  }
  pos_neg_list <- list(pos_profile_list, neg_profile_list)
  pos_neg_list
}
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
