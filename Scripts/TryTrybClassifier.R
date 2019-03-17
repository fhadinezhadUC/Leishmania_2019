TTclassifier <- function() {
  # Make two clusters of tRNA genes: Leishmania and Trypanosoma
  genefilepath <-
    "/home/fatemeh/Leishmania_2019/Leishmania_2019/Results/Integrated_Genes/integrated_tse_ara.txt"
  geneDF <-
    read.table(genefilepath, header = TRUE, colClasses = "character")
  # 4381 genes 
  # use the tRNAs with same tsefunc and arafunc and not ""
  trainset <- geneDF[(geneDF$tsefunc == geneDF$arafunc) &
           (geneDF$tsefunc != ""), ]
  # 3560 genes in trainset
  profile_list <- profile.cm(trainset)
  
}

profile.cm <- function(trainDF) {
  func_list <- split(trainDF, f = trainDF$tsefunc)
  # for each Df make a dataframe to keep the profile cm
  seqlength <- nchar(sequences[1])
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
      Afunc = 0
      Cfunc = 0
      Tfunc = 0
      Gfunc = 0
      gapfunc = 0
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
      prf[1, j] <- Afunc / nrow(currentfuncdf) # A
      prf[2, j] <- Cfunc / nrow(currentfuncdf) # C
      prf[3, j] <- Gfunc / nrow(currentfuncdf) # G
      prf[4, j] <- Tfunc / nrow(currentfuncdf) # T
      prf[5, j] <- gapfunc / nrow(currentfuncdf) # gap
    }
    # devide each column in prf by the number of sequences in the current funcclass
    # make sure sum of each column is asame
    profile_list[[i]] <- prf
    names(profile_list[i]) == names(func_list[i])
  }
  
  
  profile_list
}
