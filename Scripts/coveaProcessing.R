# This script will Read the aligned sequences in seqDB dataframe and secondary structures in SSDB and saves them in files:
# coveaDF.txt and coveaDF_SS.txt. Also, save the CS line in struct_file.txt
# Removing sites that have more than 99% gap (removing rows in seqDB dataframe that have more than 99% ".")
# Removing genes that have more than 8 gaps in their aligned sequence
# remove sequence which has N/n in their sequence
# Save the Edited covea file in EditedCovea.covea file. Also as a fasta file (with "."s replaced with "-") in EditedCovea.fasta to be converted to clustal format for tsfm
# We ended up with ? aligned genes 
coveaProcessing <- function() {
  dirpath <-
    "/home/fatemeh/Leishmania_2019/Leishmania_2019/Results/Integrated_Genes/"
  readSeqsIntoDf(dirpath)
  library(readr)
  seqDB <-
    read_csv(paste(dirpath, "coveaDF.txt", sep = ""),col_types = cols(.default = col_character()))
  SSDB <-
    read_csv(paste(dirpath, "coveaDF_SS.txt", sep = ""))
  editAlignment(seqDB, SSDB, dirpath)

}
#####################################################################################
# # of columns in the dataframe is: # of sequences + 1 CS line
# name of the columns are "CS","seqname1", "seqname2", ...
readSeqsIntoDf <- function(dirpath) {
  CS <-
    grep("#=CS +",
         readLines(paste(
           dirpath, "TryTrypHomoC.covea", sep = ""
         )),
         value = TRUE)
  CStemp <-
    unlist(strsplit(CS, split = "\\s+"))
  CS <- CStemp[seq(2, length(CStemp), 2)]
  CS <- paste(CS, collapse = '')
  CSarr <- substring(CS, seq(1, nchar(CS), 1), seq(1, nchar(CS), 1))
  
  # read the name of the sequences into variable "seqnames"
  SQs <- grep("#=SQ +",
              readLines(paste(
                dirpath, "TryTrypHomoC.covea", sep = ""
              )),
              value = TRUE)
  seqnames <- character(length = length(SQs))
  for (i in 1:length(SQs)) {
    seqnames[i] <- unlist(strsplit(SQs[i], split = " "))[2]
  }
  
  # define the main data frame as "seqDB" and assign names to it
  m <- matrix(ncol = (length(seqnames) + 1), nrow = length(CSarr))
  seqDB  <- as.data.frame(m)
  names(seqDB) <- c(seqnames, "CS")
  
  for (i in 1:length(seqnames)) {
    pat <- paste("^", seqnames[i], " +", sep = "")
    myseq <-
      grep(pattern = pat,
           readLines(
             paste(dirpath, "TryTrypHomoC.covea", sep = "")
           ),
           value = TRUE)
    temp <-
      unlist(strsplit(myseq, split = "\\s+"))
    myseq <- temp[seq(2, length(temp), 2)]
    myseq <- paste(myseq, collapse = '')
    myseqarr <-
      substring(myseq, seq(1, nchar(myseq), 1), seq(1, nchar(myseq), 1))
    seqDB[, i] <- myseqarr
  }
  seqDB$CS <- CSarr
  
  #_____________________ reading the #=SS lines into a nother data frame as SSDB ____________
  # define the main data frame as "seqDB" and assign names to it
  
  SSDB = seqDB
  
  SSs <- grep("#=SS +",
              readLines(paste(
                dirpath, "TryTrypHomoC.covea", sep = ""
              )),
              value = TRUE)
  
  # number of sections is CStemp/2 = 3
  numsec <- length(CStemp) / 2
  end <- length(SSs) / numsec
  for (i in 1:end) {
    if (numsec == 3)
      myseq <-
        paste(unlist(strsplit(SSs[i], split = "\\s+"))[2],
              unlist(strsplit(SSs[end + i], split = "\\s+"))[2],
              unlist(strsplit(SSs[(2 * end) + i], split = "\\s+"))[2],
              sep = "")
    if (numsec == 2)
      myseq <-
        paste(unlist(strsplit(SSs[i], split = "\\s+"))[2], unlist(strsplit(SSs[end +
                                                                                 i], split = "\\s+"))[2], sep = "")
    
    myseqarr <-
      substring(myseq, seq(1, nchar(myseq), 1), seq(1, nchar(myseq), 1))
    SSDB[, i] <- myseqarr
  }
  
  library(readr)
  write_csv(seqDB, path  = paste(dirpath, "coveaDF.txt", sep = ""))
  write_csv(SSDB, path  = paste(dirpath, "coveaDF_SS.txt", sep = ""))
  
}
#####################################################################################
writeCovea <- function(SSDB, seqDB, dirpath,cs) {
  coveaseqs <- character(length = ncol(seqDB) - 1)
  coveass <- character(length = ncol(seqDB) - 1)
  
  seqDB <- seqDB[,-ncol(seqDB)]
  SSDB <- SSDB[,-ncol(SSDB)]
  for (i in 1:(length(coveaseqs))) {
    # coveaseqs[i] <- paste(seqDB[, i], collapse = '')
    # coveass[i] <- paste(SSDB[, i], collapse = '')
    print(i)
    coveaseqs[i] <-
      paste((as.data.frame(seqDB[, i])[, 1]), collapse = '')
    coveass[i] <-
      paste((as.data.frame(SSDB[, i])[, 1]), collapse = '')
    
  }
  
  mynames <- names(seqDB)[!names(seqDB) %in% "CS"]
  
  library(gdata)
  write.fwf(
    data.frame(mynames,
               coveaseqs,
               coveass),
    file = paste(dirpath, "TryTrypHomoC_EditedCovea.covea", sep = ""),
    sep = "\n",
    colnames = FALSE
  )
  
  # remove "."s from sequences and write them in a fasta file to run with covea
  library(Hmisc)
  for (i in 1:length(coveaseqs)) {
    coveaseqs[i] <- translate(coveaseqs[i], "[.]", "-")
    #coveaseqs[i] <- translate(coveaseqs[i], "[n]", "-")
    #coveaseqs[i] <- translate(coveaseqs[i], "[N]", "-")
  }
  write.fwf(
    data.frame(paste(">", mynames, sep = ""),
               coveaseqs),
    file = paste(dirpath, "TryTrypHomoC_EditedCovea.fasta", sep = ""),
    sep = "\n",
    colnames = FALSE
  )
  
  #cs <- paste(seqDB$CS, collapse = '')
  write.fwf(data.frame(cs),
    file = paste(dirpath, "TryTrypHomoC_struct_file.txt", sep = ""),
    colnames = FALSE
  )
  
}
#####################################################################################
editAlignment <- function(seqDB, SSDB, dirpath) {
  # removing sites that have more than 99% gap (removing rows that have more than 99% ".")
  # seqDB = mydb
  delpos = " "
  for (i in 1:nrow(seqDB)) {
    temp <-
      data.frame(table(
        seqDB[i,] == "." |
          seqDB[i,] == "t" |
          seqDB[i,] == "g" | seqDB[i,] == "c" | seqDB[i,] == "a"
      ))
    if (length(temp[temp$Var1 == "FALSE", 2]) != 0)
    {
      if (temp[temp$Var1 == "FALSE", 2] != ncol(seqDB))
      {
        gapperc <-
          (temp[temp$Var1 == "TRUE", 2] / (temp[temp$Var1 == "TRUE", 2] + temp[temp$Var1 ==
                                                                                 "FALSE", 2])) * 100
        if (gapperc > 97)
        {
          delpos = c(delpos, i)
        }
        
      }
    }
    else{
      if (temp[temp$Var1 == "TRUE", 2] == ncol(seqDB))
        delpos = c(delpos, i)
    }
  }
  delpos <- delpos[2:length(delpos)]
  delpos = as.integer(delpos)
  seqDB <- seqDB[-delpos,]
  SSDB <- SSDB[-delpos,]
  
  cs <- paste(seqDB$CS, collapse = '')
  
  # remove sequences with more than 8 gaps!
  delpos = " "
  flag = 0
  for (i in 1:ncol(seqDB)) {
    temp <- data.frame(table(seqDB[, i] == "."))
    if (length(temp[temp$Var1 == "TRUE", 2]) != 0)
      if (temp[temp$Var1 == "TRUE", 2] > 8)
      {
        if (names(seqDB[, i]) != "CS")
          delpos = c(delpos, i)
      }
    for (x in 1:length(unlist(seqDB[, i]))) {
      if (unlist(seqDB[, i])[x] == "n" | unlist(seqDB[, i])[x] == "N")
        flag = 1
    }
    if (flag == 1)
    {
      print("shit")
      delpos = c(delpos, i)
      flag = 0
    }
  }
  if (length(delpos) > 1)
  {
    delpos <- delpos[2:length(delpos)]
    delpos = as.integer(delpos)
    seqDB <- seqDB[,-delpos]
    SSDB <- SSDB[,-delpos]
  }
  
  # remove sites that are lowercase or . in 99% of sequences
  
  # remove sequences with gap in their anticodon
  
  writeCovea(SSDB, seqDB, dirpath,cs)
  
}
#####################################################################################
# partition the fasta file based on their functional classes
# partition the EditedCovea.fasta files based on the sourceOrg
