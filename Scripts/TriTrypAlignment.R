TriTrypAlignment <- function() {
  library(gdata)
  library(readr)
  library(Hmisc)
  
  geneDF <- read.input()
  geneDF_Int_noVar <-
    remove.Vararm.Intersection(geneDF[geneDF$foundby == "both",])
  geneDF_Int_noVar_noIntron <-
    remove.Intron.Intersection(geneDF_Int_noVar)
  write.genefile(geneDF_Int_noVar_noIntron)
  # run covea
  system(
    "covea /home/fatemeh/Leishmania_2019/Leishmania_2019/Scripts/TRNA2-euk.cm /home/fatemeh/Leishmania_2019/Leishmania_2019/Results/Integrated_Genes/Integrated_Genes_NoVarIntron.fasta > /home/fatemeh/Leishmania_2019/Leishmania_2019/Results/Integrated_Genes/Trytryp_genes_NoVarIntron.covea"
  )
  # covea processing
  
  dirpath <-
    "/home/fatemeh/Leishmania_2019/Leishmania_2019/Results/Integrated_Genes/"
  covea_filename <- "Trytryp_genes_NoVarIntron.covea"
  
  readSeqsIntoDf(dirpath, covea_filename)
  
  seqDB <-
    read_csv(paste(dirpath, "coveaDF.txt", sep = ""),
             col_types = cols(.default = col_character()))
  SSDB <-
    read_csv(paste(dirpath, "coveaDF_SS.txt", sep = ""))
  
  fasta_filename <- "TriTryp_EditedCovea.fasta"
  CS_filename <- "TriTryp_structfile.txt"
  editAlignment(seqDB, SSDB, dirpath, fasta_filename, CS_filename)
  #map2sprinzle(dirpath,fasta_filename,CS_filename)
  # Use Script for clustering the 
  }

write.genefile <- function(geneDF_Int_noVar_noIntron) {
  # this function will write the genefile with removed variable arm and removed intron into a fasta file format called Integrated_Genes_NoVarIntron.fasta
  write.fwf(
    data.frame(
      paste(">", geneDF_Int_noVar_noIntron$GeneID, sep = ""),
      geneDF_Int_noVar_noIntron$GeneSeq
    ),
    file = paste(resultpath, "Integrated_Genes_NoVarIntron.fasta", sep = ""),
    sep = "\n",
    colnames = FALSE
  )
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
  
  geneseq <- character(length = length(genefile$geneid))
  geness <- character(length = length(genefile$geneid))
  geneDF_fasta <-
    data.frame(
      genefile$geneid,
      genefile$tsegeneseq,
      genefile$arageneseq,
      genefile$tsegeness,
      genefile$arageness,
      genefile$foundby
    )
  names(geneDF_fasta) <-
    c("GeneID",
      "GeneSeqTse",
      "GeneSeqAra",
      "SSTSE",
      "SSARA",
      "foundby")
  
  genefile$arageneseq <- as.character(genefile$arageneseq)
  genefile$tsegeneseq <- as.character(genefile$tsegeneseq)
  geneDF_fasta$GeneID <- as.character(geneDF_fasta$GeneID)
  geneDF_fasta$GeneSeqTse <- as.character(geneDF_fasta$GeneSeqTse)
  geneDF_fasta$GeneSeqAra <- as.character(geneDF_fasta$GeneSeqAra)
  geneDF_fasta$SSARA <- as.character(geneDF_fasta$SSARA)
  geneDF_fasta$SSTSE <- as.character(geneDF_fasta$SSTSE)
  geneDF_fasta
}

remove.Vararm.Intersection <- function(geneDF_both) {
  # This function will remove the variable arms from both secondary structure and sequence
  # For genes found by both genefinders or tseonly: we use the TSE sequences as reference which is always reported up to base 73.
  # We will count number of arms by counting number of occuances of "><" in function countarms. if we had 4 arms, we remove tha third arm starting ">" and ending with "<"
  
  for (i in 1:nrow(geneDF_both)) {
    # count number of "><"
    numarms <- countarms(geneDF_both[i,]$SSTSE)
    if (numarms == 4)
    {
      varcor <-
        removearm(geneDF_both[i,]$GeneSeqTse, geneDF_both[i,]$SSTSE)
      firstchunk <-
        substring(geneDF_both$GeneSeqTse[i], 1, varcor[1] - 1)
      SSTSE_firstchunk <-
        substring(geneDF_both$SSTSE[i], 1, varcor[1] - 1)
      
      secondchunk <-
        substring(geneDF_both$GeneSeqTse[i],
                  varcor[2] + 1,
                  nchar(geneDF_both$GeneSeqTse[i]))
      SSTSE_secondchunk <-
        substring(geneDF_both$SSTSE[i],
                  varcor[2] + 1,
                  nchar(geneDF_both$SSTSE[i]))
      
      geneDF_both$GeneSeqTse[i] <-
        paste(firstchunk, secondchunk, sep = "")
      geneDF_both$SSTSE[i] <-
        paste(SSTSE_firstchunk, SSTSE_secondchunk, sep = "")
    }
    
    if (numarms < 3)
    {
      print("We have genes with unusual structure which need to be removed!")
    }
  }
  geneDF_both
}

remove.Intron.Intersection <-  function(geneDF_Int_noVar) {
  # This function will remove the introns from both secondary structure and sequence
  # For genes found by both genefinders: we use the TSE sequences as reference which is always reported up to base 73.
  # Based on TSE manuscript, in genes found by tse nucleotides matching the "consensus" tRNA model used in Cove analysis appear in upper case,
  # while introns and other nucleotides in non-conserved positions are printed in lower-case letters. So, we will remove all the lower case letters
  # the result will be saved in a dataframe with only two columns "geneid" and "geneseq"
  noVar_noIntron_df <-
    data.frame(geneDF_Int_noVar$GeneID, geneDF_Int_noVar$GeneSeqTse)
  names(noVar_noIntron_df) <- c("GeneID", "GeneSeq")
  noVar_noIntron_df$GeneID <- as.character(noVar_noIntron_df$GeneID)
  noVar_noIntron_df$GeneSeq <-
    as.character(noVar_noIntron_df$GeneSeq)
  
  for (i in 1:nrow(geneDF_Int_noVar)) {
    noVar_noIntron_df$GeneSeq[i] <-
      gsub("[[:lower:]]", "", noVar_noIntron_df$GeneSeq[i])
  }
  
  noVar_noIntron_df
}

countarms <- function(geneSS) {
  geneSSarr <-
    substring(geneSS, seq(1, nchar(geneSS), 1), seq(1, nchar(geneSS), 1))
  counter = 0
  flag = 0
  vararmend = 0
  for (i in 1:length(geneSSarr)) {
    if (geneSSarr[i] == "<" & flag == 0)
    {
      counter <- counter + 1
      flag = 1
      
    }
    if (geneSSarr[i] == ">")
      flag = 0
  }
  
  counter
}

removearm <- function(geneSeq, geneSS) {
  flag = 0
  counter = 0
  forward = 1
  varbeg = 0
  varend = 0
  geneSSarr <-
    substring(geneSS, seq(1, nchar(geneSS), 1), seq(1, nchar(geneSS), 1))
  for (i in 1:length(geneSSarr)) {
    if (geneSSarr[i] == ">" & flag == 0)
    {
      counter <- counter + 1
      flag = 1
      if (counter == 3)
      {
        # this is the begining of the variable arm
        # count howmany > you have
        # go forward untill you see that many <
        forward = 0
        varbeg = i
        while (geneSSarr[i] != "<")
        {
          if (geneSSarr[i] == ">")
            forward = forward + 1
          i = i + 1
        }
        while (forward != 0) {
          if (geneSSarr[i] == "<")
            forward = forward - 1
          
          i = i + 1
        }
        varend = i - 1
      }
    }
    if (geneSSarr[i] == "<")
      flag = 0
  }
  # print(varend)
  # print(varbeg)
  varcor <- c(varbeg, varend)
  varcor
}

readSeqsIntoDf <- function(dirpath,covea_filename) {
  # this function will read the sequences from covea file along with their secondary structure into a dataframe
  # will write the result in coveaDF.txt and coveaDF_SS.txt
  coveafilepath <- paste(dirpath, covea_filename, sep = "")
  CS <-
    grep("#=CS +",
         readLines(coveafilepath),
         value = TRUE)
  CStemp <-
    unlist(strsplit(CS, split = "\\s+"))
  CS <- CStemp[seq(2, length(CStemp), 2)]
  CS <- paste(CS, collapse = '')
  CSarr <- substring(CS, seq(1, nchar(CS), 1), seq(1, nchar(CS), 1))
  
  # read the name of the sequences into variable "seqnames"
  SQs <- grep("#=SQ +",
              readLines(coveafilepath),
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
           readLines(coveafilepath),
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
              readLines(coveafilepath),
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
  
  write_csv(seqDB, path  = paste(dirpath, "coveaDF.txt", sep = ""))
  write_csv(SSDB, path  = paste(dirpath, "coveaDF_SS.txt", sep = ""))
  
}

writeCovea <- function(SSDB, seqDB, dirpath,cs,fasta_filename,CS_filename) {
  # this function will translates . to - and writes the output as 
  coveaseqs <- character(length = ncol(seqDB) - 1)
  coveass <- character(length = ncol(seqDB) - 1)
  
  seqDB <- seqDB[, -ncol(seqDB)]
  SSDB <- SSDB[, -ncol(SSDB)]
  for (i in 1:(length(coveaseqs))) {
    coveaseqs[i] <-
      paste((as.data.frame(seqDB[, i])[, 1]), collapse = '')
    coveass[i] <-
      paste((as.data.frame(SSDB[, i])[, 1]), collapse = '')
    
  }
  
  mynames <- names(seqDB)[!names(seqDB) %in% "CS"]
  
  # write.fwf(
  #   data.frame(mynames,
  #              coveaseqs,
  #              coveass),
  #   file = paste(dirpath, "TryTrypHomoC_EditedCovea.covea", sep = ""),
  #   sep = "\n",
  #   colnames = FALSE
  # )
  
  # remove "."s from sequences and write them in a fasta file to run with covea
  for (i in 1:length(coveaseqs)) {
    coveaseqs[i] <- translate(coveaseqs[i], "[.]", "-")
  }
  
  write.fwf(
    data.frame(paste(">", mynames, sep = ""),
               coveaseqs),
    file = paste(dirpath, fasta_filename, sep = ""),
    sep = "\n",
    colnames = FALSE
  )
  
  #cs <- paste(seqDB$CS, collapse = '')
  write.fwf(data.frame(cs),
            file = paste(dirpath, CS_filename, sep = ""),
            colnames = FALSE)
  
}

editAlignment <- function(seqDB, SSDB, dirpath,fasta_filename,CS_filename) {
  # this function will remove:
  # 1. sites that have more than 99% gap (removing rows that have more than 99% ".") 
  # 2. sequences with more than 8 gaps
  # 3. sequences with letter n or N in them
  # at the end it will call the function writeCovea to write 
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
        if (gapperc > 99)
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
  
  writeCovea(SSDB, seqDB, dirpath,cs,fasta_filename,CS_filename)
  
}

map2sprinzle <- function(dirpath,fasta_filename,CS_filename){
  # this function will read the final aigned sequences with their secondary structure
  # assignes positions to each base according to sprinzle
  seqdf <-
    read.table(paste(dirpath, fasta_filename, sep = ""))
  seq_arr <- as.character(seqdf$V1)
  headers <- seq_arr[seq(1, length(seq_arr), 2)]
  sequences <- seq_arr[seq(2, length(seq_arr), 2)]
  
  con <- file(description=paste(dirpath,CS_filename,sep = ""), open="r")
  cs <- linn <-readLines(con)
  csarr <- unlist(strsplit(cs,split = "\\s"))
  CS <- csarr[length(csarr)]
  CSARR <- substring(CS, seq(1, nchar(CS), 1), seq(1, nchar(CS), 1))
  sprinzlepos <- integer(length = length(CSARR))
  if(substr(CS,1,16)==">>>>>>>..>>>>...")
    for (i in 1:16) 
      sprinzlepos[i] <- i
  # find three most conserved positions from pos 17 to 21

  freqm <- matrix(data = 0,nrow = 4,ncol = 5)
  
  # rows: A,C,T,G
  for (i in 1:length(sequences)) {
    for (j in 1:5) {
      if(substr(sequences[i],j,j)=="A")
      freqm[1,j] <- freqm[1,j]+1
      if(substr(sequences[i],j,j)=="C")
      freqm[2,j] <- freqm[2,j]+1
      if(substr(sequences[i],j,j)=="T")
      freqm[3,j] <- freqm[3,j]+1
      if(substr(sequences[i],j,j)=="G")
      freqm[4,j] <- freqm[4,j]+1
    }
  }
#  [,1] [,2] [,3] [,4] [,5]
#  [1,]   77   65  306  732  745
#  [2,]   90 1796 1704 1455 1088
#  [3,]  392  440  701  649  668
#  [4,] 3017 1275  865  739 1075
  sprinzlepos[17] <- 18
  sprinzlepos[18] <- 19
  sprinzlepos[19] <- 20
  sprinzlepos[20] <- 20
  sprinzlepos[21] <- 21
  if(substr(CS,22,25)=="<<<<")
    for (i in 22:25) 
      sprinzlepos[i] <- i
  sprinzlepos[26] <- 26 # not paired according to sprinzle!
  
  if(substr(CS,27,43)==">>>>>.......<<<<<")
    for (i in 27:43) 
      sprinzlepos[i] <- i
  sprinzlepos[44] <- 44 # not paired according to sprinzle!
  
  # 47 is not present
  sprinzlepos[45] <- 45
  sprinzlepos[46] <- 46
  sprinzlepos[47] <- 48
  
  if(substr(CS,48,64)==">>>>>.......<<<<<")
    for (i in 48:64) 
      sprinzlepos[i] <- i + 1
  if(substr(CS,64,72)=="<<<<<<<<.")
    for (i in 65:72) 
      sprinzlepos[i] <- i + 1
  
}