# This script will take the integrated_Tse_Ara gene file (output of Integrate_Tse_Ara script),
# Extract the genefiles found by both TSE and ARA with same identity 3484
# We dismissed 13 genes which were marked as either pseudo or truncated by TSE.
# There were 3471 genes left to be processed.
# In TSE  On the sequence line, nucleotides matching the "consensus" tRNA model used in Cove analysis appear in upper case,
# while introns and other nucleotides in non-conserved positions are printed in lower-case letters.
# The script will uses TSE sequences as reference since they are all reported up to position 73.
# Make a new dataframe with columns: "GeneID", "GeneSeq", "SSTSE", "SSARA"
# will add string "_undet" at the end of the genes with mismatched identity by ARA and TSE
# remove the introns and non-conserved positions (lowercase letters) and variable arms
# output the result in file Integrated_Genes_NoVarIntron.fasta
# the output will be used for covea to align:
# covea TRNA2-euk.cm Integrated_Genes_NoVarIntron.fasta > Trytryp_genes_NoVarIntron.covea
remove.Intron.Vararm <- function() {
  library(gdata)
  resultpath <-
    "/home/fatemeh/Leishmania_2019/Leishmania_2019/Results/Integrated_Genes/"
  genefilepath <-
    "/home/fatemeh/Leishmania_2019/Leishmania_2019/Results/Integrated_Genes/integrated_tse_ara.txt"
  geneDF <-
    read.table(genefilepath, header = TRUE, colClasses = "character")
  # we use genes that are found by both tse and ara and their identity matches.
  both <- geneDF$foundby == "both"
  genefile <- geneDF[both, ]
  isvagueIdentity <- (genefile$tsefunc != genefile$arafunc)
  genefile$geneid[!isvagueIdentity] <-
    paste(genefile$geneid[!isvagueIdentity], "_", genefile$tsefunc[!isvagueIdentity], sep = "")
  genefile$geneid[isvagueIdentity] <-
    paste(genefile$geneid[isvagueIdentity], "_undet", sep = "")
  # genefile <- genefile[!isvagueIdentity, ]
  nonpseudo <- genefile$tsenote == "notfound"
  genefile <- genefile[nonpseudo,]
  
  genefile_Nointron <-
    data.frame(
      genefile$geneid,
      genefile$tsegeneseq,
      genefile$tsegeness,
      genefile$arageness
    )
  names(genefile_Nointron) <-
    c("GeneID", "GeneSeq", "SSTSE", "SSARA")
  
  genefile$arageneseq <- as.character(genefile$arageneseq)
  genefile$tsegeneseq <- as.character(genefile$tsegeneseq)
  genefile_Nointron$GeneID <- as.character(genefile_Nointron$GeneID)
  genefile_Nointron$GeneSeq <-
    as.character(genefile_Nointron$GeneSeq)
  genefile_Nointron$SSARA <- as.character(genefile_Nointron$SSARA)
  genefile_Nointron$SSTSE <- as.character(genefile_Nointron$SSTSE)
  
  # remove the variable arms first and introns second.
  genefile_Nointron <- removeVariableArm(genefile_Nointron)
  
  genefile_Nointron$SeqNointron <- genefile_Nointron$GeneSeq
  for (i in 1:nrow(genefile)) {
    genefile_Nointron$SeqNointron[i] <-
      gsub("[[:lower:]]", "", genefile_Nointron$GeneSeq[i])
  }
  
  write.fwf(
    data.frame(
      paste(">", genefile_Nointron$GeneID, sep = ""),
      genefile_Nointron$SeqNointron
    ),
    file = paste(resultpath, "Integrated_Genes_NoVarIntron.fasta", sep = ""),
    sep = "\n",
    colnames = FALSE
  )
}

removeVariableArm <- function(genefile_Nointron) {
  for (i in 1:nrow(genefile_Nointron)) {
    # count number of "><"
    numarms <- countarms(genefile_Nointron[i,]$SSTSE)
    if (numarms == 4)
    {
      varcor <-
        removearm(genefile_Nointron[i,]$GeneSeq, genefile_Nointron[i,]$SSTSE)
      firstchunk <-
        substring(genefile_Nointron$GeneSeq[i], 1, varcor[1] - 1)
      SSTSE_firstchunk <-
        substring(genefile_Nointron$SSTSE[i], 1, varcor[1] - 1)
      
      secondchunk <-
        substring(genefile_Nointron$GeneSeq[i],
                  varcor[2] + 1,
                  nchar(genefile_Nointron$GeneSeq[i]))
      SSTSE_secondchunk <-
        substring(genefile_Nointron$SSTSE[i],
                  varcor[2] + 1,
                  nchar(genefile_Nointron$SSTSE[i]))
      
      genefile_Nointron$GeneSeq[i] <-
        paste(firstchunk, secondchunk, sep = "")
      genefile_Nointron$SSTSE[i] <-
        paste(SSTSE_firstchunk, SSTSE_secondchunk, sep = "")
    }
    if (numarms < 3)
    {
      print("We have genes with unusual structure!")
    }
  }
  genefile_Nointron
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
  
  print(varend)
  print(varbeg)
  varcor <- c(varbeg, varend)
  varcor
}