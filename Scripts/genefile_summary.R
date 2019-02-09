genefile.summary <- function(){
  # read integrated_tse_ara file into a dataframe
  outputpath <- "/home/fatemeh/Leishmania_2019/Leishmania_2019/Results/Integrated_Genes/"
  genefilepath <- "/home/fatemeh/Leishmania_2019/Leishmania_2019/Results/Integrated_Genes/integrated_tse_ara.txt"
  geneDF <- read.table(genefilepath, header = TRUE,colClasses = "character")
  # make a summary table for the 4 sets of genes found by ara, tse, their intersection, and union
  summerytable <- summery_table(geneDF,outputpath)
  
  # make a plot to show number of genes annotated for each genome,
  # make a plot to show which tRNA functional classes were not annotated by both TSE and ARA for each genome.
  vidualization(geneDF,outputpath)
  
  # make a plot to show TSE VS Ara gene displacement
  endDisplacement(geneDF,outputpath)
}
vidualization <- function(integrated_tse_ara,outputpath){
  
  genometRNAdf <-
    data.frame(
      integrated_tse_ara$sourceOrg,
      integrated_tse_ara$arafunc,
      integrated_tse_ara$tsefunc,
      integrated_tse_ara$foundby
    )
  names(genometRNAdf) <- c("sourceOrg", "arafunc", "tsefunc","foundby")
  genome_list = split(genometRNAdf, f = genometRNAdf$sourceOrg)
  
  tRNAcountdf <- data.frame(names(genome_list))
  tRNAcountdf$aracounts <- 0
  tRNAcountdf$tsecounts <- 0
  tRNAcountdf$unioncounts <- 0
  tRNAcountdf$interscounts <- 0
  names(tRNAcountdf) <- c("genome", "aracounts","tsecounts","unioncounts","interscounts")
  
  for (i in 1:length(genome_list)) {
    tRNAcountdf$genome[i] <- names(genome_list[i])
    currentgenedf <- genome_list[[i]]
    # for each genome df make the union, inters, tse and ara df
    both <- currentgenedf$foundby == "both"
    tseonly <- currentgenedf$foundby == "tse"
    araonly <- currentgenedf$foundby == "ara"
    
    istse <- (tseonly | both)
    isara <- (araonly | both)
    isintersect_tse_ara <- both
    
    tse <- currentgenedf[istse, ]
    ara <- currentgenedf[isara, ]
    
    # for intersection if tse and aragorn do not match we do not consider them in our counting
    intersect_tse_ara <- currentgenedf[isintersect_tse_ara, ]
    intersect_tse_ara <- intersect_tse_ara[as.character(intersect_tse_ara$tsefunc)==as.character(intersect_tse_ara$arafunc),]
    union_tse_ara <- currentgenedf
    
    tRNAcountdf$aracounts[i] <- nrow(ara)
    tRNAcountdf$tsecounts[i] <- nrow(tse)
    tRNAcountdf$unioncounts[i] <- nrow(union_tse_ara)
    tRNAcountdf$interscounts[i] <- nrow(intersect_tse_ara)
  }
  # keep this order for after filtering vidualization
  tRNAcountdf <- tRNAcountdf[order(tRNAcountdf$interscounts), ]
  tRNAcountdf$genome <-
    factor(tRNAcountdf$genome, levels = tRNAcountdf$genome)
  # 
  # mymatrix <- rbind(tRNAcountdf$aracounts,tRNAcountdf$tsecounts)
  # colnames(mymatrix) <- tRNAcountdf$genome
  # barplot(mymatrix,beside = TRUE,col = c("blue","red"),las=2,cex.names=0.6)
  # 
  library(ggplot2)
  p <- ggplot(data = tRNAcountdf, aes(x = genome, y = interscounts)) +
    geom_bar(stat = "identity", fill = "steelblue") + geom_text(
      aes(label = interscounts),
      vjust = 1.2,
      color = "white",
      position = position_dodge(0.9),
      size = 3.5
    ) + theme(axis.text.x = element_text(angle = 90, hjust = 1),
              plot.title = element_text(hjust = 0.5)) +
    labs(title =
           "Number of tRNA genes (intersection of genes found by ARA and TSE) in TryTryp genomes",
         x =
           "Genome", y = "Number of tRNA genes")
  p
  ggsave(
    paste(outputpath, "intersecttRNAcounts.png", sep = ""),
    width = 14,
    height = 7
  )
  
  # second plot _______________________________________________________________________________
  tRNAcountdf$genome <-
    as.character(tRNAcountdf$genome, levels = tRNAcountdf$genome)
  
  func_classes <-
    c(
      "A",
      "R",
      "N" ,
      "D",
      "C" ,
      "Q",
      "E",
      "G",
      "H",
      "I",
      "L",
      "K",
      "M" ,
      "X",
      "F",
      "P",
      "S",
      "T",
      "W",
      "Y" ,
      "V",
      "Z"
    )
  tRNAcountdf$percent_21 <- 0
  tRNAcountdf$missing <- ''
  for (i in 1:length(genome_list)) {
    # remove function Z for now
    currentgenedf <- genome_list[[i]]
    both <- currentgenedf$foundby == "both"
    intersect_tse_ara <- currentgenedf[both,]
    intersect_tse_ara <- intersect_tse_ara[as.character(intersect_tse_ara$tsefunc)==as.character(intersect_tse_ara$arafunc),]
    
    countsdf <-
      as.data.frame(table(as.character(intersect_tse_ara$arafunc)))
    tRNAcountdf[tRNAcountdf$genome == names(genome_list[i]),]$missing <-
      paste(setdiff(func_classes, as.character(countsdf$Var1)), collapse = " ")
    diffset <- setdiff(func_classes, as.character(countsdf$Var1))
    fcount <- 22 - length(diffset)
    tRNAcountdf[tRNAcountdf$genome == names(genome_list[i]),]$percent_21 <-
      (fcount / 22) * 100
    # make a vector of missing functional classes
  }
  tRNAcountdf$genome <-
    factor(tRNAcountdf$genome, levels = tRNAcountdf$genome)
  p <- ggplot(data = tRNAcountdf, aes(x = genome, y = percent_21)) +
    geom_bar(stat = "identity", fill = "#56B4E9")  + theme(axis.text.x = element_text(angle = 90, hjust = 1),
                                                           plot.title = element_text(hjust = 0.5)) +
    labs(title =
           "Percentage of 21 tRNA functional classes covered by each genome considering genes found by ARA only",
         x =
           "Genome", y = "Percentage of tRNA functional classes(total=21)") + geom_text(
             aes(label =
                   missing),
             vjust = 0.5,
             angle = 90,
             hjust = -0.1
           )
  p
  ggsave(paste(outputpath, "ara_funcPerc.png", sep = ""),
         width = 14,
         height = 7)
  
  #__________________ third plot 
  inter_genes <- integrated_tse_ara[integrated_tse_ara$foundby == "both",]
  inter_diff_type <- inter_genes[inter_genes$arafunc!=inter_genes$tsefunc,]
  table(paste(inter_diff_type$arafunc,inter_diff_type$tsefunc,sep = "|"))
  
}
summery_table <- function(integrated_tse_ara,outputpath) {
  #______________________________________________________________________________________________________________________________________________________
  # making a dataframe table with the four rows for each geneset of:
  # 1. tse2
  # 2. ARA
  # 3. Union(tse2,ARA)
  # 4. Intersection(tse2,ARA)
  # with 35 columns for  T (#trnas), N (#nucleotides), N/T, length range, %G, %C, %T, %A, %intron-containing, and the 23 class frequencies including
  # IUPAC aa codes for the 20 elongators(A, C, E, ... , Y), X for initiators, Z for selenocysteine, $ for pseudogenes, ? for sup, # for stop and O for pyl
  
  both <- integrated_tse_ara$foundby == "both"
  tseonly <- integrated_tse_ara$foundby == "tse"
  araonly <- integrated_tse_ara$foundby == "ara"
  
  istse <- (tseonly | both)
  isara <- (araonly | both)
  isintersect_tse_ara <- both
  
  tse <- integrated_tse_ara[istse, ]
  ara <- integrated_tse_ara[isara, ]
  
  # for intersection if tse and aragorn do not match we do not consider them in our counting
  intersect_tse_ara <- integrated_tse_ara[isintersect_tse_ara, ]
  intersect_tse_ara <- intersect_tse_ara[intersect_tse_ara$tsefunc==intersect_tse_ara$arafunc,]
  union_tse_ara <- integrated_tse_ara
  
  # making a 4 * 36 size dataframe
  m <- matrix(nrow = 4,
              ncol = 36,
              data = 0)
  summerytable <- as.data.frame(m)
  names(summerytable) <-
    c(
      "GeneSet",
      "#tRNA",
      "#nucleotides",
      "N/T",
      "lengthrange",
      "%G",
      "%C",
      "%T",
      "%A",
      "%intronContaining",
      "A",
      "C",
      "D",
      "E",
      "F",
      "G",
      "H",
      "I",
      "K",
      "L",
      "M",
      "N",
      "P",
      "Q",
      "R",
      "S",
      "T",
      "V",
      "W",
      "Y",
      "X",
      "Z",
      "$",
      "?",
      "#",
      "O"
    )
  
  summerytable[, 1] <- c("TSE2", "ARA", "UNION", "INTERSECTION")
  summerytable[, 2] <-
    c(nrow(tse),
      nrow(ara),
      nrow(union_tse_ara),
      nrow(intersect_tse_ara))
  
  # sum of the gene length for union
  unionnucleotides <-
    sum(nchar(tse$tsegeneseq)) + sum(nchar(ara$arageneseq)) - sum(nchar(intersect_tse_ara$tsegeneseq))
  summerytable[, 3] <-
    c(sum(nchar(tse$tsegeneseq)) , sum(nchar(ara$arageneseq)) , unionnucleotides , sum(nchar(intersect_tse_ara$tsegeneseq)))
  summerytable[, 4] <-
    c(
      summerytable[1, 3] / summerytable[1, 2],
      summerytable[2, 3] / summerytable[2, 2],
      summerytable[3, 3] / summerytable[3, 2],
      summerytable[4, 3] / summerytable[4, 2]
    )
  tserange <-
    paste(min(nchar(tse$tsegeneseq)) , max(nchar(tse$tsegeneseq)) , sep = "-")
  ararange <-
    paste(min(nchar(ara$arageneseq)) , max(nchar(ara$arageneseq)) , sep = "-")
  unionrange <-
    paste(min(min(nchar(tse$tsegeneseq)), min(nchar(ara$arageneseq))) ,  max(max(nchar(tse$tsegeneseq)), max(nchar(ara$arageneseq))) , sep = "-")
  intersectionrange <-
    paste(min(nchar(intersect_tse_ara$tsegeneseq)) , max(nchar(intersect_tse_ara$tsegeneseq)) , sep =  "-")
  summerytable[, 5] <-
    c(tserange, ararange, unionrange, intersectionrange)
  
  #_____________________ A T C G percentage ______________________________
  
  tse$tsegeneseq <- tolower(tse$tsegeneseq)
  count <-
    sapply(c("a", "g", "c", "t"), function(nuc)
      str_count(tse$tsegeneseq, fixed(nuc)))
  countdf <- as.data.frame(count)
  A <- sum(countdf$a)
  G <- sum(countdf$g)
  C <- sum(countdf$c)
  TT <- sum(countdf$t)
  total <- sum(A, G, C, TT)
  summerytable$`%A`[1] <- A * 100 / total
  summerytable$`%G`[1] <- G * 100 / total
  summerytable$`%C`[1] <- C * 100 / total
  summerytable$`%T`[1] <- TT * 100 / total
  
  ara$arageneseq <- tolower(ara$arageneseq)
  count <-
    sapply(c("a", "g", "c", "t"), function(nuc)
      str_count(as.character(ara$arageneseq), fixed(nuc)))
  countdf <- as.data.frame(count)
  A <- sum(countdf$a)
  G <- sum(countdf$g)
  C <- sum(countdf$c)
  TT <- sum(countdf$t)
  total <- sum(A, G, C, TT)
  summerytable$`%A`[2] <- A * 100 / total
  summerytable$`%G`[2] <- G * 100 / total
  summerytable$`%C`[2] <- C * 100 / total
  summerytable$`%T`[2] <- TT * 100 / total
  
  intersect_tse_ara$tsegeneseq <-
    tolower(intersect_tse_ara$tsegeneseq)
  count <-
    sapply(c("a", "g", "c", "t"), function(nuc)
      str_count(as.character(intersect_tse_ara$tsegeneseq), fixed(nuc)))
  countdf <- as.data.frame(count)
  A <- sum(countdf$a)
  G <- sum(countdf$g)
  C <- sum(countdf$c)
  TT <- sum(countdf$t)
  total <- sum(A, G, C, TT)
  summerytable$`%A`[4] <- A * 100 / total
  summerytable$`%G`[4] <- G * 100 / total
  summerytable$`%C`[4] <- C * 100 / total
  summerytable$`%T`[4] <- TT * 100 / total
  
  seq <- integrated_tse_ara$tsegeneseq
  
  for (i in 1:length(seq)) {
    if (integrated_tse_ara$foundby[i] == "ara")
      seq[i] <- integrated_tse_ara$arageneseq[i]
  }
  
  seq <- tolower(seq)
  count <-
    sapply(c("a", "g", "c", "t"), function(nuc)
      str_count(as.character(seq), fixed(nuc)))
  countdf <- as.data.frame(count)
  A <- sum(countdf$a)
  G <- sum(countdf$g)
  C <- sum(countdf$c)
  TT <- sum(countdf$t)
  total <- sum(A, G, C, TT)
  summerytable$`%A`[3] <- A * 100 / total
  summerytable$`%G`[3] <- G * 100 / total
  summerytable$`%C`[3] <- C * 100 / total
  summerytable$`%T`[3] <- TT * 100 / total
  
  #_____________________ intron containing% ________________________________
  
  # araintronbegin is the location of the intron in the tRNA gene not genome!
  arahasintron <- ara$araintronbegin != "nointron"
  nrow(ara[arahasintron, ])
  
  tsehasintron <-
    (tse$tseintronbegin != "0" & tse$tseintronbegin != "notfound")
  nrow(tse[tsehasintron, ])
  
  intersecthasintron <-
    ((intersect_tse_ara$tseintronbegin == "0") &
       (intersect_tse_ara$araintronbegin == "nointron")
    )
  nrow(intersect_tse_ara[!intersecthasintron, ])
  
  
  unionhasintron <-
    (((union_tse_ara$tseintronbegin == "0") |
        (union_tse_ara$tseintronbegin == "notfound")
    ) &
      (union_tse_ara$araintronbegin == "nointron"))
  nrow(union_tse_ara[!unionhasintron,])
  
  summerytable[, 10] <-
    c(
      nrow(tse[tsehasintron,]) * 100 / nrow(tse),
      nrow(ara[arahasintron,]) * 100 / nrow(ara),
      nrow(union_tse_ara[!unionhasintron,]) * 100 / nrow(union_tse_ara),
      nrow(intersect_tse_ara[!intersecthasintron, ]) * 100 / nrow(intersect_tse_ara)
    )
  
  #______________________ $ as pseudo frequency __________________________________
  
  tsepseudo <- tse[tse$tsenote == "pseudo", ]
  arapseudo <- ara[ara$aranote == "pseudo", ]
  intersectpseudo <-
    intersect_tse_ara[intersect_tse_ara$tsenote == "pseudo" |
                        intersect_tse_ara$aranote == "pseudo",]
  unionpseudo <-
    union_tse_ara[union_tse_ara$tsenote == "pseudo" |
                    union_tse_ara$aranote == "pseudo",]
  summerytable[, 33] <-
    c(nrow(tsepseudo),
      nrow(arapseudo),
      nrow(unionpseudo),
      nrow(intersectpseudo))
  
  
  
  #___________________ Aminoacid frequency ________________________________________
  
  aminoacid <- character(length = nrow(tse))
  aminoacid <- tse$tsefunc
  # frequency 23 classes
  tseaminoacid <- data.frame(table(aminoacid))
  summerytable$A[1] <-
    tseaminoacid[tseaminoacid$aminoacid == "A",]$Freq
  summerytable$C[1] <-
    tseaminoacid[tseaminoacid$aminoacid == "C",]$Freq
  summerytable$D[1] <-
    tseaminoacid[tseaminoacid$aminoacid == "D",]$Freq
  summerytable$E[1] <-
    tseaminoacid[tseaminoacid$aminoacid == "E",]$Freq
  summerytable$F[1] <-
    tseaminoacid[tseaminoacid$aminoacid == "F",]$Freq
  summerytable$G[1] <-
    tseaminoacid[tseaminoacid$aminoacid == "G",]$Freq
  summerytable$H[1] <-
    tseaminoacid[tseaminoacid$aminoacid == "H",]$Freq
  summerytable$I[1] <-
    tseaminoacid[tseaminoacid$aminoacid == "I",]$Freq
  summerytable$K[1] <-
    tseaminoacid[tseaminoacid$aminoacid == "K",]$Freq
  summerytable$L[1] <-
    tseaminoacid[tseaminoacid$aminoacid == "L",]$Freq
  summerytable$M[1] <-
    tseaminoacid[tseaminoacid$aminoacid == "M",]$Freq
  summerytable$N[1] <-
    tseaminoacid[tseaminoacid$aminoacid == "N",]$Freq
  summerytable$P[1] <-
    tseaminoacid[tseaminoacid$aminoacid == "P",]$Freq
  summerytable$Q[1] <-
    tseaminoacid[tseaminoacid$aminoacid == "Q",]$Freq
  summerytable$R[1] <-
    tseaminoacid[tseaminoacid$aminoacid == "R",]$Freq
  summerytable$S[1] <-
    tseaminoacid[tseaminoacid$aminoacid == "S",]$Freq
  summerytable$T[1] <-
    tseaminoacid[tseaminoacid$aminoacid == "T",]$Freq
  summerytable$V[1] <-
    tseaminoacid[tseaminoacid$aminoacid == "V",]$Freq
  summerytable$W[1] <-
    tseaminoacid[tseaminoacid$aminoacid == "W",]$Freq
  summerytable$Y[1] <-
    tseaminoacid[tseaminoacid$aminoacid == "Y", ]$Freq
  summerytable$X[1] <- 
    tseaminoacid[tseaminoacid$aminoacid == "X", ]$Freq
  summerytable$Z[1] <-
    tseaminoacid[tseaminoacid$aminoacid == "Z", ]$Freq
  
  if (length(tseaminoacid[tseaminoacid$aminoacid == "#",]$Freq) != 0)
    summerytable$`#`[1] <-
    tseaminoacid[tseaminoacid$aminoacid == "#",]$Freq
  if (length(tseaminoacid[tseaminoacid$aminoacid == "?",]$Freq) != 0)
    summerytable$`?`[1] <-
    tseaminoacid[tseaminoacid$aminoacid == "?",]$Freq
  if (length(tseaminoacid[tseaminoacid$aminoacid == "O",]$Freq) != 0)
    summerytable$O[1] <-
    tseaminoacid[tseaminoacid$aminoacid == "O",]$Freq
  
  
  
  
  #____________________________________________________________________
  
  aminoacid <- character(length = nrow(ara))
  # frequency 23 classes
  aminoacid <- ara$arafunc
  tseaminoacid <- data.frame(table(aminoacid))
  summerytable$A[2] <-
    tseaminoacid[tseaminoacid$aminoacid == "A", ]$Freq
  summerytable$C[2] <-
    tseaminoacid[tseaminoacid$aminoacid == "C", ]$Freq
  summerytable$D[2] <-
    tseaminoacid[tseaminoacid$aminoacid == "D", ]$Freq
  summerytable$E[2] <-
    tseaminoacid[tseaminoacid$aminoacid == "E", ]$Freq
  summerytable$F[2] <-
    tseaminoacid[tseaminoacid$aminoacid == "F", ]$Freq
  summerytable$G[2] <-
    tseaminoacid[tseaminoacid$aminoacid == "G", ]$Freq
  summerytable$H[2] <-
    tseaminoacid[tseaminoacid$aminoacid == "H", ]$Freq
  summerytable$I[2] <-
    tseaminoacid[tseaminoacid$aminoacid == "I", ]$Freq
  summerytable$K[2] <-
    tseaminoacid[tseaminoacid$aminoacid == "K", ]$Freq
  summerytable$L[2] <-
    tseaminoacid[tseaminoacid$aminoacid == "L", ]$Freq
  summerytable$M[2] <-
    tseaminoacid[tseaminoacid$aminoacid == "M", ]$Freq
  summerytable$N[2] <-
    tseaminoacid[tseaminoacid$aminoacid == "N", ]$Freq
  summerytable$P[2] <-
    tseaminoacid[tseaminoacid$aminoacid == "P", ]$Freq
  summerytable$Q[2] <-
    tseaminoacid[tseaminoacid$aminoacid == "Q", ]$Freq
  summerytable$R[2] <-
    tseaminoacid[tseaminoacid$aminoacid == "R", ]$Freq
  summerytable$S[2] <-
    tseaminoacid[tseaminoacid$aminoacid == "S", ]$Freq
  summerytable$T[2] <-
    tseaminoacid[tseaminoacid$aminoacid == "T", ]$Freq
  summerytable$V[2] <-
    tseaminoacid[tseaminoacid$aminoacid == "V", ]$Freq
  summerytable$W[2] <-
    tseaminoacid[tseaminoacid$aminoacid == "W", ]$Freq
  summerytable$Y[2] <-
    tseaminoacid[tseaminoacid$aminoacid == "Y", ]$Freq
  summerytable$X[2] <-
    tseaminoacid[tseaminoacid$aminoacid == "X", ]$Freq
  summerytable$Z[2] <-
    tseaminoacid[tseaminoacid$aminoacid == "Z", ]$Freq
  if (length(tseaminoacid[tseaminoacid$aminoacid == "#",]$Freq) != 0)
    summerytable$`#`[2] <-
    tseaminoacid[tseaminoacid$aminoacid == "#",]$Freq
  if (length(tseaminoacid[tseaminoacid$aminoacid == "?",]$Freq) != 0)
    summerytable$`?`[2] <-
    tseaminoacid[tseaminoacid$aminoacid == "?",]$Freq
  if (length(tseaminoacid[tseaminoacid$aminoacid == "O",]$Freq) != 0)
    summerytable$O[2] <-
    tseaminoacid[tseaminoacid$aminoacid == "O",]$Freq
  
  
  #____________________________________________________________________
  # for intersection if tse and aragorn do not match we do not consider them in our counting
  aminoacid <- character(length = nrow(intersect_tse_ara))
  aminoacid <- intersect_tse_ara$tsefunc

  # there are two undet genes 
  
  tseaminoacid <- data.frame(table(aminoacid))
  summerytable$A[4] <-
    tseaminoacid[tseaminoacid$aminoacid == "A", ]$Freq
  summerytable$C[4] <-
    tseaminoacid[tseaminoacid$aminoacid == "C", ]$Freq
  summerytable$D[4] <-
    tseaminoacid[tseaminoacid$aminoacid == "D", ]$Freq
  summerytable$E[4] <-
    tseaminoacid[tseaminoacid$aminoacid == "E", ]$Freq
  summerytable$F[4] <-
    tseaminoacid[tseaminoacid$aminoacid == "F", ]$Freq
  summerytable$G[4] <-
    tseaminoacid[tseaminoacid$aminoacid == "G", ]$Freq
  summerytable$H[4] <-
    tseaminoacid[tseaminoacid$aminoacid == "H", ]$Freq
  summerytable$I[4] <-
    tseaminoacid[tseaminoacid$aminoacid == "I", ]$Freq
  summerytable$K[4] <-
    tseaminoacid[tseaminoacid$aminoacid == "K", ]$Freq
  summerytable$L[4] <-
    tseaminoacid[tseaminoacid$aminoacid == "L", ]$Freq
  summerytable$M[4] <-
    tseaminoacid[tseaminoacid$aminoacid == "M", ]$Freq
  summerytable$N[4] <-
    tseaminoacid[tseaminoacid$aminoacid == "N", ]$Freq
  summerytable$P[4] <-
    tseaminoacid[tseaminoacid$aminoacid == "P", ]$Freq
  summerytable$Q[4] <-
    tseaminoacid[tseaminoacid$aminoacid == "Q", ]$Freq
  summerytable$R[4] <-
    tseaminoacid[tseaminoacid$aminoacid == "R", ]$Freq
  summerytable$S[4] <-
    tseaminoacid[tseaminoacid$aminoacid == "S", ]$Freq
  summerytable$T[4] <-
    tseaminoacid[tseaminoacid$aminoacid == "T", ]$Freq
  summerytable$V[4] <-
    tseaminoacid[tseaminoacid$aminoacid == "V", ]$Freq
  summerytable$W[4] <-
    tseaminoacid[tseaminoacid$aminoacid == "W", ]$Freq
  summerytable$Y[4] <-
    tseaminoacid[tseaminoacid$aminoacid == "Y", ]$Freq
  summerytable$X[4] <- 
    tseaminoacid[tseaminoacid$aminoacid == "X", ]$Freq
  summerytable$Z[4] <-
    tseaminoacid[tseaminoacid$aminoacid == "Z", ]$Freq
  tseaminoacid <- data.frame(table(aminoacid))
  
  if (length(tseaminoacid[tseaminoacid$aminoacid == "#",]$Freq) != 0)
    summerytable$`#`[4] <-
    tseaminoacid[tseaminoacid$aminoacid == "#",]$Freq
  if (length(tseaminoacid[tseaminoacid$aminoacid == "?",]$Freq) != 0)
    summerytable$`?`[4] <-
    tseaminoacid[tseaminoacid$aminoacid == "?",]$Freq
  if (length(tseaminoacid[tseaminoacid$aminoacid == "O",]$Freq) != 0)
    summerytable$O[4] <-
    tseaminoacid[tseaminoacid$aminoacid == "O",]$Freq
  
  
  
  #____________________________________________________________________
  
  aminoacid <- character(length = nrow(union_tse_ara))
  aminoacid <- union_tse_ara$tsefunc
  unionidentity <- union_tse_ara$tseidentity
  for (i in 1:nrow(union_tse_ara)) {
    if (union_tse_ara$foundby[i] == "ara")
      aminoacid[i] <- union_tse_ara$arafunc[i]
  }
  # frequency 23 classes
  tseaminoacid <- data.frame(table(aminoacid))
  summerytable$A[3] <-
    tseaminoacid[tseaminoacid$aminoacid == "A", ]$Freq
  summerytable$C[3] <-
    tseaminoacid[tseaminoacid$aminoacid == "C", ]$Freq
  summerytable$D[3] <-
    tseaminoacid[tseaminoacid$aminoacid == "D", ]$Freq
  summerytable$E[3] <-
    tseaminoacid[tseaminoacid$aminoacid == "E", ]$Freq
  summerytable$F[3] <-
    tseaminoacid[tseaminoacid$aminoacid == "F", ]$Freq
  summerytable$G[3] <-
    tseaminoacid[tseaminoacid$aminoacid == "G", ]$Freq
  summerytable$H[3] <-
    tseaminoacid[tseaminoacid$aminoacid == "H", ]$Freq
  summerytable$I[3] <-
    tseaminoacid[tseaminoacid$aminoacid == "I", ]$Freq
  summerytable$K[3] <-
    tseaminoacid[tseaminoacid$aminoacid == "K", ]$Freq
  summerytable$L[3] <-
    tseaminoacid[tseaminoacid$aminoacid == "L", ]$Freq
  summerytable$M[3] <-
    tseaminoacid[tseaminoacid$aminoacid == "M", ]$Freq
  summerytable$N[3] <-
    tseaminoacid[tseaminoacid$aminoacid == "N", ]$Freq
  summerytable$P[3] <-
    tseaminoacid[tseaminoacid$aminoacid == "P", ]$Freq
  summerytable$Q[3] <-
    tseaminoacid[tseaminoacid$aminoacid == "Q", ]$Freq
  summerytable$R[3] <-
    tseaminoacid[tseaminoacid$aminoacid == "R", ]$Freq
  summerytable$S[3] <-
    tseaminoacid[tseaminoacid$aminoacid == "S", ]$Freq
  summerytable$T[3] <-
    tseaminoacid[tseaminoacid$aminoacid == "T", ]$Freq
  summerytable$V[3] <-
    tseaminoacid[tseaminoacid$aminoacid == "V", ]$Freq
  summerytable$W[3] <-
    tseaminoacid[tseaminoacid$aminoacid == "W", ]$Freq
  summerytable$Y[3] <-
    tseaminoacid[tseaminoacid$aminoacid == "Y", ]$Freq
  summerytable$X[3] <- 
    tseaminoacid[tseaminoacid$aminoacid == "X", ]$Freq
  summerytable$Z[3] <-
    tseaminoacid[tseaminoacid$aminoacid == "Z", ]$Freq
  if (length(tseaminoacid[tseaminoacid$aminoacid == "#",]$Freq) != 0)
    summerytable$`#`[3] <-
    tseaminoacid[tseaminoacid$aminoacid == "#",]$Freq
  if (length(tseaminoacid[tseaminoacid$aminoacid == "?",]$Freq) != 0)
    summerytable$`?`[3] <-
    tseaminoacid[tseaminoacid$aminoacid == "?",]$Freq
  if (length(tseaminoacid[tseaminoacid$aminoacid == "O",]$Freq) != 0)
    summerytable$O[3] <-
    tseaminoacid[tseaminoacid$aminoacid == "O",]$Freq
  
  summerytable
  
  # formating the table and saving it into a file with fixed length format
  n <- data.frame(
    "GeneSet",
    "#tRNA",
    "#nucleotides",
    "N/T",
    "lengthrange",
    "%G",
    "%C",
    "%T",
    "%A",
    "%intron",
    "A",
    "C",
    "D",
    "E",
    "F",
    "G",
    "H",
    "I",
    "K",
    "L",
    "M",
    "N",
    "P",
    "Q",
    "R",
    "S",
    "T",
    "V",
    "W",
    "Y",
    "X",
    "Z",
    "$",
    "?",
    "#",
    "O"
  )
  names(n) <-  c(
    "GeneSet",
    "#tRNA",
    "#nucleotides",
    "N/T",
    "lengthrange",
    "%G",
    "%C",
    "%T",
    "%A",
    "%intron",
    "A",
    "C",
    "D",
    "E",
    "F",
    "G",
    "H",
    "I",
    "K",
    "L",
    "M",
    "N",
    "P",
    "Q",
    "R",
    "S",
    "T",
    "V",
    "W",
    "Y",
    "X",
    "Z",
    "$",
    "?",
    "#",
    "O"
  )
  
  names(summerytable) <- names(n)
  for (i in 1:ncol(summerytable)) {
    summerytable[, i] <- as.character(summerytable[, i])
  }
  for (i in 6:ncol(summerytable)) {
    summerytable[, i] <- format(as.numeric(summerytable[, i]), digits = 4)
  }
  summerytable[, 4] <-
    format(as.numeric(summerytable[, 4]), digits = 4)
  
  for (i in 1:ncol(summerytable)) {
    summerytable[, i] <- as.character(summerytable[, i])
  }
  
  write.fwf(
    rbind(n, summerytable),
    width = c(
      12,
      7,
      12,
      7,
      11,
      5,
      5,
      5,
      5,
      10,
      4,
      3,
      3,
      3,
      3,
      3,
      3,
      3,
      3,
      3,
      3,
      3,
      3,
      3,
      3,
      3,
      3,
      3,
      3,
      3,
      3,
      3,
      3,
      3,
      3,
      3
    ),
    colnames = FALSE,
    file = paste(outputpath,"Gene_Summary_Table.txt",sep = "")
  )
  summerytable
  #paste(summerytable[4,],collapse = " & ")
}
endDisplacement <- function(integrated_tse_ara,outputpath){
  library(ggplot2)
  library(dplyr)
  isboth <- integrated_tse_ara$foundby == "both"
  tseonly <- integrated_tse_ara$foundby == "tse"
  araonly <- integrated_tse_ara$foundby == "ara"
  
  istse <- (tseonly | isboth)
  isara <- (araonly | isboth)
  
  tse <- integrated_tse_ara[istse,]
  ara <- integrated_tse_ara[isara,]
  
  # for intersection if tse and aragorn do not match we do not consider them in our counting
  both <- integrated_tse_ara[isboth,]
  #both <- both[both$tsefunc == both$arafunc, ]
  union_tse_ara <- integrated_tse_ara
  
  Fiveprimend = integer(length = nrow(both))
  Threeprimend = integer(length = nrow(both))
  
  for (i in 1:nrow(both)) {
    if (both$direction[i] == "+")
    {
      Fiveprimend[i] = as.integer(both$tsebegin[i]) - as.integer(both$arabegin[i])
      Threeprimend[i] = as.integer(both$tseend[i]) - as.integer(both$araend[i])
    }
    if (both$direction[i] == "-")
    {
      Threeprimend[i] = (as.integer(both$tsebegin[i]) - as.integer(both$arabegin[i])) *
        -1
      Fiveprimend[i] = (as.integer(both$tseend[i]) - as.integer(both$araend[i])) *
        -1
    }
  }
  
  bad <- is.na(Fiveprimend)
  Fiveprimend <- Fiveprimend[!bad]
  #Difference between tse1end and tse2end
  
  bad <- is.na(Threeprimend)
  Threeprimend <- Threeprimend[!bad]
  
  par(mfrow = c(1, 1))
  X = Fiveprimend
  Y = Threeprimend
  data <- data.frame(X, Y)
  data = group_by(data, X, Y)
  data = summarize(data, frequency = n())
  data$frequency = data$frequency #/ sims
  xbreaks <- as.integer(levels(factor(X)))
  ybreaks <- as.integer(levels(factor(Y)))
  total <- sum(data$frequency)
  plottitle <-
    paste(
      "Empirical joint distribution of end-displacements in \n",
      nrow(both),
      " tRNA gene models found by both Aragorn and tRNAscan-SE 2.0",
      sep = ""
    )
  
  p <- ggplot(data = data, aes(X, Y)) +
    geom_tile(aes(fill = frequency), color = "white") + geom_text(aes(label = frequency)) +
    scale_fill_gradient(low = "lightblue", high = "darkred", name = "Frequency") +
    labs(title = plottitle) +
    theme(plot.margin = unit(c(1.8, .5, 1.75, 1.55), "cm")) + theme(plot.title = element_text(
      color = "#383838",
      face = "bold",
      size = 17,
      vjust = 4,
      hjust = 0.5
    )) +
    theme(
      axis.title.x = element_text(
        color = "#383838",
        face = "bold",
        size = 13,
        vjust = -1.5
      ),
      axis.title.y = element_text(
        color = "#383838",
        face = "bold",
        size = 13,
        vjust = 2
      )
    ) +
    theme(
      axis.text.x = element_text(face = "bold",
                                 color = "#383838",
                                 size = 10),
      axis.text.y = element_text(face = "bold",
                                 color = "#383838",
                                 size = 10)
    ) +
    xlab("5-prime displacement (tSE2-Ara)") + ylab("3-prime displacement (tSE2-Ara)") + scale_y_continuous(breaks = ybreaks) +
    scale_x_continuous(breaks = xbreaks) 
  ggsave(paste(outputpath, "EndDisplacement.png", sep = ""),
         width = 14,
         height = 7)
}


#D|I L|? L|E L|M N|Y O|M S|R W|G 
#5   3   1   9  11   2   1   3 