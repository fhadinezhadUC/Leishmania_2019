# read all the genomes from file
# for each genome count number of sequences by counting number of headlines
# plot a diagram with x= name of the genome, seqcountcomparison = (number of sequence/length of sequence)*scale factor of (1/max(seqcountcomparison))
genomecomparison <- function() {
  library(gsubfn)
  library(ggplot2)
  filenames <-
    list.files(
      "/home/fatemeh/Leishmania_Aug2018/tritrypdb.org/",
      pattern = "*.fasta",
      full.names = TRUE
    )
  seqcountcomp <- integer(length = length(filenames))
  genomelength <- integer(length = length(filenames))
  for (i in 1:length(filenames)) {
    headers <- read.pattern(filenames[i], "^>+.*")
    genomelen = 0
    myseqlen = 0
    maxnum = 0
    # extract the length of the sequence from this line by taking the field after "to"
    # sum the length of all sequences as the length of the genome
    for (r in 1:nrow(headers)) {
      fields <- unlist(strsplit(as.character(headers[r, 1]), " "))
      for (j in 1:length(fields)) {
        if (length(grep("length=*", fields[j])) != 0)
        {
          myseqlen <- substring(fields[j], 8)
          myseqlen <- as.integer(myseqlen)
          break
        }
      }
      genomelen <- genomelen + myseqlen
    }
    genomelength[i] <- genomelen
    seqcountcomp[i] <- nrow(headers) / genomelen
    if (maxnum < seqcountcomp[i])
      maxnum <- seqcountcomp[i]
  }
  score <- maxnum / seqcountcomp
  # the higher the bar the better genome is sequenced.
  genomenames <- filenames
  for (i in 1:length(filenames)) {
    temp <- unlist(strsplit(filenames[i], split = "/"))
    temp <- temp[length(temp)]
    temp <- gsub(temp, pattern = "_Genome.fasta", replacement = "")
    temp <- gsub(temp, pattern = "TriTrypDB-38_", replacement = "")
    genomenames[i] <- temp
  }
  
  info <- data.frame(genomenames, score)
  names(info) <- c("genome", "score")
  
  info <- info[order(info$score), ]
  info$genome <-
    factor(info$genome, levels = info$genome)
  library(ggplot2)
  xname = paste(round(maxnum, digits = 5),
                "/ (Number of sequences/length of genome)")
  p <- ggplot(data = info, aes(x = genome, y = score)) +
    geom_bar(stat = "identity", fill = "steelblue") + theme(axis.text.x = element_text(angle = 90, hjust = 1),
                                                            plot.title = element_text(hjust = 0.5)) +
    labs(title =
           "TryTryp Version 38 Genome Comparison",
         x =
           "Genome", y = xname)
  p
  ggsave(
    paste(
      "/home/fatemeh/Leishmania_2019/Leishmania_2019/Document/Figures/",
      "GenomeComparisonV38.png",
      sep = ""
    ),
    width = 14,
    height = 7
  )
  
  
}