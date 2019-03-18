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