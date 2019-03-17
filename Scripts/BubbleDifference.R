BubbleDifference <- function() {
  TableDir <-
    "/home/fatemeh/Leishmania_2019/Leishmania_2019/Results/tsfmInput-output/output2/Logos/"
  dirpath <- TableDir
  clusterdir <- list.dirs(path = dirpath, recursive = FALSE)
  #_______________________________________________________________
  # run for the first cluster
  splitedpath <- unlist(strsplit(clusterdir[1], split = "/"))
  clus_name <- splitedpath[length(splitedpath)]
  if (clus_name == "HOMO")
    next
  tablepath <-
    paste(TableDir,
          clus_name,
          "/Bubble/",
          clus_name,
          "_Table.txt",
          sep = "")
  tabledf <- read.table(tablepath, header = TRUE)
  tabledf <- calculatediff(tabledf)
  tabledf$cluster <- rep(clus_name, nrow(tabledf))
  #_________________________________________________________________
  # run for next clusters and bind them to the main tabledf
  for (i in 2:length(clusterdir)) {
    splitedpath <- unlist(strsplit(clusterdir[i], split = "/"))
    clus_name <- splitedpath[length(splitedpath)]
    if (clus_name == "HOMO")
      next
    tablepath <-
      paste(TableDir,
            clus_name,
            "/Bubble/",
            clus_name,
            "_Table.txt",
            sep = "")
    df <- read.table(tablepath, header = TRUE)
    df <- calculatediff(df)
    df$cluster <- rep(clus_name, nrow(df))
    
    tabledf <- rbind(tabledf, df)
  }
  tabledf$importance <- tabledf$fbits * tabledf$fht
  visualize.diff(tabledf)
}


calculatediff <- function(df) {
  gains <- (df$gainbits * df$gainfht)
  convs <- (df$convbits * df$convfht)
  colormap <- function (g, c) {
    y <- rep(0, length(g))
    y[g <  0.48            & c < 0.44]             <- "#D1E5F0"
    y[g >= 0.48 & g < 0.95 & c < 0.44]             <- "#92C5DE"
    y[g >= 0.95            & c < 0.44]             <- "#4393C3"
    y[g <  0.48            & c >= 0.44 & c < 0.70] <- "#F4A582"
    y[g >= 0.48 & g < 0.95 & c >= 0.44 & c < 0.70] <- "#D6604D"
    y[g >= 0.95            & c >= 0.44 & c < 0.70] <- "#B2182B"
    y[g <  0.48            & c >= 0.70]            <- "#A6DBA0"
    y[g >= 0.48 & g < 0.95 & c >= 0.70]            <- "#5AAE61"
    y[g >= 0.95            & c >= 0.70]            <- "#1B7837"
    y
    
  }
  colors <- colormap(gains, convs)
  df$colors <- colors
  df
}
visualize.diff <- function(df) {
  library(magrittr)
  library(ggplot2)
  library(plyr)
  library(ggrepel)
  library(gridExtra)
  # make another column to number the clusters
  clusternames <-
    as.character(as.data.frame(table(df$cluster))$Var1)
  clusterIndex <-
    data.frame(clusternames, seq(1, length(clusternames), 1))
  names(clusterIndex) <- c("name", "index")
  df$clusterIndex <- 0
  for (i in 1:nrow(df)) {
    df$clusterIndex[i] <-
      clusterIndex[clusterIndex$name == df$cluster[i],]$index
  }
  
  # Create a column 'smartlabels' with importance value that are less than or equal to 0.5 as blank ("") and use that in the label argument
  df <- df %>%
    mutate(smartlabels = replace(clusterIndex, importance <= 0.5, ""))
  
  # remove rows that their importance is zero
  # df <- df[df$importance != 0, ]
  
  set.seed(42)
  filepath <-
    "/home/fatemeh/Leishmania_2019/Leishmania_2019/Results/tsfmInput-output/output2/BubblePlots/"

  for (class in levels(df$aa)) {
    filenm <- paste(filepath, "Differences_", class, ".pdf", sep = "")
    plotlist <- list()
    count = 1
    model_df <- df[df$aa == class,]
    for (stat in levels(df$state)) {
      ylabel <- paste("Functional Importance for state", stat)
      model_state_df <- model_df[model_df$state == stat,]
      # Sort by coord and cluster
      df_sorted <- arrange(model_state_df, coord, cluster)
      # Calculate the cumulative sum of importance for each coord
      df_cumsum <- ddply(df_sorted, "coord",
                         transform, label_ypos = cumsum(importance))
      
      p <-
        ggplot(data = df_cumsum, aes(x = as.factor(df_cumsum$coord), y = importance)) +
        geom_bar(stat = "identity",
                 fill = df_cumsum$colors,
                 colour = "black",width = 0.95) +
        geom_text(
          aes(y = label_ypos, label = smartlabels),
          vjust = 1.5,
          color = "black",
          size = 3.5
        ) + xlab("coordinate") + ylab(ylabel) +
        theme(
          axis.text.x = element_text(
            angle = 90,
            hjust = 1,
            vjust = 0.5,
            size = 4.5
          ),
          axis.title.x = element_blank()
        ) +
        scale_fill_brewer(palette = "Paired")+ ylim(0, 30) #+ theme_minimal()
      plotlist[[count]] <-  p
      count <- count + 1
    }
    ml <- do.call("grid.arrange", c(plotlist, ncol = 2))
    ggsave(ml,
           width = 20,
           height = 15,
           filename = filenm)
  }
  # it gave an error whrn I removed lines with zero importance
  #do.call("grid.arrange", c(plotlist, ncol=3))
  #ml <- multiplot(plotlist,ncol=1)
  #ml
}
