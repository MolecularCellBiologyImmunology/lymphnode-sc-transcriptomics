library(dplyr)

Douwe <- read.csv("LNS_ALL.csv",sep=",",header=TRUE)
Douwe[Douwe == 0] <- NA
Douwe$Average <- rowMeans(Douwe[,-1], na.rm = TRUE)
Average1 <- data.frame("geneid" = Douwe$geneid, "Douwe AVG" = Douwe$Average)

Freiburg <- read.csv("LNS_W_ALL_GENECOUNTS_UNFILTERED.csv",sep=",",header=TRUE)
Freiburg[Freiburg == 0] <- NA
Freiburg$Average <- rowMeans(Freiburg[,-1], na.rm = TRUE)
Average2 <- data.frame("geneid" = Freiburg$GENEID, "Freiburg AVG" = Freiburg$Average)

Average <- full_join(Average1, Average2, by = "geneid")
