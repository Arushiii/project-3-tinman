library(tidyverse)
library(dplyr)
library(gridExtra)


star1 <- as.data.frame(t(read.delim("star_results2/SRR1177960_star_Log.final.out")))
star2 <- as.data.frame(t(read.delim("star_results2/SRR1177967_star_Log.final.out")))
star3 <- as.data.frame(t(read.delim("star_results2/SRR1177968_star_Log.final.out")))
star4 <- as.data.frame(t(read.delim("star_results2/SRR1177971_star_Log.final.out")))
star5 <- as.data.frame(t(read.delim("star_results2/SRR1177984_star_Log.final.out")))
star6 <- as.data.frame(t(read.delim("star_results2/SRR1177985_star_Log.final.out")))
star7 <- as.data.frame(t(read.delim("star_results2/SRR1177986_star_Log.final.out")))
star8 <- as.data.frame(t(read.delim("star_results2/SRR1178023_star_Log.final.out")))
star9 <- as.data.frame(t(read.delim("star_results2/SRR1178049_star_Log.final.out")))

##################### could not get this for-loop to work#####################
# stars <- list(star1,star2,star3,star4,star5,star6,star7,star8,star9)
#
# for (val in stars)
# {
#   names(val) <- as.matrix(val[1, ])
#   val <- val[-1, ]
#   val[] <- lapply(val, function(x) type.convert(as.character(x)))
# }
##############################################################################

names(star1) <- as.matrix(star1[1, ])
star1 <- star1[-1, ]
star1[] <- lapply(star1, function(x) type.convert(as.character(x)))
star1 <- star1 %>% mutate(sample="SRR1177960") 

names(star2) <- as.matrix(star2[1, ])
star2 <- star2[-1, ]
star2[] <- lapply(star2, function(x) type.convert(as.character(x)))
star2 <- star2 %>% mutate(sample="SRR1177967") 

names(star3) <- as.matrix(star3[1, ])
star3 <- star3[-1, ]
star3[] <- lapply(star3, function(x) type.convert(as.character(x)))
star3 <- star3 %>% mutate(sample="SRR1177968") 

names(star4) <- as.matrix(star4[1, ])
star4 <- star4[-1, ]
star4[] <- lapply(star4, function(x) type.convert(as.character(x)))
star4 <- star4 %>% mutate(sample="SRR1177971")

names(star5) <- as.matrix(star5[1, ])
star5 <- star5[-1, ]
star5[] <- lapply(star5, function(x) type.convert(as.character(x)))
star5 <- star5 %>% mutate(sample="SRR1177984") 

names(star6) <- as.matrix(star6[1, ])
star6 <- star6[-1, ]
star6[] <- lapply(star6, function(x) type.convert(as.character(x)))
star6 <- star6 %>% mutate(sample="SRR1177985") 

names(star7) <- as.matrix(star7[1, ])
star7 <- star7[-1, ]
star7[] <- lapply(star7, function(x) type.convert(as.character(x)))
star7 <- star7 %>% mutate(sample="SRR1177986") 

names(star8) <- as.matrix(star8[1, ])
star8 <- star8[-1, ]
star8[] <- lapply(star8, function(x) type.convert(as.character(x)))
star8 <- star8 %>% mutate(sample="SRR1178023") 

names(star9) <- as.matrix(star9[1, ])
star9 <- star9[-1, ]
star9[] <- lapply(star9, function(x) type.convert(as.character(x)))
star9 <- star9 %>% mutate(sample="SRR1178049") 


allstar1 <- full_join(star1,star2)
allstar2  <- full_join(star3,star4)
allstar3  <- full_join(star5,star6)
allstar4 <- full_join(star7,star8)
allstar_a <- full_join(allstar1,allstar2)
allstar_b <- full_join(allstar3,allstar4)
allstar <- full_join(allstar_a,allstar_b)
allstar <- full_join(star9, allstar)
allstar_sum <- allstar %>% select(sample, 7, 22, 27,28,29)
allstar_sum
star_table <- write.table(allstar_sum, "./star_table.txt", sep="\t")

png("star_table.png", height = 25*nrow(allstar_sum), width = 240*ncol(allstar_sum))
grid.table(allstar_sum)
dev.off()

allstar %>% mutate("read mean" =as.logical("Number of input reads |"))
