I_ LIBRARY 
library(qiime2R)
library(dplyr)
library(tidyr)
library(phyloseq)
library(microViz)
library(tidyverse)
library(microbiome)
library(vegan)
library(picante)
library(ALDEx2)
library(ggplot2)
library(dendextend)
library(MicrobiotaProcess)
#####################################################################

#Donnée
MD_Global = read.table(file = "MD_env_VFA_Health.txt", header = TRUE, sep = "\t", row.names = NULL, stringsAsFactors = FALSE)
Count_Genus_Global = read.table(file = "Count_Genus_Rumen.txt", header = TRUE, sep = "\t", row.names = NULL, stringsAsFactors = FALSE)
Count_Phylum_Global = read.table(file = "Count_Phylum_Rumen.txt", header = TRUE, sep = "\t", row.names = NULL, stringsAsFactors = FALSE)

#Phylumdata 
Phylum_env<-full_join(MD_Global,Count_Phylum_Global)
Phylum_env<- na.omit(Phylum_env)



LDAPhylum_Mixed_T3 = read.table(file = "Count_T3_Mixed_LEFSE.txt", header = TRUE, sep = "\t", row.names = NULL, stringsAsFactors = FALSE)
LDAPhylum_Control_T3 = read.table(file = "Count_T3_Control_LEFSE.txt", header = TRUE, sep = "\t", row.names = NULL, stringsAsFactors = FALSE)
LDAPhylum_Dam_T10 = read.table(file = "Count_T10_Dam_LEFSE.txt", header = TRUE, sep = "\t", row.names = NULL, stringsAsFactors = FALSE)
LDAPhylum_Mixed_T10 = read.table(file = "Count_T10_Mixed_LEFSE.txt", header = TRUE, sep = "\t", row.names = NULL, stringsAsFactors = FALSE)
LDAPhylum_Control_T10 = read.table(file = "Count_T10_Control_LEFSE.txt", header = TRUE, sep = "\t", row.names = NULL, stringsAsFactors = FALSE)
LDAPhylum_Mixed_T13 = read.table(file = "Count_T13_Mixed_LEFSE.txt", header = TRUE, sep = "\t", row.names = NULL, stringsAsFactors = FALSE)
LDAPhylum_Control_T13 = read.table(file = "Count_T13_Control_LEFSE.txt", header = TRUE, sep = "\t", row.names = NULL, stringsAsFactors = FALSE)

#Phylum_Time13
Phylum_Mixed_3<-full_join(MD_Global,LDAPhylum_Mixed_T3)
Phylum_Control_3<-full_join(MD_Global,LDAPhylum_Control_T3)
Phylum_Dam_10<-full_join(MD_Global,LDAPhylum_Dam_T10)
Phylum_Mixed_10<-full_join(MD_Global,LDAPhylum_Mixed_T10 )
Phylum_Control_10<-full_join(MD_Global,LDAPhylum_Control_T10)
Phylum_Mixed_13<-full_join(MD_Global,LDAPhylum_Mixed_T13)
Phylum_Control_13<-full_join(MD_Global,LDAPhylum_Control_T13)

Phylum_Mixed_3 <- subset(Phylum_Mixed_3, Phylum_Mixed_3$Time_cat == "3")
Phylum_Mixed_3<-Phylum_Mixed_3 [,-1:-27] # supp 
Phylum_Control_3 <- subset(Phylum_Control_3, Phylum_Control_3$Time_cat == "3")
Phylum_Control_3<-Phylum_Control_3 [,-1:-27] # supp 
Phylum_Dam_10 <- subset(Phylum_Dam_10, Phylum_Dam_10$Time_cat == "10")
Phylum_Dam_10<-Phylum_Dam_10 [,-1:-27] # supp 
Phylum_Control_10 <- subset(Phylum_Control_10, Phylum_Control_10$Time_cat == "10")
Phylum_Control_10<-Phylum_Control_10 [,-1:-27] # supp 
Phylum_Mixed_10 <- subset(Phylum_Mixed_10, Phylum_Mixed_10$Time_cat == "10")
Phylum_Mixed_10<-Phylum_Mixed_10 [,-1:-27] # supp 
Phylum_Control_13 <- subset(Phylum_Control_13, Phylum_Control_13$Time_cat == "13")
Phylum_Control_13<-Phylum_Control_13 [,-1:-27] # supp 
Phylum_Mixed_13 <- subset(Phylum_Mixed_13, Phylum_Mixed_13$Time_cat == "13")
Phylum_Mixed_13<-Phylum_Mixed_13 [,-1:-27] # supp 



testRes_Phylum_Mixed_3 = cor.mtest(Phylum_Mixed_3, conf.level = 0.95)
testRes_Phylum_Control_3 = cor.mtest(Phylum_Control_3, conf.level = 0.95)
testRes_Phylum_Dam_10 = cor.mtest(Phylum_Dam_10, conf.level = 0.95)
testRes_Phylum_Control_10 = cor.mtest(Phylum_Control_10, conf.level = 0.95)
testRes_Phylum_Mixed_10 = cor.mtest(Phylum_Mixed_10, conf.level = 0.95)
testRes_Phylum_Control_13 = cor.mtest(Phylum_Control_13, conf.level = 0.95)
testRes_Phylum_Mixed_13 = cor.mtest(Phylum_Mixed_13, conf.level = 0.95)



M_Phylum_Mixed_3<-cor(Phylum_Mixed_3)
M_Phylum_Control_3<-cor(Phylum_Control_3)
M_Phylum_Dam_10<-cor(Phylum_Dam_10)
M_Phylum_Control_10<-cor(Phylum_Control_10)
M_Phylum_Mixed_10<-cor(Phylum_Mixed_10)
M_Phylum_Control_13<-cor(Phylum_Control_13)
M_Phylum_Mixed_13<-cor(Phylum_Mixed_13)


par(mfrow=c(2,3))
corrplot(M_Phylum_Mixed_3, p.mat = testRes_Phylum_Mixed_3, method = 'circle', type = 'lower', insig='blank',sig.level=0.05,
         number.cex = 0.8, diag=FALSE,tl.col = 'black')

corrplot(Phylum_Control_3, p.mat = testRes_Phylum_Control_3, method = 'circle', type = 'lower', insig='blank',sig.level=0.05,
         number.cex = 0.8, diag=FALSE,tl.col = 'black')

par(mfrow=c(2,3))

corrplot(M_Phylum_Control_10, p.mat = testRes_Phylum_Control_10, method = 'circle', type = 'lower', insig='blank',sig.level=0.05,
         number.cex = 0.8, diag=FALSE,tl.col = 'black')

corrplot(M_Phylum_Mixed_10, p.mat = testRes_Phylum_Mixed_10, method = 'circle', type = 'lower', insig='blank',sig.level=0.05,
         number.cex = 0.8, diag=FALSE,tl.col = 'black')

corrplot(M_Phylum_Dam_10, p.mat = testRes_Phylum_Dam_10, method = 'circle', type = 'lower', insig='blank',sig.level=0.05,
         number.cex = 0.8, diag=FALSE,tl.col = 'black')

corrplot(M_Phylum_Control_13, p.mat = testRes_Phylum_Control_13, method = 'circle', type = 'lower', insig='blank',sig.level=0.05,
         number.cex = 0.8, diag=FALSE,tl.col = 'black')

corrplot(M_Phylum_Mixed_13, p.mat = testRes_Phylum_Mixed_13, method = 'circle', type = 'lower', insig='blank',sig.level=0.05,
         number.cex = 0.8, diag=FALSE,tl.col = 'black')




#Genusdata 
count_LDAMIXED13 = read.table(file = "Count_Genus_Marker_TIME13_mixed.txt", header = TRUE, sep = "\t", row.names = NULL, stringsAsFactors = FALSE)
count_LDAcontrol13 = read.table(file = "Count_Genus_Marker_TIME13_control.txt", header = TRUE, sep = "\t", row.names = NULL, stringsAsFactors = FALSE)
count_LDADam13 = read.table(file = "Count_Genus_Marker_TIME13_Dam.txt", header = TRUE, sep = "\t", row.names = NULL, stringsAsFactors = FALSE)
count_LDADam10 = read.table(file = "Count_Genus_Marker_TIME10_Dam.txt", header = TRUE, sep = "\t", row.names = NULL, stringsAsFactors = FALSE)
count_LDADMIXED10 = read.table(file = "Count_Genus_Marker_TIME10_Mixed.txt", header = TRUE, sep = "\t", row.names = NULL, stringsAsFactors = FALSE)
count_LDcontrol10= read.table(file = "Count_Genus_Marker_TIME10_control.txt", header = TRUE, sep = "\t", row.names = NULL, stringsAsFactors = FALSE)


#Genus_Time13
Genus_Mixed_13<-full_join(MD_Global,count_LDAMIXED13)
Genus_control_13<-full_join(MD_Global,count_LDAcontrol13)
Genus_Dam_13<-full_join(MD_Global,count_LDADam13)
Genus_Dam_10<-full_join(MD_Global,count_LDADam10)
Genus_Mixed_10<-full_join(MD_Global,count_LDADMIXED10)
Genus_control_10<-full_join(MD_Global,count_LDcontrol10)


Genus_Mixed_13 <- subset(Genus_Mixed_13, Genus_Mixed_13$Time_cat == "13")
Genus_Mixed_13<-Genus_Mixed_13 [,-1:-27] # supp 

Genus_control_13 <- subset(Genus_control_13, Genus_control_13$Time_cat == "13")
Genus_control_13<-Genus_control_13 [,-1:-27] # supp 

Genus_Dam_13 <- subset(Genus_Dam_13, Genus_Dam_13$Time_cat == "13")
Genus_Dam_13<-Genus_Dam_13 [,-1:-27] # supp 

Genus_Dam_10 <- subset(Genus_Dam_10, Genus_Dam_10$Time_cat == "10")
Genus_Dam_10<-Genus_Dam_10 [,-1:-27] # supp 

Genus_Mixed_10 <- subset(Genus_Mixed_10, Genus_Mixed_10$Time_cat == "10")
Genus_Mixed_10<-Genus_Mixed_10 [,-1:-27] # supp 

Genus_control_10 <- subset(Genus_control_10, Genus_control_10$Time_cat == "10")
Genus_control_10<-Genus_control_10 [,-1:-27] # supp 

###########################################################################

#Marker LEFSE : 
MARKER13PHYLUM<-full_join(MD_Global,Count_LDAT13Phylum)
MARKER10PHYLUM<- full_join(MD_Global,Count_LDAT10Phylum)
MARKER3PHYLUM<- full_join(MD_Global,Count_LDAT3Phylum)

#Time13

MARKT13_Mixed <- subset(MARKER13PHYLUM, MARKER13PHYLUM$Lot_x_Time == "panache_13")
MARKT13_Mixed<-MARKT13_Mixed [,-1:-27] # supp 
MARKT13_DAM <- subset(MARKER13PHYLUM, MARKER13PHYLUM$Lot_x_Time == "mere_13")
MARKT13_DAM<-MARKT13_DAM [,-1:-27] # supp 
MARKT13_Control <- subset(MARKER13PHYLUM, MARKER13PHYLUM$Lot_x_Time == "temoin_13")
MARKT13_Control<-MARKT13_Control [,-1:-27] # supp 

#Time3

MARKT3_Mixed <- subset(MARKER3PHYLUM, MARKER13PHYLUM$Lot_x_Time == "panache_3")
MARKT3_Mixed<-MARKT3_Mixed [,-1:-27] # supp 
MARKT3_DAM <- subset(MARKER3PHYLUM, MARKER13PHYLUM$Lot_x_Time == "mere_3")
MARKT3_DAM<-MARKT3_DAM [,-1:-27] # supp 
MARKT3_Control <- subset(MARKER3PHYLUM, MARKER13PHYLUM$Lot_x_Time == "temoin_3")
MARKT3_Control<-MARKT3_Control [,-1:-27] # supp 


str(MARKT13)
str(MARKT10)
str(MARKT3)

#################################################################
### Correlation indicator OTUs
library(corrplot)
# mat : is a matrix of data
# ... : further arguments to pass to the native R cor.test function
cor.mtest <- function(mat, ...) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat<- matrix(NA, n, n)
  diag(p.mat) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j], ...)
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
    }
  }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  p.mat
}
##############################################################
testRes_Mixed13 = cor.mtest(Genus_Mixed_13, conf.level = 0.95)
testRes_Control13 = cor.mtest(Genus_control_13, conf.level = 0.95)
testRes_Dam13 = cor.mtest(Genus_Dam_13, conf.level = 0.95)

testRes_Mixed10 = cor.mtest(Genus_Mixed_10, conf.level = 0.95)
testRes_Control10 = cor.mtest(Genus_control_10, conf.level = 0.95)
testRes_Dam10 = cor.mtest(Genus_Dam_10, conf.level = 0.95)

M_MixedG13<-cor(Genus_Mixed_13)
M_Control13<-cor(Genus_control_13)
M_Dam13<-cor(Genus_Dam_13)

M_MixedG10<-cor(Genus_Mixed_10)
M_Control10<-cor(Genus_control_10)
M_Dam10<-cor(Genus_Dam_10)

par(mfrow=c(1,1))

corrplot(M_MixedG13, p.mat = testRes_Mixed13, method = 'circle', type = 'lower', insig='blank',sig.level=0.05,
         number.cex = 0.8, diag=FALSE,tl.col = 'black')

corrplot(M_Control13, p.mat = testRes_Control13, method = 'circle', type = 'lower', insig='blank',sig.level=0.05,
         number.cex = 0.8, diag=FALSE,tl.col = 'black')

corrplot(M_Dam13, p.mat = testRes_Dam13, method = 'circle', type = 'lower', insig='blank',sig.level=0.05,
         number.cex = 0.8,diag=FALSE,tl.col = 'black')

corrplot(M_MixedG10, p.mat = testRes_Mixed10, method = 'circle', type = 'lower', insig='blank',sig.level=0.05,
         number.cex = 0.8, diag=FALSE,tl.col = 'black')

corrplot(M_Control10, p.mat = testRes_Control10, method = 'circle', type = 'lower', insig='blank',sig.level=0.05,
         number.cex = 0.8, diag=FALSE,tl.col = 'black')

corrplot(M_Dam10, p.mat = testRes_Dam10, method = 'circle', type = 'lower', insig='blank',sig.level=0.05,
         number.cex = 0.8,diag=FALSE,tl.col = 'black')

#############################################################

############################################################################



corrplot(M_Dam13[1:8,9:29], p.mat = testRes_Dam13[1:8,9:29], method = 'circle', type = 'lower', insig='blank',sig.level=0.05,
         number.cex = 0.8, diag=FALSE)














