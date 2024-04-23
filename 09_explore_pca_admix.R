###############
### Pcangsd ###
###############

# import bams used
bamUse <-read.table("bamList_625_Ahya", header=FALSE, sep=",", stringsAsFactors=FALSE)

colnames(bamUse)[1] <- c("Code") 

# any duplicated samples?
bamUse[duplicated(bamUse$Code),]

# if no metadata, then 
dat <- bamUse

# add color
library(RColorBrewer)
library(scales)

col <-brewer.pal(10, "BrBG")

library(data.table)
dat$col <- "black"

###############
## PCA
###############

library(RcppCNPy)

# Reads in estimated covariance matrix
cov <- as.matrix(read.table("./data/pca_outQ5_ECT_625_chrs.cov")) 

dim(cov)

e <- eigen(cov)

(e$values[1]/sum(e$values))

# proportion of total variation for each PC
for (i in e$values) {
  print(i/sum(e$values))
}

# first four PC axes
e$values[1]/sum(e$values)*100
e$values[2]/sum(e$values)*100
e$values[3]/sum(e$values)*100
e$values[4]/sum(e$values)*100

### plot

par(mfrow=c(1,1))

plot(e$vectors[,1:2],lwd=2,  xlab=paste0(c("PC1" , signif(e$values[1]/sum(e$values)*100, 3), sep="")), 
     ylab=paste0(c("PC2" , signif(e$values[2]/sum(e$values)*100, 3), sep="")), 
     pch=16, col=alpha(dat$col, 1), cex=1.5, main="")

# color by cluster 
# "Group2" = "#00a1d5","Group3" = "#b24745","Group4" = "#df8f44","Group1" = "#374e55"

# Ahya 625 assign group colors
dat[which((e$vectors[,1] > -0.03)), ]$col="#374e55"
dat[which((e$vectors[,2] > 0.11)), ]$col="#00a1d5"
dat[which((e$vectors[,2] < 0.11 & e$vectors[,2] > 0.05)), ]$col="#b24745"
dat[which((e$vectors[,2] < 0.0 & e$vectors[,1] < -0.06)), ]$col="#df8f44"

# replot with colors
plot(e$vectors[,1:2],lwd=2,  xlab="PC1  (4.64%)", ylab="PC2  (1.14%)", 
     pch=16, col=alpha(dat$col, 0.7), cex=2, main="")

################
## ADMIXUTRE
################

# assign groups
dat$Group <- NA
head(dat)

# Ahya 625
dat[which((e$vectors[,1] > -0.03)), ]$Group=1
dat[which((e$vectors[,2] > 0.11)), ]$Group=2
dat[which((e$vectors[,2] < 0.11 & e$vectors[,2] > 0.05)), ]$Group=3
dat[which((e$vectors[,2] < 0.0 & e$vectors[,1] < -0.06)), ]$Group=4

# plot first
# the input here is a matrix where rows will be the qvalues and the columns the samples

par(mfrow=c(4,1), mar = c(4,4,2,2))

# Ahya 625  
K <- as.matrix(read.table("./data/pca_outQ5_ECT_625_chrs_2.admix.2.Q"))
dim(K)

# transverse-> order -> transverse back 
x <- t(as.matrix(K))
xOrd<- x[ ,order(dat$Group)]
xOrd<-t(xOrd)

# plot
barplot(t(as.matrix(xOrd)), border= NA,
        col = c("#00a1d5","#374e55"), ylab="Ancestry Proportion")

# order meta data by Group
datGroup <- dat[order(dat$Group),]

# make matrix of group names
group_mat <- as.matrix(unique(datGroup$Group))

# loop to add lines and text
for(i in 1:length(group_mat)) {
  
  group <- group_mat[i]
  pop <- which(datGroup$Group==group)
  
  # add lines
  lines(x = c(max(pop)-0.4, min(pop)+0.4) * 1.2, y = c(-0.005, -0.005))
  
  # median for middle of sample range
  text(x = median(c(max(pop)-0.5, min(pop)+0.5) * 1.2) - 1, 
       y = -0.15, adj=0, labels = group, 
       srt = 270, xpd = TRUE, cex=1.0)
  
}

K <- as.matrix(read.table("./data/pca_outQ5_ECT_625_chrs_3.admix.3.Q"))
dim(K)

# transverse-> order -> transverse back 
x <- t(as.matrix(K))
xOrd<- x[ ,order(dat$Group)]
xOrd<-t(xOrd)

# plot
barplot(t(as.matrix(xOrd)), border= NA,
        col = c("#00a1d5","#b24745","#374e55"), ylab="Ancestry Proportion")

# order meta data by Group
datGroup <- dat[order(dat$Group),]

# make matrix of group names
group_mat <- as.matrix(unique(datGroup$Group))

# loop to add lines and text
for(i in 1:length(group_mat)) {
  
  group <- group_mat[i]
  pop <- which(datGroup$Group==group)
  
  # add lines
  lines(x = c(max(pop)-0.4, min(pop)+0.4) * 1.2, y = c(-0.005, -0.005))
  
  # median for middle of sample range
  text(x = median(c(max(pop)-0.5, min(pop)+0.5) * 1.2) - 1, 
       y = -0.15, adj=0, labels = group, 
       srt = 270, xpd = TRUE, cex=1.0)
  
}


K <- as.matrix(read.table("./data/pca_outQ5_ECT_625_chrs_4.admix.4.Q"))
dim(K)

# transverse-> order -> transverse back 
x <- t(as.matrix(K))
xOrd<- x[ ,order(dat$Group)]
xOrd<-t(xOrd)

# plot
barplot(t(as.matrix(xOrd)), border= NA,
        col = c("#b24745","#df8f44","#00a1d5","#374e55"), ylab="Ancestry Proportion")

# order meta data by Group
datGroup <- dat[order(dat$Group),]

# make matrix of group names
group_mat <- as.matrix(unique(datGroup$Group))

# loop to add lines and text
for(i in 1:length(group_mat)) {
  
  group <- group_mat[i]
  pop <- which(datGroup$Group==group)
  
  # add lines
  lines(x = c(max(pop)-0.4, min(pop)+0.4) * 1.2, y = c(-0.005, -0.005))
  
  # median for middle of sample range
  text(x = median(c(max(pop)-0.5, min(pop)+0.5) * 1.2) - 1, 
       y = -0.15, adj=0, labels = group, 
       srt = 270, xpd = TRUE, cex=1.0)
  
}

K <- as.matrix(read.table("./data/pca_outQ5_ECT_625_chrs_5.admix.5.Q"))
dim(K)

# transverse-> order -> transverse back 
x <- t(as.matrix(K))
xOrd<- x[ ,order(dat$Group)]
xOrd<-t(xOrd)

# plot
barplot(t(as.matrix(xOrd)), border= NA,
        col = c("#374e55","#df8f44","#b24745","#00a1d5","black"), ylab="Ancestry Proportion")

# order meta data by Group
datGroup <- dat[order(dat$Group),]

# make matrix of group names
group_mat <- as.matrix(unique(datGroup$Group))

# loop to add lines and text
for(i in 1:length(group_mat)) {
  
  group <- group_mat[i]
  pop <- which(datGroup$Group==group)
  
  # add lines
  lines(x = c(max(pop)-0.4, min(pop)+0.4) * 1.2, y = c(-0.005, -0.005))
  
  # median for middle of sample range
  text(x = median(c(max(pop)-0.5, min(pop)+0.5) * 1.2) - 1, 
       y = -0.15, adj=0, labels = group, 
       srt = 270, xpd = TRUE, cex=1.0)
  
}
