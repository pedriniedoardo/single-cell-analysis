
setwd("~/R/SingleCellqPCR")

dat <- read.csv("data mock vs 48h.csv", header=TRUE, skip=11)

#96 cell expressions S01-S96 of 48 genes assys A01-A96, (2 technical replicates Ayy, Azz per gene for each cell Sxx),
#of single(nCells=1) and 100 cells(nCells=100)

#mock_xxxx = Control, 48h_xxx = treated, _low = low MDR (staining), _high = high MDR expression each of the 96 of the single cells has 48 primers with two technical repeats, which Ct values going to be averaged

colnames(dat) <- c("Cell_ID","Treatment_MDR","Sample_Type","nCells","Gene" ,"Exp_Type", 
                   "Ct", "Calibrated", "Ct_Quality","Ct_Outcome","Treshold",         
                   "In_Range","Out_Range","Peak_Ratio",
                   "Comments")
#reverse the order of Cell_IDs easier keeping track 
dat <- dat[order(dat$Cell_ID, dat$Gene),]


LOD <- 27
#explore and add coulms to data frame with package dplyr
library(dplyr)

#Ct < LOD and "Pass", stays the same

#Ct < LOD and "Fail", convert to Ct=-1 an Log2Ex=mean(gene exp of treatment group)

#Ct > LOD  and "Pass" or "Fail", convert to Ct = LOD an Log2Ex = 0

############################### extract the data for the analysis ################ Human RNA(all transcripts) and empty(no transcripts) are the controls, used to  figure out limit of detection(?) remove from df for analysis
df <- dat[!dat$Treatment_MDR %in% c("empty", "All Human RNA"),]
#drop unsued Treatments factor levels
df$Treatment_MDR <- droplevels(df$Treatment_MDR)
#remove their Cell Id levels too
df$Cell_ID  <- droplevels(df$Cell_ID)

########################## Pre-processing data ################################

#Handling Missing (reaction that did not give rise to any Cq values) and Offscale data (reaction that give rise to Cq value too high to be truste) acording to:

# ---------------------- Fluidigm manual ----------------------------


#Chapter 3, p.33, see Examples 1, define data as
#missing due to reaction failure, failed curves or melt curves should be restored based on replicate information
#extreem values such as 999 its not a missing data point, but means there is no expression (?), with either Pass or Fail Ct reaction outcome

#extract the variables used for the data analysis 
#keep <- c("Cell_ID", "Treatment_MDR","nCells","Gene","Ct")
#df <- df[keep]

# ------------------------- qPCR work flow for single cell ---------------
#sec 2.6 and fig.3

#################################### remodeling the data frame #################

#extract cell number (the first 3 symbols) from the Cell_ID,rename the Cell_ID labels and drop the unsused assay number
df$Cell_ID <- substr(df$Cell_ID, start=1, stop=3)

#Average genes Log2Ex values of technical replicates taking their mean (which corespons to the geometric mean of the linear expression values)
df <- ddply(df, .(Cell_ID, Gene), mutate, AvgLog= mean(Log2Ex, na.rm="TRUE"))

df <- df %>% distinct(Cell_ID, Gene)

#Gene expressions in linear scale 2^Log2Ex!!!!!!!!!!!!!!!!!!!!!!!

data <- select(df, Cell_ID, Treatment_MDR, nCells, Gene, AvgLog)

hc <- filter(data, nCells == 100)
hc$nCells <- NULL
hc$Cell_ID <- as.factor(hc$Cell_ID)

sc <- filter(data, nCells == 1)
sc$nCells <- NULL
sc$Cell_ID <- as.factor(sc$Cell_ID)

colnames(sc) <- c("Id","Treatment", "Gene", "Log2Ex")
#make a new colum with type of sample: controls (mock_low and mock_high) and treated (48h_high, 48_low)
sc$Type <- sc$Treatment
hc$Type <- hc$Treatment
#combine levels
levels(sc$Type) <- c("Treated","Treated","Control","Control")
levels(hc$Type) <- c("Treated","Treated","Control","Control")

#Histograms of Expression levels of 80 cells expresing gene in logarithmic and linear scale
#generate plot per gene
gene <- "ABCC4"
ggplot(filter(sc, Gene == "ABCC4"), aes(Log2Ex, color=Treatment, fill=Treatment)) + geom_histogram(bins=40) + 
geom_freqpoly(binwidth=0.5)  +
ylab("No of cells") + labs(title=gene)

ggplot(filter(sc,Gene == gene),aes(Log2Ex, ..density..)) + geom_density(aes(fill=Type, position="stack"), alpha=0.5) + 
 ylab("No of cells")

ggplot(filter(sc, Gene == "ABCC4"), aes(Log2Ex, fill=Treatment)) + geom_histogram(bins=40) + facet_wrap(~Treatment) +
#geom_freqpoly(binwidth=0.5)  +
ylab("No of cells") + labs(title=gene)

library(dplyr)
#generate plots for 8 genes at a time
genes <- levels(sc$Gene)[1:8]
data.set <- filter(sc,Gene %in% genes)

ggplot(data.set,aes(Log2Ex,fill=Type)) +
geom_histogram(bins=40) + facet_wrap(~ Gene, scales="free") + 
theme(strip.background=element_rect(fill="white")) +
ylab("No of cells") + theme_bw()

ggplot(data.set, aes(Gene, Log2Ex, fill=Type)) + geom_violin(scale = "width", trim=TRUE, adjust=0.35,position=position_dodge(width=0.65))
#theme(legend.position="none") +
  
ggplot(data.set, aes(Gene, Log2Ex, fill=Treatment)) + geom_violin(scale = "width", trim=TRUE, adjust=0.35,position=position_dodge(width=0.65), alpha=0.75)  
#Use a smaller bandwidth for closer density fit (default is 1)
  

#transform data for PCA: variables(Genes) as columns, samples(Cells) as rows 
require(reshape2)
scge <- dcast(sc, Id + Treatment + Type ~ Gene, value.var="Log2Ex")

#To normalized remove Genes(columns) with constant expressions or NA values

zero.var.genes <- names(which(apply(hcge[ ,-(1:3)], 2, var, na.rm=TRUE) == 0))
na.genes <- names(which(is.na(apply(hcge[ ,-(1:3)], 2, var, na.rm=TRUE))))

remove.cols <- c(zero.var.genes, na.genes)
hcge <- hcge[ , !names(hcge) %in% remove.cols]

hcge <- hcge [,sapply(hcge, function(c) var(c, na.rm=TRUE) !=0)]

# ----------------------- PCA ---------------------
#princomp gives the scores
pca<- prcomp(scge[,-(1:3)], center=TRUE, scale.=TRUE)
names(pca)
head(pca$sdev)
#use summary to extract the eigenvalues and variances from an object of class prcomp
summary(pca)

# Variance of the principal components scree plots
screeplot(pca.scge, npcs = nrow(pca$rotation), type = "l", pch=19, main="Scree Plot PCA for single cells"); abline(h=1, col="red", lwd=3)

eig <- (pca$sdev)^2
#variance in percentage
variance <- eig*100/sum(eig)
#cumulative variances with cumulative sums
cumvar <- cumsum(variance)


eig.var.df <- data.frame(eig = eig, variance = variance, cumvariance = cumvar)
head(eig.var.df)

#Scree plot using base graphics
barplot(eig.var.df[, 2], names.arg=1:nrow(eig.var.df), 
        main = "Variances",
        xlab = "Principal Components",
        ylab = "Percentage of variances",
        col ="steelblue")
# Add connected line segments to the plot
lines(x = 1:nrow(eig.var.df), 
      eig.var.df[, 2], 
      type="b", pch=19, col = "red")


#scree plots using package factoextra eigenvalues
library(factoextra)
eig.val <- get_eigenvalue(pca)
head(eig.val)
fviz_screeplot(pca, ncp=11)

# A good dimension reduction is achieved when the the first few PCs account for a large proportion of the variability (80-90%)
eig.var.df

# ------------------- pca plots --------------------------------

pca.data <- cbind(scge[,1:3], pca$x)

pca1.2 <- ggplot(data = pca.data, aes(x=PC1, y=PC2))
pca1.3 <- ggplot(data = pca.data, aes(x=PC1, y=PC3))
pca1.4 <- ggplot(data = pca.data, aes(x=PC1, y=PC4))

pca2.3 <- ggplot(data = pca.data, aes(x=PC2, y=PC3))
pca2.4 <- ggplot(data = pca.data, aes(x=PC2, y=PC4))

pca2.3 <- ggplot(data = pca.data, aes(x=PC3, y=PC4))

#summary(pca1.2)
pca1.2 + geom_point(aes(colour=Treatment, shape=Type), size=3, alpha = 0.75) +
  labs(title="The first two PCA projections explaining only 77% of cell variability") + 
  scale_color_manual(values = c("orangered3","orange2","deepskyblue4","deepskyblue3")) +
  scale_shape_manual(values = c(17,16)) +
  theme_bw() +
  #geom_vline(xintercept=0, col="darkred") +
  #geom_hline(yintercept=0, col="darkred") 

  
library(rgl)
plot3d(pca$x[,1:3], col=pca$Treatment)  

ed <- dist(pca.data[,4:45], method="euclidean")
cols.treat <- c(rep("red",20), rep("orange", 20), rep("blue", 20),  rep("lightblue", 20)) 

library("MASS")
sam.2d <-  sammon(ed, trace=FALSE)
plot(sam.2d$points, col=cols.treat, pch=18, cex=1.5,
     xlab="Dimension 1", ylab="Dimension 2",
     main="Non-metric MDS based on Euclidean Distance bw samples in PCA space") 

#to see in interactive 3D use rgl package
library("rgl")
sam.3d <-  sammon(ed, k=3, trace=FALSE)

plot3d(sam.3d$points, col=cols.treat, type='s', size=2, 
       xlab="Dimension 1", ylab="Dimension 2", zlab="Dimension 3",
       main="3D Non-metric MSD based of PCA coordinates")

text3d(sam.3d$points, col=cols.treat, cex=1.25, adj = 1.5) #text=pca.data$Treatment
#rgl.snapshot("persp3dd.png","png")

#save the 3d plot as pdf file and coppy into PPT with snapshot
rgl.postscript("pca3d.pdf","pdf", drawText=TRUE) 

# ---------------- Gene Loadings ------------------------------
gene.loadings <- data.frame(Genes, pca$rotation, stringsAsFactors = FALSE)
gl <- ggplot(data=gene.loadings, aes(x=PC1, y=PC2, label=Genes))

gl + geom_point(aes(colour=factor(Genes)), size=2) +
geom_text(aes(colour=factor(Genes)), vjust=1.7, hjust=-0.1, size=2.5) +
theme_bw() + 
geom_vline(xintercept=0, col="red") +
geom_hline(yintercept=0, col="red") +
labs(title="Loading of the 48 genes") +  theme_bw() 

############### bar charts with PC faceting of Gene Loadins
gene.rotations <- data.frame(Genes=rownames(pca$rotation), pca$rotation)
m.gene.rotations <- melt(gene.rotations)


