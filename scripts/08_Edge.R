library(edgeR)
library(limma)
library(RColorBrewer)
library(mixOmics)
library(HTSFilter)
library(ggplot2)
library(ggrepel)
library(dplyr)


# set the directory from which files are imported
directory <- "/srv/data/my_shared_data_folder/ace2covid/results/counts/"
dir(directory)

# read files
rawCountTable <- read.table(paste0(directory,"reads_count.txt"), header=TRUE,
                            row.names=1, sep = "\t")
sampleInfo <- read.table(paste0(directory,"Metadata.txt"), header=TRUE,
                         row.names=1, sep = "\t")

head(rawCountTable)
nrow(rawCountTable)

rawCountTable <- rawCountTable[,match(rownames(sampleInfo),
                                      colnames(rawCountTable))]
head(rawCountTable)
colnames(rawCountTable)
rownames(sampleInfo)

# create the â€˜conditionâ€™ column
condition = c("african", "african", "african", "african", "african", "african", "african", 
              "african", "african", "african", "african", "african", "african", "african",
              "non-african", "non-african", "non-african", "non-african", "non-african", "non-african",
              "non-african", "non-african", "non-african", "non-african", "non-african", "non-african",
              "non-african", "non-african", "non-african", "non-african", "non-african", "non-african", 
              "non-african", "non-african", "non-african", "non-african", "non-african", "non-african",
              "non-african", "non-african", "non-african", "non-african", "non-african", "non-african",
              "non-african", "non-african", "non-african")

sampleInfo = data.frame(sampleInfo, condition)

#sampleInfo = subset(sampleInfo, select= -condition)
sampleInfo$condition

# create a DGEList data object
dgeFull <- DGEList(rawCountTable, group=sampleInfo$condition)
dgeFull


# add the sample information object in the DGEList data
dgeFull$sampleInfo <- sampleInfo
dgeFull
head(dgeFull$sampleInfo)

# preparing the data object for the analysis 
# select the subset paired-end samples from degFull
dge <- DGEList(dgeFull$counts)
dge$sampleInfo <- dgeFull$sampleInfo[dgeFull$sampleInfo$type=="paired-end",]

# data exploration and quality assessment
# extract pseudo-counts (ie \(\log_2(K+1)\))
pseudoCounts <- log2(dge$counts+1)
head(pseudoCounts)

# histogram for pseudo-counts
hist(pseudoCounts[,"SRR13081202"])

# boxplot for pseudo-counts
boxplot(pseudoCounts, col="cyan")

# MA-plots between samples


# MDS for pseudo-counts (using limma package)
plotMDS(pseudoCounts)

# heatmap for pseudo-counts (using mixOmics package)
sampleDists <- as.matrix(dist(t(pseudoCounts)))
sampleDists

cimColor <- colorRampPalette(rev(brewer.pal(9, "Blues")))(16)
#cim(sampleDists, col=cimColor, symkey=FALSE)

# Differential expression analysis
# remove genes with zero counts for all samples

dge <- DGEList(dge$counts[apply(dge$counts, 1, sum) != 0, ],
               group=sampleInfo$condition)
dge$sampleInfo <- dge$sampleInfo
head(dge$counts)

# estimate the normalization factors
dge <- calcNormFactors(dge, method="TMM")
dge$samples$group <- relevel(dge$samples$group, ref="non-african") #sets non-african as the control group
dge$samples$group

# estimate common and tagwise dispersion
dge <- estimateCommonDisp(dge)
dge <- estimateTagwiseDisp(dge)
dge

# perform an exact test for the difference in expression between the conditions african and non-african
dgeTest <- exactTest(dge)
dgeTest

# Independant filtering
#  remove low expressed genes

filtData <- HTSFilter(dge)$filteredData

dgeTestFilt <- exactTest(filtData)
dgeTestFilt

# Diagnostic plot for multiple testing
# plot a histogram of unadjusted p-values
hist(dgeTest$table[,"PValue"], breaks=50)


# plot a histogram of unadjusted p-values after filtering
hist(dgeTestFilt$table[,"PValue"], breaks=50)


# Inspecting the results
# extract a summary of the differential expression statistics
resNoFilt <- topTags(dgeTest, n=nrow(dgeTest$table))
head(resNoFilt)

#resFilt <- topTags(dgeTestFilt, n=nrow(dgeTest$table))
#head(resFilt)

#rownames(resFilt)

#resFilt$genelabels <- factor(resFilt$gene, levels = c("H2BC5", "CCNF", "H2BC21", "H3C4"))
#resFilt$
# compare the number of differentially expressed genes with and without filtering
# before independent filtering
sum(resNoFilt$table$FDR < 0.05)

# after independent filtering
sum(resFilt$table$FDR < 0.05)

# extract and sort differentially expressed genes
sigDownReg <- resNoFilt$table[resNoFilt$table$FDR<0.05,]
sigDownReg <- sigDownReg[order(sigDownReg$logFC),]
head(sigDownReg)

sigUpReg <- sigDownReg[order(sigDownReg$logFC, decreasing=TRUE),]
head(sigUpReg)

# write the results in csv files
write.csv(sigDownReg, file="sigDownReg.csv")
write.csv(sigUpReg, file="sigUpReg.csv")

# Interpreting the DE analysis results
# create a MA plot with 5% differentially expressed genes

plotSmear(dgeTestFilt,
          de.tags = rownames(resNoFilt$table)[which(resNoFilt$table$FDR<0.05)])

# create a Volcano plot

# create a column with thresholds
significance = -log10(resNoFilt$table$FDR) < 1
resNoFilt$table$significance <- significance
dim(resNoFilt$table)

## setting the values
resNoFilt$table$diffexpressed <- "No"

# if logFC > 0.1 and FDR < 0.05, set as "UP"
resNoFilt$table$diffexpressed[resNoFilt$table$logFC > 0.1 & resNoFilt$table$FDR < 0.05] <- "Up"

# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
resNoFilt$table$diffexpressed[resNoFilt$table$logFC < -1 & resNoFilt$table$FDR < 0.05] <- "Down"

# set different colors
mycolors <- c("cornflowerblue", "black", "firebrick")
names(mycolors) <- c("Down", "No", "Up")

dev.off()


input <- cbind(gene=rownames(resNoFilt$table), resNoFilt$table) #convert the rownames to a column
volc = ggplot(input, aes(logFC, -log10(PValue))) + #volcanoplot with logFC versus pvalue
  geom_point() +
  theme_classic()+
  geom_vline(xintercept=c(-1.0, 0.8), col="orange", linetype="dashed") +
  geom_hline(yintercept=-log10(0.005), col="orange" , linetype="dashed")+
  geom_point(aes(col=diffexpressed)) + #add points colored by significance
  scale_color_manual(values= mycolors) + 
  ggtitle("african vs non-african") #e.g. 'Volcanoplot DESeq2'
volc+geom_text_repel(data=head(input, 1500), aes(label=gene), size = 3) #adding text for the top 20 genes
#ggsave("Volcanoplot.jpeg", device="jpeg") #In case you want to easily save to disk


# transform the normalized counts in log-counts-per-million
y <- cpm(dge, log=TRUE, prior.count = 1)
head(y)

# select 1% differentially expressed genes and produce a heatmap
selY <- y[rownames(resNoFilt$table)[resNoFilt$table$FDR<0.05 & 
                                    abs(resNoFilt$table$logFC)>1.5],]

cimColor <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)[255:1]
cim(selY, margins = c(10, 5))
