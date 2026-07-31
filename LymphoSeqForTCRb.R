library(LymphoSeq)
library(ggplot2)
library(vegan)
library(dplyr)
library(ggpubr)


prefix = "IpiNivo_TCRb"
TCRB.path <- "GBM_IpiNivo_TCR" #put all the .tsv files in one folder
TCRB.list <- readImmunoSeq(path = TCRB.path)

names(TCRB.list) #each data frame listed in the TCRB.list object is named according the ImmunoSEQ file names

#load a metadata file for sample subsetting
TCRB.metafile <- "IpiNivo_TCRb.Metadata.txt"
TCRB.metadata <- read.delim(TCRB.metafile,sep="\t",header=T,row.names=1)
TCRB.metadata <- TCRB.metadata[TCRB.metadata$TCRbID %in% names(TCRB.list),]

####Extracting productive sequences################
###################################################
productive.TRB.aa <- productiveSeq(file.list = TCRB.list, aggregate = "aminoAcid", prevalence = TRUE) 
#Alternatively, you may set aggregate to "nucleotide"
productive.TRB.nt <- productiveSeq(file.list = TCRB.list, aggregate = "nucleotide",prevalence = FALSE)

head(productive.TRB.aa[[1]][,3])
head(productive.TRB.nt[[1]][,3]) 
#no amino acid sequences given means nucleotide sequence is out of frame and does not produce a productively transcribed amino acid sequence.
#If an asterisk (*) appears in the amino acid sequences, this would indicate an early stop codon.


####Create a table of summary statistics##########
##################################################
summaryStat <- clonality(file.list = TCRB.list)
#Both Gini coefficient and clonality are reported on a scale from 0 to 1 where 0 indicates all sequences have the same frequency and 1 indicates the repertoire is dominated by a single sequence.
summaryStat <- summaryStat %>% arrange(factor(samples, levels = names(TCRB.list)))

#we can also include the normalized shannon diversity  index, which is also the 1-clonality
diversity_list_shannon = c()
norm_diversity_list_shannon = c()
diversity_list_simpson = c()
richness_list = c()

for(i in 1:length(TCRB.list)){
  table = TCRB.list[[i]][,c(4)]
  H_shannon = diversity(table,index = "shannon", MARGIN = 1, base = exp(1))
  H_simpson = diversity(table,index = "simpson", MARGIN = 1, base = exp(1))
  
  S = specnumber(table)
  norm_shannon = H_shannon/log(S)
  
  diversity_list_shannon <- c(diversity_list_shannon, H_shannon)
  norm_diversity_list_shannon <- c(norm_diversity_list_shannon,norm_shannon)
  diversity_list_simpson <- c(diversity_list_simpson, H_simpson)
  richness_list <- c(richness_list, S)
}

summaryStat$Diversity_shannon = diversity_list_shannon
summaryStat$Diversity_shannon_norm = norm_diversity_list_shannon
summaryStat$Diversity_simpson = diversity_list_simpson
summaryStat$Richness = richness_list

summaryStat$RNAseqID = TCRB.metadata$RNAseqID

tOut = paste(prefix,"SummaryStat","txt",sep=".")
write.table(summaryStat, file=tOut, sep="\t", row.names=TRUE, col.names=TRUE)


####Visualizing repertoire diversity/clonality#############
#################################################
#plot boxplot of selected summary stats ##########################################
summaryStat$Treatment = TCRB.metadata$Pre.surgical.Arm
rownames(summaryStat) = summaryStat$samples

summaryStat$Treatment = factor(summaryStat$Treatment, levels=c(3,2,1))

varlist = c("clonality",
            "giniCoefficient",
            "Diversity_shannon",
            "Diversity_shannon_norm",
            "Diversity_simpson",
            "Richness")


colors = c("black", "black", "black")
colors_fill = c("grey", "grey", "grey")
ncol = 1  
nrow1 = length(varlist)/ncol
width1 = ncol*400*ncol
height1 = nrow1*400

outfile = paste(prefix,"Boxplot.SummaryStat.pdf",sep=".") 
pdf(outfile, width=width1/72, height=height1/72)

panels=list()
l=length(colnames(expr_gene))

for(i in 1:length(varlist)){
  var = varlist[i]
  data = summaryStat[,c(var, "Treatment")]
  colnames(data)=c("normExpr","condition")
  
  
  y = ggplot(data=data, aes(x=condition,y=normExpr,color=condition, fill=condition)) +
    #geom_violin(trim = FALSE)+
    geom_boxplot(fill="white", width=0.8, colour="black")+
    geom_dotplot(binaxis='y', stackdir='center',dotsize = 1)+
    #geom_jitter(width=0.3, alpha=0.3)+
    scale_fill_manual(values=colors_fill)+
    scale_color_manual(values=colors)+
    theme(
      axis.text.x = element_text(size=18),
      axis.text.y = element_text(size=18),
      axis.ticks.x = element_blank(),
      axis.title = element_blank(),
      legend.position = "none",
      plot.title=element_text(face = "italic",size=16))+
    ggtitle(var)
  
  panels[[i]]=y
  
}
ggarrange(plotlist=panels,ncol=ncol,nrow=nrow1)
dev.off()


#calculate p-values #################
tk <- kruskal.test(clonality ~ Treatment, data = summaryStat)
print(tk$p.value)

tk <- kruskal.test(giniCoefficient ~ Treatment, data = summaryStat)
print(tk$p.value)

tk <- kruskal.test(Diversity_shannon ~ Treatment, data = summaryStat)
print(tk$p.value)

tk <- kruskal.test(Diversity_simpson ~ Treatment, data = summaryStat)
print(tk$p.value)

tk <- kruskal.test(Richness ~ Treatment, data = summaryStat)
print(tk$p.value)


for(i in 1:length(varlist)){
  var = varlist[i]
  
  sub1 <- subset(summaryStat,Treatment==1)
  sub2 <- subset(summaryStat,Treatment==2)
  t <- t.test(sub1[,var],sub2[,var],paired=F)
  Delta = median(sub2[,var])-median(sub1[,var])
  p = t$p.value
  
  print(var)
  print("2vs1")
  print(p)
  print(Delta)
  
  sub1 <- subset(summaryStat,Treatment==2)
  sub2 <- subset(summaryStat,Treatment==3)
  t <- wilcox.test(sub1[,var],sub2[,var],paired=F)
  Delta = median(sub2[,var])-median(sub1[,var])
  p = t$p.value
  
  print(var)
  print("3vs2")
  print(p)
  print(Delta)
  
  sub1 <- subset(summaryStat,Treatment==1)
  sub2 <- subset(summaryStat,Treatment==3)
  t <- wilcox.test(sub1[,var],sub2[,var],paired=F)
  Delta = median(sub2[,var])-median(sub1[,var])
  p = t$p.value
  
  print(var)
  print("3vs1")
  print(p)
  print(Delta)

}


##find the top clones for each sample ################################################
######################################################################################
#The uniqueSeqs function creates a data frame of all unique, productive sequences and reports 
#the total count in all samples.
unique.seqs <- uniqueSeqs(productive.aa = productive.TRB.aa)
head(unique.seqs)
sequence.matrix <- seqMatrix(productive.aa = productive.TRB.aa, sequences = unique.seqs$aminoAcid)
head(sequence.matrix)
top.freq <- topFreq(productive.aa = productive.TRB.aa, percent = 0.01)
head(top.freq)


#look for sequences that is shared and expanded in ipi+nivo compared to other arms
## 1) Precompute sample columns and arm indices ONCE
sequence.matrix_sub <- sequence.matrix[sequence.matrix$aminoAcid %in% top.freq$aminoAcid,]
sample_cols <- 3:ncol(sequence.matrix_sub)
freq_mat <- as.matrix(sequence.matrix_sub[, sample_cols, drop = FALSE])

arm1_ids   <- TCRB.metadata$TCRbID[TCRB.metadata$Pre.surgical.Arm == 1]
arm1_idx   <- which(colnames(freq_mat) %in% arm1_ids)

arm2_ids   <- TCRB.metadata$TCRbID[TCRB.metadata$Pre.surgical.Arm == 2]
arm2_idx   <- which(colnames(freq_mat) %in% arm2_ids)

arm3_ids   <- TCRB.metadata$TCRbID[TCRB.metadata$Pre.surgical.Arm == 3]
arm3_idx   <- which(colnames(freq_mat) %in% arm3_ids)

other_idx <- arm3_idx

## (optional sanity check)
# stopifnot(length(arm1_idx) > 0, length(other_idx) > 0)

## 2) Vectorized-ish Wilcoxon over rows (apply is much faster than reshape+merge per loop)
pval <- apply(freq_mat, 1, function(v) {
  x <- v[arm2_idx]
  y <- v[other_idx]
  # Use asymptotic p-value for speed
  wilcox.test(x, y, alternative = "greater", exact = FALSE, correct = FALSE)$p.value
})

## 3) Assemble results without growing vectors in a loop
antigen <- top.freq$antigen[match(sequence.matrix_sub$aminoAcid, top.freq$aminoAcid)]

out <- data.frame(
  aminoAcid = sequence.matrix_sub$aminoAcid,
  antigen = antigen,
  numberSamples = sequence.matrix_sub$numberSamples,
  pval = pval,
  stringsAsFactors = FALSE
)

out$Padj <- p.adjust(out$pval,method = "fdr")

tOut = paste(prefix,"ExpandedTCR","Arm2vs3","txt",sep=".")
write.table(out, file=tOut, sep="\t", row.names=FALSE, col.names=TRUE)

sequence_list <- c("CASSLGSYEQYF")
idx_list <- which(sequence.matrix$aminoAcid %in% sequence_list)

seq_freq_tran <- as.data.frame(t(sequence.matrix[idx_list,3:ncol(sequence.matrix)]))
colnames(seq_freq_tran) <- sequence_list
seq_freq_tran$TCRbID <- rownames(seq_freq_tran)

merged_df <- merge(
  seq_freq_tran,
  TCRB.metadata,
  by = "TCRbID",
  all.x = TRUE
)


merged_df$Pre.surgical.Arm = factor(merged_df$Pre.surgical.Arm, levels=c(3,2,1))

colors = c("black", "black", "black")
colors_fill = c("grey", "grey", "grey")
ncol = 1  
nrow1 = length(sequence_list)/ncol
width1 = ncol*400*ncol
height1 = nrow1*400

outfile = paste(prefix,"Boxplot","ExpandedTCR","Arm1vs3","pdf",sep=".") 
pdf(outfile, width=width1/72, height=height1/72)

panels=list()
#l=length(colnames(expr_gene))

for(i in 1:length(sequence_list)){
  var = sequence_list[i]
  data = merged_df[,c(var, "Pre.surgical.Arm","TCRbID")]
  colnames(data)=c("normExpr","condition","TCRbID")
  
  
  y = ggplot(data=data, aes(x=condition,y=normExpr,color=condition, fill=condition)) +
    #geom_violin(trim = FALSE)+
    geom_boxplot(fill="white", width=0.8, colour="black", outlier.size=0, outlier.alpha = 0)+
    geom_dotplot(binaxis='y', stackdir='center',dotsize = 1)+
    #geom_jitter(width=0.3, alpha=0.3)+
    scale_fill_manual(values=colors_fill)+
    scale_color_manual(values=colors)+
    theme(
      axis.text.x = element_text(size=18),
      axis.text.y = element_text(size=18),
      axis.ticks.x = element_blank(),
      axis.title = element_blank(),
      legend.position = "none",
      plot.title=element_text(face = "italic",size=16))+
    ggtitle(var)
  
  panels[[i]]=y
  
}
ggarrange(plotlist=panels,ncol=ncol,nrow=nrow1)
dev.off()



#One very useful thing to do is merge the output of seqMatrix and topFreq to get 
#the summary table of each sequence's frequency in each sample
top.freq <- topFreq(productive.aa = productive.TRB.aa, percent = 0)
top.freq.matrix <- merge(top.freq, sequence.matrix, by = "aminoAcid")
head(top.freq.matrix)

tOut = paste(prefix,"FreqMatrix","AA","txt",sep=".")
write.table(top.freq.matrix, file=tOut, sep="\t", row.names=FALSE, col.names=TRUE)


####Comparing samples##################################
#######################################################

## check the similarity of the top 50 sequences for each sample
productive.TRB.aa.top <- lapply(
  productive.TRB.aa,
  function(df) {
    df[order(df$frequencyCount, decreasing = TRUE), ][1:50, ]
  }
)

similarity.matrix <- similarityMatrix(productive.seqs = productive.TRB.aa.top) # similarity scores involves weighting each shared sequence in the two distributions by the geometric mean of the frequency of each sequence in the two distributions.

diag(similarity.matrix) <- 0

## plot the similarity matrix for each treatment arm
treatment = 1
samples = TCRB.metadata[TCRB.metadata$Pre.surgical.Arm == treatment,]$TCRbID


outfile = paste(prefix,"SampleSimilarity.AA","Arm",treatment,"top50","pdf",sep=".")
pdf(outfile,width = 850/72, height = 800/72)
pairwisePlot(matrix = similarity.matrix[rev(samples),rev(samples)])+ 
  ggplot2::scale_fill_gradient(low = "#deebf7", high = "#3182bd") + 
  ggplot2::labs(fill = "Similarity score")
dev.off()




