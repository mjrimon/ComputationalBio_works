#!/home/miniconda3/bin/Rscript


library("DESeq2")
#library("RColorBrewer")
#library("cluster")
#library("ggplot2")

# Create an object with the directory containing HTseq counts:
directory <- "~/results/genome_analysis/DEA"


# We have to ensure same object in all vector, therefore we write down de names
sampleFiles <- c( "hightemp01.count", "hightemp02.count","normal01.count","normal02.count")
# create a vector of sample names. Ensure these are in the same order as the "sampleFiles" object!
sampleNames <- c( "hightemp01", "hightemp02","normal01","normal02")
# next line get filenames from sampleNames
#sampleFiles <- sapply(sampleNames, paste0, '.count', USE.NAMES = FALSE)
# create a vector of conditions. again, mind that they are ordered correctly!
sampleCondition <- c("hightemp","hightemp","normal", "normal")
# create a vector of conditions. again, mind that they are ordered correctly!
replicate <- c("Rep1","Rep2","Rep1","Rep2")
# now create a data frame from these three vectors.
sampleTable <- data.frame(
    sampleName = sampleNames,
    fileName = sampleFiles,
    condition = sampleCondition,
    replicate = replicate)

#sampleTable

## Make DESeq2 object from counts and metadata
ddsHTSeq <- DESeqDataSetFromHTSeqCount( sampleTable = sampleTable,
                                        directory = directory,
                                        design = ~condition)
#ddsHTSeq
#ddsHTSeq$condition
#specifying the reference level:
ddsHTSeq$condition <- relevel(ddsHTSeq$condition, ref = "normal")
#ddsHTSeq$condition

#Here we perform a minimal pre-filtering to keep only rows that have at least 10 reads total.
# what does expression look like across genes?
# sum counts for each gene across samples
sumcounts <- rowSums(counts(ddsHTSeq))
# take the log
logsumcounts <- log(sumcounts,base=10)
# plot a histogram of the log scaled counts
pdf(paste0(directory, "/hist_log_count.pdf"))
plot_without_filter <- hist(logsumcounts,breaks=100, main="Histogram of the log scaled counts")
dev.off()
# get genes with summed counts greater than 10; remove lowly expressed genes
keep <- sumcounts > 10
# keep only the genes for which the vector "keep" is TRUE
ddsHTSeq_filter <- ddsHTSeq[keep,]
sumcounts <- rowSums(counts(ddsHTSeq_filter))
# take the log
logsumcounts <- log(sumcounts,base=10)
# plot a histogram of the log scaled counts
pdf(paste0(directory, "/hist_log_count_with_filter.pdf"))
plot_with_filter <- hist(logsumcounts,breaks=100,main="Histogram of the log scaled counts \n keep only the genes with summed counts greater than 10" )
dev.off()


dds <- DESeq(ddsHTSeq_filter)
pdf(paste0(directory, "/dispEsts.pdf"))
plotDispEsts(dds, main="Dispersion plot")
dev.off()

# get results table
# Get basic statistics about the number of significant genes
# treatedvscontrol
#res <- results(dds, pAdjustMethod="BH")
#res_DF <- as.data.frame(res)
resultsNames(dds)
#summary(res)

# check out the first few lines
#head(res)
#mcols(res, use.names = T)
#resultsNames(dds)

##Create normalized read counts
#normalized_counts <- counts(dds, normalized=TRUE)
  #head(normalized_counts)
#normalized_counts_mad <- apply(normalized_counts, 1, mad)
#normalized_counts <- normalized_counts[order(normalized_counts_mad, decreasing=T), ]
  #as.data.frame(normalized_counts)


#pdf("hist_norm_count.pdf")
#plot_with_filter <- hist(normalized_counts,main="Histogram of the normalized counts \n from genes with summed counts greater than 10" )
#dev.off()



#DESeq get results table
Res_A_X_total <- results(dds, name="condition_hightemp_vs_normal", pAdjustMethod="BH")
summary(Res_A_X_total)
Res_A_X_total <- Res_A_X_total[order(Res_A_X_total$padj),]
Res_A_X_total <- data.frame(Res_A_X_total)
#head(Res_A_X_total)
Res_A_X_total <- tibble::rownames_to_column(Res_A_X_total, var = "ensembl_gene_id")
#nrow(Res_A_X_total)

#head(Res_A_X_total)

#rownames(Res_A_X_total)
#Export output files
# In then questions we need to know padj < 0.01
Res_A_X_total$sig <- ifelse(Res_A_X_total$padj < 0.01, "yes", "no")
Res_A_X_tota_0.01 <-subset(Res_A_X_total, padj < 0.01)
#write.table(normalized_counts,"Condition_hightemp_vs_normal_normcounts.csv")
write.table(Res_A_X_total, paste0(directory, "/Condition_hightemp_vs_normal_total.csv"), sep=";", dec=",")
#write.table(Res_A_X_tota_0.05, "Condition_hightemp_vs_normal_result_padj_0.05.csv")
write.table(Res_A_X_tota_0.01, paste0(directory, "/Condition_hightemp_vs_normal_result_padj_0.01.csv"), sep=";", dec=",")


# Annotate significant genes

genomeFile <- '~/final_project/genome.gff'
gff <- read.delim(genomeFile, header=F, comment.char="#") 

colnames(gff) = c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attributes")

# Function to get named attributs from gff attributes
getAttributeField <- function (x, field, attrsep = ";") {
  s = strsplit(x, split = attrsep, fixed = TRUE)
  sapply(s, function(atts) {
    a = strsplit(atts, split = "=", fixed = TRUE)
    m = match(field, sapply(a, "[", 1))
    if (!is.na(m)) {
      rv = a[[m]][2]
    }
    else {
      rv = as.character(NA)
    }
    return(rv)
  })
}


gff$ensembl_gene_id <- getAttributeField(gff$attributes, "ID")
gff$Name <- getAttributeField(gff$attributes, "Name")
gff$Product <- getAttributeField(gff$attributes, "product")

# Only need ensembl_gene_id, name and product from gff to annotate Res_A_X
gff <- gff[c("ensembl_gene_id", "Name", "Product")]

# join columns by ensembl_gene_id
Ann_Res_A_X_tota_0.01 <- merge(Res_A_X_tota_0.01, gff, by='ensembl_gene_id')

Ann_Res_A_X_tota_0.01 <- Ann_Res_A_X_tota_0.01[order(Ann_Res_A_X_tota_0.01$padj, -Ann_Res_A_X_tota_0.01$baseMean),]
# select only asked columns
Ann_Res_A_X_tota_0.01 <- Ann_Res_A_X_tota_0.01[c("ensembl_gene_id","Product","Name","pvalue","padj","log2FoldChange")]

write.table(Ann_Res_A_X_tota_0.01, paste0(directory, "/Annotated_Condition_hightemp_vs_normal_result_padj_0.01.csv"), sep=";", dec=",")

Ann_Res_A_X_tota_0.01



### PlotMA {.tabset .tabset-fade .tabset-pills}

#In DESeq2, the function plotMA shows the log2 fold changes attributable to a given variable over the mean of normalized counts for all the samples in the DESeqDataSet. Points will be colored red if the adjusted p value is less than 0.05. Points which fall out of the window are plotted as open triangles pointing either up or down.
pdf(paste0(directory, "/plotMA.pdf"))

plotMA(dds , alpha = 0.05, main=("Standard Bland-Altman- MA Plot"))

dev.off()

