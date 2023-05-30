#GSE104075
#Differential gene-expression analysis using DEseq2
#authors: Chengqi Wang, J. Oberstaller

# Load count data from featureCounts
input <- "vehicle_drug_feature_counts.txt"
df <- read.table(input,header = T, sep = '\t', row.names = 1)
df <- df[,6:9]

# Check the distribution of counts per sample
par(mfrow = c(2,2), mar = c(4,4,1,1))
print( apply(df, 2, function(x) quantile(as.numeric(x))) )
hist.plots <- apply(df, 2, function(x) {hist(log2(x), breaks = 100)})

## -------------------------------------------------------------------------------------------------------
#Low-expressed genes can make for noisy data. Filtering out the lowest-expressed genes (bottom 25%),which corresponds to less than ~16 (or 2^4) reads. [16-based on 25th percentile values for all 4 samples]
idExp <- apply(df, 1, function(x) any(x > 16))
dfExp <- df[idExp, ]
print(dim(dfExp))

## -------------------------------------------------------------------------------------------------------
#Filtered data is ready for normalization and expression-analysis with DESeq2.
#Building DEseq class
library('DESeq2')
table2 <- data.frame(name = colnames(dfExp),
                     condition = c('control','control', 'treatment', 'treatment'))
dds <- DESeqDataSetFromMatrix(dfExp, colData=table2, design= ~ condition)
dds <- DESeq(dds)

## -------------------------------------------------------------------------------------------------------
#normalized reads count
norCounts <- counts(dds, normalized=TRUE)

#DEseq results
res <- results(dds)

#Getting differentially expressed genes and taking a look at our results object
res

# res is a dataframe object, so we can check out metadata for what the columns mean
mcols(res, use.names=TRUE)

## -------------------------------------------------------------------------------------------------------
#extracting genes with adjusted P-value < 0.01
resSig <- res[ which(res$padj < 0.01), ]
dim(resSig)

## -------------------------------------------------------------------------------------------------------
# sort by log2FoldChange to get the significant genes with the strongest down-regulation
head( resSig[ order( resSig$log2FoldChange ), ] ) #This allows you to quickly inspect the genes or features with the smallest fold changes in expression among the significantly differentially expressed ones.
# and strongest up-regulation
tail( resSig[ order( resSig$log2FoldChange ), ] )


## -------------------------------------------------------------------------------------------------------
resSig.sorted.df <- as.data.frame(resSig[ order( resSig$log2FoldChange ), ]) #contains the sorted significant results from resSig. This allows you to perform various data frame operations or export the sorted results as a separate data frame for further analysis or visualization.

write.table(resSig.sorted.df, file = "resSig.sorted.tab.txt",sep = "\t", row.names = TRUE,quote = FALSE)

##heatmap plot
##plot heatmap of normalized read-counts
library("RColorBrewer")
library("gplots")
sigNorData <- norCounts[which(res$padj < 0.01),]
hmcol <-  colorRampPalette(brewer.pal(9, "GnBu"))(100)

## expression heatmap
heatMapDF_nor <- t( apply(sigNorData, 1, function(x){(x-mean(x))/sd(x)}) )
colnames(heatMapDF_nor) <- c('control1','control2','treat1'  , 'treat2')

pdf('heatmap.pdf', height = 10, width = 10)
heatmap.2(heatMapDF_nor, col = hmcol, trace="none", margin=c(10, 10),labRow = F)
dev.off()


## -------------------------------------------------------------------------------------------------------
#Volcano plots
# Generate a data frame with a column assigning color to significant deferentially expressed genes
res_plot      <- data.frame( res )
res_plot$col  <- 'gray40'

# setting cutoffs for significantly up (red) and down (blue) genes; everything else is gray
res_plot$col[res_plot$log2FoldChange > 1 & res_plot$padj < 0.01] <- 'red'
res_plot$col[res_plot$log2FoldChange < -1 & res_plot$padj < 0.01] <- 'cornflowerblue'

pdf('volcano.pdf', height = 10, width = 10)
plot(res_plot$log2FoldChange,
     -log10(res_plot$padj),
     col = res_plot$col, pch = 19, xlab = 'log2(fold change)',
     ylab = '-log10(p-adj)'
)
dev.off()


## -------------------------------------------------------------------------------------------------------
# Functional enrichment analysis

## load GO gaf file (GO annotation)
geneGO <- read.delim2('PlasmoDB-48_Pfalciparum3D7_GO.gaf.txt', header = F, sep = '\t')
geneGO <- geneGO[,c(2,5)]
library(GO.db)
go_db <- Term(GOTERM)
go_On <- Ontology(GOTERM)
go_inf<- data.frame(ontology = go_On,
                    term     = go_db)
BiocManager::install("GOfuncR")

packages <- c("BiocManager","gplots","knitr","RColorBrewer")
if (length(setdiff(packages, rownames(installed.packages()))) > 0) {
  install.packages(setdiff(packages, rownames(installed.packages())),repos='https://cloud.r-project.org') 
}
# clean up leftover variables from above
rm("packages", "bioconductor.packages")


## -------------------------------------------------------------------------------------------------------
## build function for go enrichment - the go_enrichment function takes gene IDs for the query set and background set, along with gene-GO term annotations. It performs GO enrichment analysis by calculating counts, conducting Fisher's exact test, and returning a data frame with the results, including adjusted p-values.
go_enrichment <- function(geneID_query,
                          geneID_background,
                          geneAnnotation) 
{
  geneQueryInf      <- geneAnnotation[geneAnnotation[,1] %in% geneID_query,]
  geneBackgroundInf <- geneAnnotation[geneAnnotation[,1] %in% geneID_background,]
  queryGOterm <- unique(geneQueryInf[,2]) #Determine unique query GO terms:
  
  goResult <- c() #Perform GO enrichment analysis
  aa <- sapply(queryGOterm, function(id) #initialize an empty vector goResult and iterate over each query GO term using sapply.
    #Calculate counts and perform Fisher's exact test - These lines calculate the number of genes in the query set and background set associated with each query GO term. Then, it calculates the number of genes not associated with the query GO term in both sets. Finally, it performs Fisher's exact test to assess the enrichment of the query GO term.
  {
    numQuery        <- length( which(geneQueryInf[,2] == id) )
    numBackground   <- length( which(geneBackgroundInf[,2] == id))
    
    numQuery_no     <- length(geneID_query) - numQuery
    numBackground_no<- length(geneID_background) - numBackground
    
    #print(c(numBackground, numBackground_no))
    fishTest <- fisher.test(rbind( c(numQuery, numQuery_no),
                                   c(numBackground, numBackground_no) ),
                            alternative = 'greater')
    #Store the results for each query GO term in the goResult matrix. The rbind function appends the results for each GO term as a new row in goResult
    infReturn <- c(numQuery,
                   numQuery_no,
                   numBackground,
                   numBackground_no, fishTest$p.value)
    goResult <<- rbind(goResult, infReturn)
  })
  #Assign row names and column names to the goResult matrix, convert it to a data frame, and calculate adjusted p-values (padj) using the false discovery rate (FDR) method. Finally, the function returns the goResult data frame.
  rownames(goResult) <- queryGOterm
  colnames(goResult) <- c('#QueryWithGOterm',
                          '#QueryWithoutGOterm',
                          '#BackgroundWithGOterm',
                          '#BackgroundWithoutGOterm', 'pvalue')
  goResult <- data.frame(goResult)
  goResult$padj <- p.adjust(goResult$pvalue, method = 'fdr')
  return(goResult)
}


## -------------------------------------------------------------------------------------------------------
goEnrichment <- go_enrichment(rownames(resSig)[resSig$log2FoldChange > 1], #select the row names (gene IDs) from the resSig object where the log2FoldChange value is greater than 1. 
                              rownames(df), #represents the background set of gene IDs, typically including all genes analyzed in the study.
                              geneGO) #the gene-GO annotation data
idMatch <- match(rownames(goEnrichment), rownames(go_inf))
goEnrichment <- data.frame(goEnrichment,
                           go_inf[idMatch, ])
goEnrichment$term[goEnrichment$padj < 0.1] #provide a subset of enriched GO terms that are considered statistically significant based on the chosen threshold



## -------------------------------------------------------------------------------------------------------
### bar plot for go enrichment

goSig   <- goEnrichment[goEnrichment$padj < 0.01,]
goSig$expection <- goSig$X.BackgroundWithGOterm/(goSig$X.BackgroundWithGOterm +
                                                   goSig$X.BackgroundWithoutGOterm) * (goSig$X.QueryWithGOterm + 
                                                                                         goSig$X.QueryWithoutGOterm)

goSig <- goSig[order(goSig$X.QueryWithGOterm), ]
goSigDraw <- t( goSig[,c(1,9)] )
colnames(goSigDraw) <- goSig$term

pdf('GO_barplot.pdf', height = 10, width = 10)
par(mar = c(4,20,1,1))
barplot(goSigDraw, horiz = T, las = 1)
dev.off()

