# GSE104075 Differential Gene Expression Analysis using DESeq2

This code performs differential gene expression analysis on the dataset GSE104075 using the DESeq2 package in R. It includes data preprocessing, normalization, differential expression analysis, visualization, and functional enrichment analysis.

## Getting Started
Install the required R packages: DESeq2, RColorBrewer, gplots, GO.db, and GOfuncR. 
   
Code Explanation
The code can be divided into the following sections:

1) Loading and Preprocessing Data:
Reads the count data from vehicle_drug_feature_counts.txt and selects relevant columns.
Performs basic quality checks on the count data.

2) Filtering Low-Expressed Genes:
Removes genes with low expression (bottom 25%) based on a threshold of less than 16 reads.

3) Normalization and Differential Expression Analysis:
Constructs a DESeqDataSet object from the filtered data.
Applies normalization and performs differential expression analysis using DESeq2.

4) Extracting Significantly Differentially Expressed Genes:
Filters the results based on adjusted p-values and selects genes with p-value < 0.01.

5) Visualization:
Generates a heatmap of normalized read counts.
Creates a volcano plot to visualize the significance and fold changes of differentially expressed genes.

6) Functional Enrichment Analysis:
Performs Gene Ontology (GO) enrichment analysis using the go_enrichment function.
Matches enriched GO terms with the GO term information from the annotation file.
Outputs a subset of enriched GO terms with adjusted p-values < 0.1.

7) Bar Plot for GO Enrichment:
Creates a bar plot of enriched GO terms, highlighting the enrichment levels.

Output Files
The code generates the following output files in the same directory:

resSig.sorted.tab.txt: A tab-separated file containing the sorted significant results of differentially expressed genes.
heatmap.pdf: A PDF file containing the heatmap plot of normalized read counts.
volcano.pdf: A PDF file containing the volcano plot of differentially expressed genes.
GO_barplot.pdf: A PDF file containing the bar plot of enriched GO terms.
