#!/usr/bin/env Rscript

# Load necessary libraries
suppressMessages(library("ggpubr"))

# Misc helper functions 
err <- function(...){cat(sprintf(...), sep='\n', file=stderr())}
fatal <- function(...) {err(...); quit(status = 1)}

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
    fatal("Usage: eda_plots.R <path_to_blastn_results> <output_directory>")
}
blastn_file <- args[1]
output_directory <- args[2]

# Read in blastn results where the query
# is each gene of interest (BlaZ, MecA, MecI)
# and the subject is the genomic + plasmid
# sequences of all Staphylococcus A. strains
# Here is the expected format blastn output file:
# 1  qseqid                                    BlaZ
# 2  sseqid  NZ_CP017095.1_SAMN04966147_plasmid|...
# 3  pident                                 100.000
# 4  length                                     846
# 5  mismatch                                     0
# 6  gapopen                                      0
# 7  qstart                                       1
# 8  qend                                       846
# 9  sstart                                   11212
# 10 send                                     12057
# 11 evalue                                     0.0
# 12 bitscore                                  1563
blastn_results <- read.table(
    file = blastn_file,
    sep = "\t",
    header = TRUE,
    stringsAsFactors = TRUE,
    comment.char=""
)
# Add new column for query alignment length
blastn_results$qlen <- (blastn_results$qend - blastn_results$qstart) + 1
# Add new column for subject alignment length
blastn_results$slen <- (blastn_results$send - blastn_results$sstart) + 1
blastn_results$abslen <- abs((blastn_results$send - blastn_results$sstart) + 1)
# Create a density plot of percent identity
# for each query gene
important_features <- c(
    "pident" = "Percentage Identitical matches", 
    "length" = "Alignment Length", 
    "mismatch" = "Number of mismatches",
    "gapopen" = "Number of gap openings",
    "qlen" = "Query/gene length (i.e qend-qstart)",
    "slen" = "Subject/strain length (i.e send-sstart)",
    "abslen" = "Subject/strain abs length (i.e |send-sstart|)"
)
# Gene  len
# BlaZ	846
# MecA	2007
# MecI	372
gene_length <- c("BlaZ"=846, "MecA"=2007, "MecI"=372)

for (feature in names(important_features)){
    # Create a folder for each feature
    outdir <- file.path(output_directory, "figs", feature)
    dir.create(outdir, showWarnings = FALSE)
    setwd(outdir)
    # Create a feature box plot with all genes
    fig <- ggboxplot(
        blastn_results,
        x = "qseqid",
        y = feature,
        fill = "qseqid",
        palette = c("#00AFBB", "#E7B800", "#FC4E07"),
        xlab = paste0(important_features[feature]),
    )
    # Write feature boxplot to file
    ggexport(
        fig,
        filename = paste("staphylococcus_aureus_",feature,"-all-genes-boxplot.png", sep=""),
        width = 1800,
        height = 1000,
        res = 300
    )

    for (gene in levels(blastn_results$qseqid)){
        # Get max value
        max_val <- max(blastn_results[blastn_results$qseqid==gene,feature])
        print(paste(0.95*max_val, max_val))
        # Create a feature histogram for each gene
        fig <- gghistogram(
            blastn_results[blastn_results$qseqid==gene,], 
            x = feature,
            y = "count",
            add = "mean",
            rug = TRUE,
            bins = gene_length[gene],
            fill = "qseqid",
            palette = c("#00AFBB"),
            xlab = paste0(important_features[feature]),
        ) + coord_cartesian(xlim = c(max_val * 0.95, max_val))
        # Write histogram to file
        ggexport(
            fig,
            filename = paste("staphylococcus_aureus_",feature,"-",gene,"-gene-histogram.png", sep=''),
            width = 1800,
            height = 800,
            res = 300
        )
    }
}
