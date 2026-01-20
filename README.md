<div align="center">
   
  <h1>NCBR-533 🔬</h1>
  
  **_Data and script for processing Staphylococcus aureus data-mining project_**

</div>

# Overview

This repository contains data and scripts for identifying antibiotic resistance genes in all publicly available Staphylococcus aureus genomic and plasmid sequences. 

Any data accompanying this project can be stored in the `data/` directory, while any scripts used to process the data can be stored in the `scripts/` directory.

> [!NOTE]  
> _**By default**, any data or files added to the `data/` directory are ignored by git due to our `.gitignore`, so you can store large files or data here without worrying about them being uploaded to Github._ If you would to upload a small data file (<5MB) to Github, you can stage the file using force, `-f`, option to the `git add`:  
> ```bash
> # Stage file for commit
> git add -f data/counts.tsv
> # Commit the file to history
> git commit -m "Adding small counts matrix"
> ```


## Installation

To install the repository locally, please run the following commands:

```bash
# Clone the github repository
# and change your working directory
git clone https://github.com/OpenOmics/NCBR-533.git
cd NCBR-533/
```

## Setup your environment

To setup your environment and download any missing packages, please run the following commands:

```bash
# Install any missing or
# required python packages
# in a virtual environment
python3 -m venv .venv
source .venv/bin/activate
pip install -U pip
pip install -r requirements.txt

# Install any missing R packages
./packages.R
```


## Reproduce the analyses

To reproduce any downstream analyses, please run the following commands below.

```bash
# Collapse and filter blastn results,
# the blastn db contains 130012 genomic
# and plasmid accession id representing
# all the Staphylococcus aureus genomic
# sequences on NCBI and all plasmid
# sequences on PLSDB. Filter thresholds:
# - Alignment length fraction: 0.95
# - Alignment percent identity: 95.0
./scripts/parse_blastn_results.py \
  -i data/staphylococcus_aureus_genes_blastn_results-with-header.tsv \
  -s data/blastn_seqids.txt \
  -l 0.95 -p 95.0 \
  -o results/blastn_filtered
# Get a high-level overview of the
# collapsed/filtered results
./scripts/get_overview.sh \
  0.95 95.0 \
  results/blastn_filtered_collapsed_results.tsv \
  data/staphylococcus_aureus_query_genes_of_interest.fa \
> results/overview.tsv
# Parse out each gene's results
# to create a seperate sheet in
# the excel spreadsheet, only
# parse out genomic alignments
while read gene; do
  echo "# Parsing ${gene} from collapsed results";
  awk -F '\t' -v GENE="${gene}" \
    'NR==1 || ($1==GENE && $2 ~ /_genomic$/) {print}' results/blastn_filtered_collapsed_results.tsv \
  > results/${gene}_blastn_results.tsv;
done < <(grep '^>' data/staphylococcus_aureus_query_genes_of_interest.fa | sed 's/^>//g' | awk '{print $1}')
# Create an excel spreadsheet with
# the overview, presences/absences
# file, and each genes blastn results
awk -F '\t' 'NR==1 || $1 ~ /_genomic$/ {print}'  results/blastn_filtered_gene_presence_absence.tsv > results/gene_detection.tsv
./scripts/files2spreadsheet.py \
  -i results/overview.tsv results/gene_detection.tsv results/*_blastn_results.tsv \
  -o results/NCBR-533_Staphylococcus_aureus-genomic_7gene_presence_0.95-length_95-pidentity.xlsx -a

# Collapse and filter blastn results,
# the blastn db contains 130012 genomic
# and plasmid accession id representing
# all the Staphylococcus aureus genomic
# sequences on NCBI and all plasmid
# sequences on PLSDB. Filter thresholds:
# - Alignment length fraction: 0.95
# - Alignment percent identity: 80.0
./scripts/parse_blastn_results.py \
  -i data/staphylococcus_aureus_genes_blastn_results-with-header.tsv \
  -s data/blastn_seqids.txt \
  -l 0.95 -p 80.0 \
  -o results/relaxed/blastn_filtered
# Get a high-level overview of the
# collapsed/filtered results
./scripts/get_overview.sh \
  0.95 80.0 \
  results/relaxed/blastn_filtered_collapsed_results.tsv \
  data/staphylococcus_aureus_query_genes_of_interest.fa \
> results/relaxed/overview.tsv
# Parse out each gene's results
# to create a seperate sheet in
# the excel spreadsheet, only
# parse out genomic alignments
while read gene; do
  echo "# Parsing ${gene} from collapsed results";
  awk -F '\t' -v GENE="${gene}" \
    'NR==1 || ($1==GENE && $2 ~ /_genomic$/) {print}' results/relaxed/blastn_filtered_collapsed_results.tsv \
  > results/relaxed/${gene}_blastn_results.tsv;
done < <(grep '^>' data/staphylococcus_aureus_query_genes_of_interest.fa | sed 's/^>//g' | awk '{print $1}')
# Create an excel spreadsheet with
# the overview, presences/absences
# file, and each genes blastn results
awk -F '\t' 'NR==1 || $1 ~ /_genomic$/ {print}'  results/relaxed/blastn_filtered_gene_presence_absence.tsv > results/relaxed/gene_detection.tsv
./scripts/files2spreadsheet.py \
  -i results/relaxed/overview.tsv results/relaxed/gene_detection.tsv results/relaxed/*_blastn_results.tsv \
  -o results/NCBR-533_Staphylococcus_aureus-genomic_7gene_presence_0.95-length_80-pidentity.xlsx -a
```

## Methods

### Staphylococcus aureus gene detection analyses

To systematically assess the presence of key antibiotic resistance and regulatory genes in *Staphylococcus aureus*, all publicly available *S. aureus* genomes deposited in the National Center for Biotechnology Information (NCBI) database<sup>1</sup> were downloaded along with the sequences of the following key genes of interest: **BlaZ**, **MecA**, **MecI**, **Stk1**, **Stp1**, **BlaR1**, and **BlaI**. The goal was to determine the presence or absence of these genes across the comprehensive dataset of *S. aureus* genomes. In total, 127894 unique accession identifiers corresponding to *S. aureus* genome assemblies were retrieved from NCBI. These genome assemblies were combined to build a custom BLAST nucleotide database using makeblastdb from the BLAST+ suite<sup>2</sup> (version 2.15.0+).

The full-length nucleotide sequences of the seven target genes of interest were queried against the custom *S. aureus* BLAST database using blastn with the `-max_target_seqs` option to ensure the alignment results of each query-target pair was reported. To account for the possibility of multiple hits per accession (e.g., due to fragmented assemblies or gene copy number variation), the blastn results were first collapsed by accession identifier. For each gene-accession pair, only the single top-scoring hit was retained, prioritizing alignments with the greatest percent identity and alignment length. Hits were then filtered using stringent criteria to minimize false positives: only alignments with *(i)* an alignment length fraction ≥ 0.95 (i.e., ≥ 95% of the gene length aligned) and *(ii)* percent nucleotide identity ≥ 95% were considered indicative of gene presence. The alignment length fraction was calculated as the aligned region divided by the full length of the query gene. A binary gene presence matrix (containing 127894 accession identifiers by the 7 genes of interest) was created, where presence of a gene was encoded as `1` and the absence of a gene was encoded as a `0`.

## References

<sup>**1.** National Center for Biotechnology Information (NCBI)[Internet]. Bethesda (MD): National Library of Medicine (US), National Center for Biotechnology Information; [1988] – [cited 2025 June 01]. Available from: https://www.ncbi.nlm.nih.gov/</sup>    

<sup>**2.** Camacho, C., Coulouris, G., Avagyan, V. et al. BLAST+: architecture and applications. BMC Bioinformatics 10, 421 (2009). https://doi.org/10.1186/1471-2105-10-421</sup>    
