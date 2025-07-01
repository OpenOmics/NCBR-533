<div align="center">
   
  <h1>NCBR-533 🔬</h1>
  
  **_Data and script for processing NCBR-533_**

</div>

# Overview

This repository contains data and scripts for analyzing NCBR-533. 

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

To install the repository locally, you can use the following command:

```bash
# Clone the github repository
# and change your working directory
git clone https://github.com/OpenOmics/NCBR-533.git
cd NCBR-533/
```

## Setup your environment

To setup your environment and download any missing packages, you can use the following command:

```bash
# Install any missing or 
# required python packages
# in a virtual environment
python -m venv .venv
source .venv/bin/activate
pip install -U pip
pip install -r requirements.txt

# Install any missing R packages
./packages.R
```


## Reproduce the analyses

This is where you can add any steps to reproduce the analyses. For example, you can add the following command to run the script:

```bash
# Collapse and filter blastn results,
# the blastn db contains 130012 genomic
# and plasmid accession id representing
# all the Staphylococcus aureus genomic
# sequences on NCBI and all plasmid
# sequences on PLSDB
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
```
