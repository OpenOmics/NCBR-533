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

# Create a strain-specific FASTA file
# alignment sequence for each gene.
# Filter thresholds:
# - Alignment length fraction: 0.95
# - Alignment percent identity: 95.0
mkdir -p results/fasta_sequences/stringent
# Group1: MSSA (mecA negative strains)
awk -F '\t' -v OFS='\t' 'NR>1 && $3=="0" {print}' results/gene_detection.tsv \
  | cut -f1 \
> results/MSSA_seqids.txt
# Group 2: MRSA-mecI-positive (mecA and mecI positive strains)
awk -F '\t' -v OFS='\t' 'NR>1 && ($3=="1" && $4=="1") {print}' results/gene_detection.tsv \
  | cut -f1 \
> results/MRSA-mecI-positive_seqids.txt
# Group 3: MRSA (mecA positive but mecI negative strains)
awk -F '\t' -v OFS='\t' 'NR>1 && ($3=="1" && $4=="0") {print}' results/gene_detection.tsv \
  | cut -f1 \
> results/MRSA-mecI-negative_seqids.txt
# Extract gene alignment sequence of
# each strain of interest
echo '#!/usr/bin/env bash' > ${PWD}/results/fasta_sequences/stringent/extract_sequences.swarm
for strain in "MSSA" "MRSA-mecI-positive" "MRSA-mecI-negative"; do
  echo "# ${strain} strains" >> ${PWD}/results/fasta_sequences/stringent/extract_sequences.swarm
  for gene in "BlaI" "BlaR1" "BlaZ" "MecA" "MecI" "Stk1" "Stp1"; do
    echo "${PWD}/scripts/extract_blastn_alignment_sequence.py --input-blast-result ${PWD}/results/${gene}_blastn_results.tsv --input-strain-ids ${PWD}/results/${strain}_seqids.txt --input-strain-fasta ${PWD}/data/staphylococcus_aureus_genomic.fa --output-strain-fasta ${PWD}/results/fasta_sequences/stringent/${strain}/${strain}_${gene}_aligned_sequence.fa" >> ${PWD}/results/fasta_sequences/stringent/extract_sequences.swarm
  done
done && chmod +x ${PWD}/results/fasta_sequences/stringent/extract_sequences.swarm
swarm --logdir "${PWD}/results/fasta_sequences/stringent/" --module 'python/3.12' -t 2 -g 128 -b 5 -f ${PWD}/results/fasta_sequences/stringent/extract_sequences.swarm

# Filter thresholds:
# - Alignment length fraction: 0.95
# - Alignment percent identity: 80.0
mkdir -p results/fasta_sequences/relaxed
# Group1: MSSA (mecA negative strains)
awk -F '\t' -v OFS='\t' 'NR>1 && $3=="0" {print}' results/relaxed/gene_detection.tsv \
  | cut -f1 \
> results/relaxed/MSSA_seqids.txt
# Group 2: MRSA-mecI-positive (mecA and mecI positive strains)
awk -F '\t' -v OFS='\t' 'NR>1 && ($3=="1" && $4=="1") {print}' results/relaxed/gene_detection.tsv \
  | cut -f1 \
> results/relaxed/MRSA-mecI-positive_seqids.txt
# Group 3: MRSA (mecA positive but mecI negative strains)
awk -F '\t' -v OFS='\t' 'NR>1 && ($3=="1" && $4=="0") {print}' results/relaxed/gene_detection.tsv \
  | cut -f1 \
> results/relaxed/MRSA-mecI-negative_seqids.txt
# Extract gene alignment sequence of
# each strain of interest
echo '#!/usr/bin/env bash' > ${PWD}/results/fasta_sequences/relaxed/extract_sequences.swarm
for strain in "MSSA" "MRSA-mecI-positive" "MRSA-mecI-negative"; do
  echo "# ${strain} strains" >> ${PWD}/results/fasta_sequences/relaxed/extract_sequences.swarm
  for gene in "BlaI" "BlaR1" "BlaZ" "MecA" "MecI" "Stk1" "Stp1"; do
    echo "${PWD}/scripts/extract_blastn_alignment_sequence.py --input-blast-result ${PWD}/results/relaxed/${gene}_blastn_results.tsv --input-strain-ids ${PWD}/results/relaxed/${strain}_seqids.txt --input-strain-fasta ${PWD}/data/staphylococcus_aureus_genomic.fa --output-strain-fasta ${PWD}/results/fasta_sequences/relaxed/${strain}/${strain}_${gene}_aligned_sequence.fa" >> ${PWD}/results/fasta_sequences/relaxed/extract_sequences.swarm
  done
done && chmod +x ${PWD}/results/fasta_sequences/relaxed/extract_sequences.swarm
swarm --logdir "${PWD}/results/fasta_sequences/stringent/" --module 'python/3.12' -t 2 -g 128 -b 5 -f ${PWD}/results/fasta_sequences/relaxed/extract_sequences.swarm
```

## Methods

### Staphylococcus aureus gene detection analyses

To systematically assess the presence of key antibiotic resistance and regulatory genes in *Staphylococcus aureus*, all publicly available *S. aureus* genomes deposited in the National Center for Biotechnology Information (NCBI) database<sup>1</sup> were downloaded along with the sequences of the following key genes of interest: **BlaZ**, **MecA**, **MecI**, **Stk1**, **Stp1**, **BlaR1**, and **BlaI**. The goal was to determine the presence or absence of these genes across the comprehensive dataset of *S. aureus* genomes. In total, 127894 unique accession identifiers corresponding to *S. aureus* genome assemblies were retrieved from NCBI. These genome assemblies were combined to build a custom BLAST nucleotide database using makeblastdb from the BLAST+ suite<sup>2</sup> (version 2.15.0+).

The full-length nucleotide sequences of the seven target genes of interest were queried against the custom *S. aureus* BLAST database using blastn with the `-max_target_seqs` option to ensure the alignment results of each query-target pair was reported. To account for the possibility of multiple hits per accession (e.g., due to fragmented assemblies or gene copy number variation), the blastn results were first collapsed by accession identifier. For each gene-accession pair, only the single top-scoring hit was retained, prioritizing alignments with the greatest percent identity and alignment length. Hits were then filtered using stringent criteria to minimize false positives: only alignments with *(i)* an alignment length fraction ≥ 0.95 (i.e., ≥ 95% of the gene length aligned) and *(ii)* percent nucleotide identity ≥ 95% were considered indicative of gene presence. The alignment length fraction was calculated as the aligned region divided by the full length of the query gene. A binary gene presence matrix (containing 127894 accession identifiers by the 7 genes of interest) was created, where presence of a gene was encoded as `1` and the absence of a gene was encoded as a `0`.

## References

<sup>**1.** National Center for Biotechnology Information (NCBI)[Internet]. Bethesda (MD): National Library of Medicine (US), National Center for Biotechnology Information; [1988] – [cited 2025 June 01]. Available from: https://www.ncbi.nlm.nih.gov/</sup>    

<sup>**2.** Camacho, C., Coulouris, G., Avagyan, V. et al. BLAST+: architecture and applications. BMC Bioinformatics 10, 421 (2009). https://doi.org/10.1186/1471-2105-10-421</sup>    