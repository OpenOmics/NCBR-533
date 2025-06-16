#!/bin/bash
set -eu

# Thresholds
LEN_THRESHOLD="${1:-"0.95"}"  # Default: 0.95
PID_THRESHOLD="${2:-"95.00"}"  # Default: 95
# File to parse
BLAST_RESULTS="$3"


# Get size of blastdb
# ntargets=$(grep '^>' ../blastdb/staphylococcus_aureus_combined.fa | wc -l)
ntargets="8292985"   # pre-computed, takes forever to run
# nt_genomic=$(grep '^>' ../blastdb/staphylococcus_aureus_combined.fa | awk -F '|' '{print $1}' | grep '_genomic$' | wc -l)
nt_genomic="8290867" # pre-computed, takes forever to run
# nt_plasmid=$(grep '^>' ../blastdb/staphylococcus_aureus_combined.fa | awk -F '|' '{print $1}' | grep '_plasmid$' | wc -l)
nt_plasmid="2118"    # pre-computed, takes forever to run

# Get number of accession_ids in the blastn results
# n_acc=$(grep '^>' ../blastdb/staphylococcus_aureus_combined.fa | awk -F '|' '{print $1}' | sort | uniq | wc -l)
n_acc="130012"         # pre-computed, takes forever to run
# n_acc_genomic=$(grep '^>' ../blastdb/staphylococcus_aureus_combined.fa | awk -F '|' '{print $1}' | grep '_genomic$' | sort | uniq | wc -l)
n_acc_genomic="127894" # pre-computed, takes forever to run
# n_acc_plasmid=$(grep '^>' ../blastdb/staphylococcus_aureus_combined.fa | awk -F '|' '{print $1}' | grep '_plasmid$' | sort | uniq | wc -l)
n_acc_plasmid="2118"   # pre-computed, takes forever to run

echo -e "Query\tNHits_len${LEN_THRESHOLD}_pident${PID_THRESHOLD}\tNAccessionID_genomic+plasmid\tNAccessionID_genomic\tNAccessionID_plasmid\t%AccessionIDtoTargets"
while read gene length; do
    # Number hits (genomic+plasmid)
    len_filter=$(echo "$length*${LEN_THRESHOLD}" | bc)
    n_total_hits=$(awk -F '\t' -v GENE="$gene" -v LEN_FILTER="$len_filter" -v PID_FILTER="${PID_THRESHOLD}" '$1==GENE && $4+0.0>=LEN_FILTER && $3+0.0>=PID_FILTER {print}' "${BLAST_RESULTS}" | wc -l)
    # Number genomic hits
    n_genomic_hits=$(awk -F '\t' -v GENE="$gene" -v LEN_FILTER="$len_filter" -v PID_FILTER="${PID_THRESHOLD}" '$1==GENE && $4+0.0>=LEN_FILTER && $3+0.0>=PID_FILTER {print}' "${BLAST_RESULTS}" | cut -f2 | awk -F '|' '{print $1}' | grep '_genomic$' | wc -l)
    # Number plasmid hits
    n_plasmid_hits=$(awk -F '\t' -v GENE="$gene" -v LEN_FILTER="$len_filter" -v PID_FILTER="${PID_THRESHOLD}" '$1==GENE && $4+0.0>=LEN_FILTER && $3+0.0>=PID_FILTER {print}' "${BLAST_RESULTS}" | cut -f2 | awk -F '|' '{print $1}' | grep '_plasmid$' | wc -l)
    # Report both genomic and plasmid hits
    pct=$(echo -e "${n_total_hits}\t${ntargets}" | awk -F '\t' '{printf "%.2f\n", (($1/$2)*100.0)}')
    pct2=$(echo -e "${n_total_hits}\t${n_acc}" | awk -F '\t' '{printf "%.2f\n", (($1/$2)*100.0)}')
    echo -e "${gene}_genomic+plasmid_hits\t${n_total_hits}\t${n_acc}\t${n_acc_genomic}\t${n_acc_plasmid}\t${pct2}"
    # Report genomic hits
    pct=$(echo -e "${n_genomic_hits}\t${nt_genomic}" | awk -F '\t' '{printf "%.2f\n", (($1/$2)*100.0)}')
    pct2=$(echo -e "${n_genomic_hits}\t${n_acc_genomic}" | awk -F '\t' '{printf "%.2f\n", (($1/$2)*100.0)}')
    echo -e "${gene}_genomic_hits\t${n_genomic_hits}\t${n_acc}\t${n_acc_genomic}\t${n_acc_plasmid}\t${pct2}"
    # Report plasmid hits
    pct=$(echo -e "${n_plasmid_hits}\t${nt_plasmid}" | awk -F '\t' '{printf "%.2f\n", (($1/$2)*100.0)}')
    pct2=$(echo -e "${n_plasmid_hits}\t${n_acc_plasmid}" | awk -F '\t' '{printf "%.2f\n", (($1/$2)*100.0)}')
    echo -e "${gene}_plasmid_hits\t${n_plasmid_hits}\t${n_acc}\t${n_acc_genomic}\t${n_acc_plasmid}\t${pct2}"
done < <(paste - - < staphylococcus_aureus_query_genes_of_interest.fa | awk -F '\t' -v OFS='\t' '{print $1,length($2)}' | awk -v OFS='\t' '{print $1,$NF}' | sed 's/^>//g')
