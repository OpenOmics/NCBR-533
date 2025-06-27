#!/bin/bash
#SBATCH --cpus-per-task=12 
#SBATCH --mem=256g
#SBATCH --time=10-00:00:00
#SBATCH --parsable
#SBATCH -J "find_genes_staph_blastn"
#SBATCH --mail-type=BEGIN,END,FAIL
set -eu

# Load dependencies
module purge
module load blast/2.15.0+

# Functions
function timestamp() { date +"%Y-%m-%d_%H-%M-%S"; }

# Run query gene list against the blastn 
# Staph strain genome + plasmid database
ntargets=$(grep '^>' ../blastdb/staphylococcus_aureus_combined.fa | wc -l)
ntargets=$((ntargets + 1))   # +1 to ensure we get all targets
echo "[$(timestamp)] Started blastn gene search against Staph genome + plasmid database (contains ${ntargets} strains)..."

blastn -query staphylococcus_aureus_query_genes_of_interest.fa \
    -db ../blastdb/staphylococcus_aureus_db \
    -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore" \
    -num_threads 12 \
    -out staphylococcus_aureus_genes_blastn_results.tsv \
    -max_target_seqs ${ntargets}

echo "[$(timestamp)] Finished running blastn against Staph genome + plasmid database. Results saved to staphylococcus_aureus_genes_blastn_results.tsv"
