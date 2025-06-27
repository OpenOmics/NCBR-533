#!/bin/bash
#SBATCH --cpus-per-task=4 
#SBATCH --mem=72g
#SBATCH --time=9-00:00:00
#SBATCH --parsable
#SBATCH -J "mk_staph_blastndb"
#SBATCH --mail-type=BEGIN,END,FAIL
set -eu

# Load dependencies
module purge
module load blast/2.15.0+

# Functions
function timestamp() { date +"%Y-%m-%d_%H-%M-%S"; }

# makeblastdb needs a single fasta file
# as input, need to concat the genomic 
# and plasmid fasta files
echo "[$(timestamp)] Started to concatenate genomic and plasimd fasta files"
# Update the sequence identifiers
# replace spaces with underscores
# so the entire seqid is retained
cat \
    staphylococcus_aureus_genomic.fa \
    staphylococcus_aureus_plasmids.fa \
| sed '/^>/s/ /_/g' > staphylococcus_aureus_combined.fa

# Create a blastn database from the combined
# genomic and plasimd sequence fasta file
echo "[$(timestamp)] Started to build Staphylococcus aureus combined blastn db"
makeblastdb \
    -in staphylococcus_aureus_combined.fa \
    -dbtype nucl \
    -out "staphylococcus_aureus_db" \
    -title "Staphylococcus_aureus_genomic_and_plasmid_sequences" \
    -logfile create_blastn_db.log \
    -max_file_sz 4GB

echo "[$(timestamp)] Done creating the Staphylococcus aureus blastn db"
