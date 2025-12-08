#!/bin/bash

####pairwise_mafft####
##This generates a sequence alignment for each pair of TEs that are 75% similar. The output is then used to make codon alignments and get substitution scores using Run_pal2nal.sh

FASTA="rvt_aa_75.fasta"
PAIRS="species-unique_filtered_pim_pairs.tsv"
OUTDIR="pairwise_alignments_mafft"

mkdir -p "$OUTDIR"

# Load sequences into associative array for fast lookup
declare -A seqs
while read -r line; do
    if [[ "$line" =~ ^\> ]]; then
        id=${line#>}
        seqs["$id"]=""
        current_id="$id"
    else
        seqs["$current_id"]+=$line
    fi
done < <(grep -v "^$" "$FASTA")

# Function to write FASTA file for a pair
write_pair_fasta() {
    local qid=$1
    local tid=$2
    local outfile=$3

    echo ">$qid"
    echo "${seqs[$qid]}"
    echo ">$tid"
    echo "${seqs[$tid]}" > "$outfile"
}
