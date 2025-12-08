####Run_pal2nal####
##this takes the output for MAFFT_TE-pairs.sh and generates codon alignments which pal2nal needs to calculate substitution scores

#!/bin/bash

mkdir -p codon_alignments_rerun

for aln in pairwise_alignments_mafft/*.aln; do
    # Get the base name (e.g., AldAff_CR1-1_RVT_vs_AlbGor_CR1-1_RVT)
    base=$(basename "$aln" .aln)

    # Set paths
    pep_file="/home/noz/Downloads/pairwise_alignments_mafft/${base}.aln"
    nt_file="/home/noz/Downloads/pair_fastas/${base}.fasta"
    out_file="/home/noz/Downloads/codon_alignments_rerun/${base}_codon_alignment.phy"

    # Run pal2nal if files exist
    if [[ -f "$pep_file" && -f "$nt_file" ]]; then
        perl pal2nal.v14/pal2nal.pl "$pep_file" "$nt_file" -output fasta -nogap > "$out_file"
    else
        echo "Missing: $pep_file or $nt_file"
    fi
done
