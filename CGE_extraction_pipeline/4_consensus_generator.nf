// params.species = "setmeup"
params.genomefolder = "${projectDir}/genomes/${params.species}/ncbi_dataset/data"
params.genomestring = file("$projectDir/genomes/${params.species}/ncbi_dataset/data/*/*.fna")
params.genome = "${params.genomestring[0]}"
params.orf2folder = "$projectDir/analysis/${params.species}/ORF2"
params.repeats = "$projectDir/repeats/Jockey_RVT_over200_85cluster.fa"
params.dbname = "${params.species}"
params.conda_env = '/hpcfs/users/a1749756/myconda/envs/bedtools'
params.mchelperfolder = "$projectDir/analysis/${params.species}/repeatcuration/MCHelper"
params.extendfolder = "$projectDir/analysis/${params.species}/repeatcuration/extension"

Channel.fromPath("${params.extendfolder}/*")
    .set { folder_ch }

process CONDA{

    conda params.conda_env
    script:
    """
    source activate ${params.conda_env}
    """
}


process CONFIRMJOCKEYSTATUS {
    input:
    path orf2folder

    output:
    path orf2folder 

    script:
    """
    grep ">" ${params.mchelperfolder}/curated_sequences_NR.fa | cut -c 2-15 | sort | uniq > ${params.mchelperfolder}/curated_sequences_NR.txt

    while IFS= read -r line
    do
        for word in \$line
        do
            if grep -qF "\$word" ${params.orf2folder}/85/${params.species}_85_jockeyRMSKlist.txt
            then
                echo "\$word" >> ${params.mchelperfolder}/Jockey_curated_sequences_NR.txt
            else
                echo "\$word" >> ${params.mchelperfolder}/notJockey_curated_sequences_NR.txt
            fi
        done
    done < ${params.mchelperfolder}/curated_sequences_NR.txt
    """
}


process FOLDERSETUP {

    input:
    path orf2folder

    output:
    path orf2folder 

    tag "makin folders :)"

    script:
    """
    # Check if the folder exists
    if [ ! -d "${projectDir}/analysis/${params.species}/repeatcuration/extension" ]; then
        # If it doesn't exist, create it
        mkdir -p "${projectDir}/analysis/${params.species}/repeatcuration/extension" 
        seqkit split -i -O ${params.extendfolder} ${params.mchelperfolder}/curated_sequences_NR.fa
        rename 's/curated_sequences_NR.part_//' ${params.extendfolder}/*
        find ${params.extendfolder} -type f -exec bash -c 'file={}; dir=\${file%.*}; mkdir -p \$dir; mv \$file \$dir' \\;
    fi
    """
}


process  EXTENSION{
    input:
    path extend

    output:
    path extend
    
    script:
    """   
for dir in ${params.extendfolder}/*
do
    if [[ -d \$dir ]]; then
        sequence=\$(basename \$dir)
        if [[ -f \$dir/\$sequence.fa ]]; then
            seq_length=\$(seqkit fx2tab -l -i \$dir/\$sequence.fa | cut -f 4)
            if (( seq_length >= 8000 )); then
                bash make_fasta_from_blast.sh ${params.genome} \$dir/\$sequence.fa 0 0
                cp \$sequence.fa.blast.bed.fa \$dir
            elif (( seq_length < 8000 )); then
                half_diff=\$(( (8000 - seq_length) / 2 ))
                bash make_fasta_from_blast.sh ${params.genome} \$dir/\$sequence.fa 0 \$half_diff 
                cp \$sequence.fa.blast.bed.fa \$dir
            fi           
        fi
    fi
done    
"""
}

process  MAFFT{
    cpus { 14 }
    time { '18h' }
    memory { '120.GB' }


    input:
    path extend

    output:
    path extend
    
    script:
    """   
for dir in ${params.extendfolder}/*
do
    if [[ -d \$dir ]]; then
        sequence=\$(basename \$dir)
        if [[ -f \$dir/\$sequence.fa ]]; then
          count=\$(grep ">" \$dir/\$sequence.fa.blast.bed.fa | wc -l)
            if (( count > 100 )); then
                bash ready_for_MSA.sh \$dir/\$sequence.fa.blast.bed.fa 100 75
                mafft --auto --parttree --thread 8 \$dir/\$sequence.fa.blast.bed.fa.rdmSubset.fa > \$dir/\$sequence.fa.blast.bed.fa.unc1.maf
            else
                mafft --auto --parttree --thread 8 \$dir/\$sequence.fa.blast.bed.fa > \$dir/\$sequence.fa.blast.bed.fa.unc1.maf
            fi              
        fi
    fi
done    
"""
}



/*

To rerun on a single sequence, use the following. droNov_TE00001#CLASSI__LINE__CR1.fa is used as an example:

count=\$(grep ">" \$dir/\$sequence.fa.blast.bed.fa | wc -l)
            if (( count >= 100 )); then
                bash ready_for_MSA.sh \$dir/\$sequence.fa.blast.bed.fa 100 75
                mafft --thread 8 \$dir/\$sequence.fa.blast.bed.fa.rdmSubset.fa > \$dir/\$sequence.fa.blast.bed.fa.unc1.maf
            else
                mafft --thread 8 \$dir/\$sequence.fa.blast.bed.fa > \$dir/\$sequence.fa.blast.bed.fa.unc1.maf
            fi

for dir in ~/*
do
    if [[ -d $dir ]]; then
        sequence=$(basename $dir)
        if [[ -f $dir/$sequence.fa ]]; then
            seq_length=$(seqkit fx2tab -l -i $dir/$sequence.fa | cut -f 4)
            if (( seq_length < 5500 )); then
                # Store the sequence length in a variable
                short_sequence_length=$seq_length
            fi
        fi
    fi
done


/*
process  {
    input:
    path orf2folder

    output:
    path orf2folder
    
    script:
    """
    seqkit fx2tab -l -i droNov_TE00001#CLASSI__LINE__CR1/droNov_TE00001#CLASSI__LINE__CR1.fa | cut -f 4
    bash ready_for_MSA.sh droNov_TE00001#CLASSI__LINE__CR1.fa 
    bash make_fasta_from_blast.sh ../../../../../genomes/droNov/ncbi_dataset/data/GCF_003342905.1/GCF_003342905.1_droNov1_genomic.fna droNov_TE00001#CLASSI__LINE__CR1.fa 0 500
    grep ">" droNov_TE00001#CLASSI__LINE__CR1.fa.blast.bed.fa | wc -l
    # if bigger than 100, ready_for_MSA.sh 
    bash ready_for_MSA.sh droNov_TE00001#CLASSI__LINE__CR1.fa.blast.bed.fa 100 75
    mafft --thread 8 droNov_TE00001#CLASSI__LINE__CR1.fa.blast.bed.fa.rdmSubset.fa > droNov_TE00001#CLASSI__LINE__CR1.unc1.maf
    # manually define 3 prime end 
    mafft --thread 8 droNov_TE00001#CLASSI__LINE__CR1.cur1.maf > droNov_TE00001#CLASSI__LINE__CR1.cur1.maf.realigned
    singularity exec /hpcfs/users/a1749756/tools/docker/tcoffee_Version_13.46.0.919e8c6b.sif t_coffee -other_pg seq_reformat \
    -in droNov_TE00001#CLASSI__LINE__CR1.cur1.maf.realigned -action +rm_gap 85 > droNov_TE00001#CLASSI__LINE__CR1.cur1.tc
    cons -sequence droNov_TE00001#CLASSI__LINE__CR1.cur1.tc -outseq droNov_TE00001#CLASSI__LINE__CR1.cur1.cons -plurality 10
    """
}
*/


workflow {
    CONDA
    confirmjs_ch = CONFIRMJOCKEYSTATUS(params.orf2folder)
    foldersetup_ch = FOLDERSETUP(confirmjs_ch)
    firstextend_ch= EXTENSION(foldersetup_ch)
    mafft_ch= MAFFT(firstextend_ch)

}

/*   
 foldersetup_ch = FOLDERSETUP()
    PARTEST(foldersetup_ch)
    }
*/
