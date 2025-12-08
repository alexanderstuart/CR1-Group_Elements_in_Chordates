// params.species = "setmeup"
params.genomefolder = "${projectDir}/genomes/${params.species}/ncbi_dataset/data"
params.genomestring = file("$projectDir/genomes/${params.species}/ncbi_dataset/data/*/*.fna")
params.genome = "${params.genomestring[0]}"
params.orf2folder = "$projectDir/analysis/${params.species}/ORF2"
params.repeats = "$projectDir/repeats/CGE_RVT_over200_85cluster.fa"
params.dbname = "${params.species}"
params.conda_env = '/hpcfs/users/a1749756/myconda/envs/bedtools'
params.mchelper3folder = "$projectDir/analysis/${params.species}/repeatcuration/MCHelper3"
params.extendfolder = "$projectDir/analysis/${params.species}/repeatcuration/extension"
params.curatedtefolder = "$projectDir/analysis/${params.species}/repeatcuration/${params.species}_curatedTEs"
params.cronefolder = "$projectDir/analysis/${params.species}/repeatcuration/${params.species}_curatedTEs/${params.species}_CR1s"
params.orffolder = "$projectDir/analysis/${params.species}/repeatcuration/${params.species}_curatedTEs/${params.species}_ORFs"
params.teaidfolder = "$projectDir/analysis/${params.species}/repeatcuration/${params.species}_curatedTEs/${params.species}_TEAid_plots"

Channel.fromPath("${params.extendfolder}/*")
    .set { folder_ch }


process CONDA{

    conda params.conda_env
    script:
    """
    source activate ${params.conda_env}
    """
}

process CHECKFILES {

    executor 'local'

    input:
    path extendfolder

    output:
    path extendfolder

    script:
    """
    # Debug print to check the value of extendfolder
    echo extendfolder: ${params.extendfolder}}

    # Enable nullglob option
    shopt -s nullglob

    # Get the list of .cur1.maf files in the extendfolder/* directory
    files=(${params.extendfolder}/*/*.cur1.maf)

    # Check if the files array is empty
    if [ \${#files[@]} -eq 0 ]; then
        exit 1
    else
        # If it's not, iterate over the files and check if they are empty
        for file in \${files[@]}; do
            if [ ! -s "\$file" ]; then
                exit 1
            fi
        done

    fi
     # Create a check file to use as process output
    touch ${params.extendfolder}/checkfile.txt
    """ 

}


process  MAFFT{
    cpus { 8 }
    time { '72h' }
    memory { '100.GB' }


    input:
    path extendfolder
    
    output:
    path extendfolder    
    script:
    """   
for dir in ${params.extendfolder}/*
do
    if [[ -d \$dir ]]; then
        sequence=\$(basename \$dir)
        if [[ -f \$dir/\$sequence.fa ]]; then
            # Find the cur1.maf file if it exists
            cur1_maf_file=\$(find \$dir -maxdepth 1 -name '*.cur1.maf' -type f -print -quit)

            if [[ -n "\$cur1_maf_file" ]]; then
                mafft --thread 8 "\$cur1_maf_file" > \$dir/\$sequence.cur1.maf.realigned
                count=\$(grep ">" \$dir/\$sequence.cur1.maf.realigned | wc -l)
                rounded_tcoffee_val=\$(awk -v count="\$count" 'BEGIN { printf "%.0f", count * 0.6 }')
                singularity exec --bind /scratchdata1:/scratchdata1 /hpcfs/users/a1749756/tools/docker/tcoffee_Version_13.46.0.919e8c6b.sif t_coffee -other_pg seq_reformat \
                    -in \$dir/\$sequence.cur1.maf.realigned -action +rm_gap "\$rounded_tcoffee_val" > \$dir/\$sequence.cur1.tc

                rounded_plurval=\$(awk -v count="\$count" 'BEGIN { printf "%.0f", count * 0.15 }')
                cons -sequence "\$dir/\$sequence.cur1.tc" -outseq "\$dir/\$sequence.cur1.cons" -plurality "\$rounded_plurval"
            fi
        fi
    fi
done


"""
}

process  RENAMING{


    input:
    path extendfolder
    output:
    path extendfolder    
    script:
"""
for dir in ${params.extendfolder}/*; do
    if [[ -d \$dir ]]; then
        sequence=\$(basename \$dir)
        if [[ -f \$dir/\$sequence.fa ]]; then
            filename=\$(basename "\$dir/\$sequence.fa" .cons)
            filename=\${filename%%#*}
            
            if [ -f "\$dir/\$sequence.cur1.cons" ]; then
            seqkit replace -p "n" -r "" -s < \$dir/\$sequence.cur1.cons | seqkit seq -w 60 |
                sed "s/>EMBOSS_001/>\$filename/" | sed "s/#.*//" > "\$dir/\${filename}_cons.fa"

                sequencelength=\$(grep -v ">" "\$dir/\${filename}_cons.fa" | \
                    tr -d '\n' | wc | awk '{print int(\$3*0.8 + 0.5)}')
                
                echo "#qseqid sseqid pident length mismatch qstart qend sstart send sstrand" > \$dir/\${filename}_blastn.out
                blastn -query "\$dir/\${filename}_cons.fa" -db ${params.genome} \
                    -outfmt "6 qseqid sseqid pident length mismatch qstart qend sstart send sstrand" \
                    -evalue 1e-20 | awk -v "ml=\$sequencelength" '{OFS="\t"; if (\$4 > ml) {print \$0}}' >> \$dir/\${filename}_blastn.out
                        
                awk '{OFS="\t"; if (\$1~/^#/ || \$3 <= 80) {} else { if (\$10~/plus/) {print \$2, \$8, \$9, \$1, \$3, "+"} \
                    else {print \$2, \$9, \$8, \$1, \$3, "-"}}}' < \$dir/\${filename}_blastn.out > \$dir/\${filename}_blastn.bed
                
                bedtools getfasta -fi ${params.genome} -fo \$dir/\${filename}_reps.fa -bed \$dir/\${filename}_blastn.bed -s 
            else
                echo "File \$dir/\$sequence.cur1.cons does not exist. Skipping..."
            fi
        fi
    fi
done
"""

}

process  GETORFS{
    executor 'local'

    input:
    path extendfolder
    output:
    path extendfolder    
    script:
    """   
for dir in ${params.extendfolder}/*
do
    if [[ -d \$dir ]]; then
        sequence=\$(basename \$dir)
        if [[ -f \$dir/\$sequence.fa ]]; then
            filename=\$(basename "\$dir/\$sequence.fa" .cons)
            filename=\${filename%%#*}
            cur1_maf_file=\$(find \$dir -maxdepth 1 -name '*.cur1.maf' -type f -print -quit)
            if [[ -n "\$cur1_maf_file" ]]; then
            representative_file=\$(find \$dir -maxdepth 1 -name '*_reps.fa' -type f -print -quit)
                if [[ -s "\$representative_file" ]]; then

                getorf -sequence \$dir/\${filename}_reps.fa -minsize 300 -outseq \$dir/\${filename}_orfs.fa -reverse N

                cd-hit -i \$dir/\${filename}_orfs.fa -o \$dir/\${filename}_clustered_orfs.fasta -c 0.75 -aL 0.8 -g 1 -d 0
                
                clstr_file="\$dir/\${filename}_clustered_orfs.fasta.clstr"
                awk '
                /^>Cluster/ {
                    if (NR > 1) clusters[cluster] = count;
                    cluster = \$2; count = 0;
                }
                !/^\\>/ { count++; }
                END {
                    clusters[cluster] = count;
                    for (cluster in clusters) print cluster, clusters[cluster];
                }' "\$clstr_file" | sort -k2,2nr > \$dir/\${filename}_sorted_orf_clusters.txt
                
                third_largest_size=\$(awk '\$2 >= 6 {print \$2}' \$dir/\${filename}_sorted_orf_clusters.txt | uniq -c | awk '{count+=\$1; if (count>=3) {print \$2; exit}}')
                awk -v third_largest_size="\$third_largest_size" '\$2 >= 6 && \$2 >= third_largest_size {print \$1}' \$dir/\${filename}_sorted_orf_clusters.txt > \$dir/\${filename}_largest_orf_clusters.txt

                clstr2txt.pl \$clstr_file > \$dir/\${filename}_clusterinfo.txt
                awk 'NR==FNR {headers[\$0]; next} \$2 in headers' \$dir/\${filename}_largest_orf_clusters.txt \
                \$dir/\${filename}_clusterinfo.txt | awk '\$5 == 1' | awk '{print \$1}' > \$dir/\${filename}_potential_ORFs.txt
                seqkit grep -f \$dir/\${filename}_potential_ORFs.txt \$dir/\${filename}_clustered_orfs.fasta > \$dir/\${filename}_potential_ORFs.fa
                fi
            fi
        fi
    fi
done
"""
}


process MOVEFORTRANSFER{

    executor 'local'

    input:
    path extendfolder
    output:
    path extendfolder    
    script:
    """   
    if [ ! -d "params.curatedtefolder" ]; then
        # If it doesn't exist, create it
        mkdir -p "${params.curatedtefolder}" 
        mkdir -p "${params.cronefolder}" 
        mkdir -p "${params.orffolder}" 
    fi

    find "${params.cronefolder}" -type f -print0 | xargs -0 rm -f
    find "${params.orffolder}" -type f -print0 | xargs -0 rm -f

    cp ${params.extendfolder}/*/*_potential_ORFs.fa ${params.orffolder}
    cp ${params.extendfolder}/*/*_cons.fa ${params.cronefolder}
    cp ${params.extendfolder}/*/*_reps.fa ${params.cronefolder}
    """
}


process RUNTEAID {
    executor 'local'

    input:
    path extendfolder

    output:
    path extendfolder

    script:
    """
    if [ ! -d "${params.teaidfolder}" ]; then
        mkdir -p "${params.teaidfolder}"
    fi

    for FILE in \$(ls ${params.cronefolder}/*_cons.fa)
    do
        /hpcfs/users/a1749756/tools/TE-Aid/TE-Aid -q ${params.cronefolder}/\$(basename \$FILE) -g ${params.genome} \
        -o ${params.teaidfolder}
    done
    """
}



workflow {
    CONDA
    checkfiles_ch = CHECKFILES(params.orf2folder)
    mafft_ch = MAFFT(checkfiles_ch)
    renaming_ch = RENAMING(mafft_ch)
    getorf_ch = GETORFS(renaming_ch)
    transfer_ch = MOVEFORTRANSFER(getorf_ch)
    teaid_ch = RUNTEAID(transfer_ch)

}


