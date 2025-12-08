params.species = "setmeup"
params.genomefolder = "${projectDir}/genomes/${params.species}/ncbi_dataset/data"
params.genomestring = file("$projectDir/genomes/${params.species}/ncbi_dataset/data/*/*.fna")
params.genome = "${params.genomestring[0]}"
params.orf2folder = "$projectDir/analysis/${params.species}/ORF2"
params.repeats = "$projectDir/repeats/CGE_RVT_over200_80cluster.fa"
params.dbname = "${params.species}"
params.conda_env = '/hpcfs/users/a1749756/myconda/envs/bedtools'
params.clstr2txt = "${projectDir}/tools/clstr2txt.pl"



process CONDA{

    conda params.conda_env
    script:
    """
    source activate ${params.conda_env}
    """
}

process FOLDERSETUP {

    tag "set up folders"

    script:
    """
    # Check if the folder exists
    if [ ! -d "${projectDir}/analysis/${params.species}/ORF2" ]; then
        # If it oesn't exist, create it
        mkdir -p "${projectDir}/analysis/${params.species}/ORF2"
    fi
    """
}

process MAKEBLASTDB {
    tag "makeblastdb on ${params.species}"

    //queue 'icelake'

    //publishDir params.genomefolder, mode:'copy'

    input:
    path genomefolder

    output:
    path genomefolder
    
    script:
    """
    makeblastdb -in ${params.genome} -dbtype nucl
    """
}

process TBLASTN {

    cpus { 8 }
    time { '72h' }
    memory { '120.GB' }


    tag "blast CGE elements against ${params.species}"

    input:
    path genome
    path orf2folder

    output:
    path orf2folder
    
    script:
    """
    tblastn -query ${params.repeats} -db ${params.genome} -out ${params.orf2folder}/${params.species}_raw.txt \
    -outfmt "6 sacc sstart send" -evalue 0.0001 
    """
}

process PREPAREBEDFILE {
    tag "Orient and size filter blast hits"

    input:
    path orf2folder

    output:
    path orf2folder

    script:
    """
    awk '{ if (\$2 > \$3) { t = \$2; \$2 = \$3; \$3 = t } print }' ${params.orf2folder}/${params.species}_raw.txt \
    | awk -F '\\t' '{ if ((\$3 - \$2) > 400) print }' \
    | awk -F '\\t' '{ if ((\$3 - \$2) < 850) print }' \
    | sort -k1,1 -k2,2n | bedtools merge \
    | awk -F '\\t' '{ if ((\$3 - \$2) < 850) print }' \
    > ${params.orf2folder}/${params.species}_sizefiltered.bed
    """
}

process EXTRACTSSEQS {
    tag "Extract sequences from genome"

    input:
    path genome
    path orf2folder

    output:
    path orf2folder
    
    script:
    """
    bedtools getfasta -fi ${params.genome} \
    -bed ${params.orf2folder}/${params.species}_sizefiltered.bed \
     -fo ${params.orf2folder}/${params.species}_sizefiltered.fa 

    """
}

process CLUSTERING {
    tag "Cluster at 80%"

    cpus { 8 }
    time { '48h' }
    memory { '64.GB' }

    input:
    path orf2folder

    output:
    path orf2folder
    
    script:
    """
    cd-hit-est -i ${params.orf2folder}/${params.species}_sizefiltered.fa \
    -o ${params.orf2folder}/${params.species}_sizefiltered_subfams80 -c 0.8 -d 0
    """
}

process PARSING {
    tag "Remove clusters with less than 5 representatives"

    input:
    path orf2folder

    output:
    path orf2folder
    
    script:
    """    
    perl ${params.clstr2txt} ${params.orf2folder}/${params.species}_sizefiltered_subfams80.clstr| awk -F '\\t' '\$3 > 5 && \$5 == 1 {print \$1}' > ${params.orf2folder}/${params.species}_tmp
    seqkit grep -r -f ${params.orf2folder}/${params.species}_tmp ${params.orf2folder}/${params.species}_sizefiltered_subfams80 | \
    awk -F '>' '/^>/ { \$1 = ">${params.species}_TEXXXX" ++i + 100000 "_" } 1 ' OFS="" | \
    sed 's/TEXXXX1/TE/' | sed '/^>/ s/\$/#LINE/' \
    > ${params.orf2folder}/${params.species}_sizefiltered_subfams80_parsed.fa
    """
}

process CHECK_FILE {
    tag "Check if file is empty"

    input:
    val parsing_output
    path orf2folder

    script:
    """
    if [ -s ${params.orf2folder}/${params.species}_sizefiltered_subfams80_parsed.fa ]
    then
        echo "${params.species}" >> ${projectDir}/other/completed_runs.txt
    else
        echo "Error: File is empty" >&2
        exit 1
    fi
    """
}

workflow {
    CONDA
    FOLDERSETUP()
    makeblastdb_ch = MAKEBLASTDB(params.genome)
    tblastn_ch = TBLASTN(makeblastdb_ch, params.repeats)
    preparebed_ch = PREPAREBEDFILE(tblastn_ch)
    extractseqs_ch = EXTRACTSSEQS(preparebed_ch, params.genome)
    cluster_ch = CLUSTERING(extractseqs_ch)
    parsing_ch = PARSING(cluster_ch)
    CHECK_FILE(parsing_ch, params.orf2folder)
}
