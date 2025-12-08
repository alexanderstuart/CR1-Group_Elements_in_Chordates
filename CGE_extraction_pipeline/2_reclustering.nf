params.species = "setmeup"
params.divergence = "85"
params.genomefolder = "${projectDir}/genomes/${params.species}/ncbi_dataset/data"
params.genomestring = file("$projectDir/genomes/${params.species}/ncbi_dataset/data/*/*.fna")
params.genome = "${params.genomestring[0]}"
params.orf2folder = "$projectDir/analysis/${params.species}/ORF2"
params.divergencefolder = "${params.orf2folder}/${params.divergence}"
params.repeats = "$projectDir/repeats/CGE_RVT_over200_85cluster.fa"
params.dbname = "${params.species}"
params.conda_env = '/hpcfs/users/a1749756/myconda/envs/bedtools'
params.clstr2txt = "${projectDir}/tools/clstr2txt.pl"
params.eukrepeats = "$projectDir/repeats/RepBase2022_RMSK_Eukaryota_withTacu.fa"
params.allelements = "${params.divergencefolder}/${params.species}_sizefiltered_named_subfams${params.divergence}_parsed.fa"
params.rmskout = "${params.orf2folder}/RepeatMasker/${params.species}_sizefiltered_named_subfams${params.divergence}_parsed.fa.out"
params.cgelist = "${params.divergencefolder}/${params.species}_${params.divergence}_cgelist.txt"


process CONDA{

    conda params.conda_env
    script:
    """
    source activate ${params.conda_env}
    """
}

process RENAMING {
    input:
    path orf2folder

    output:
    path orf2folder
    
    script:
    """
    cat ${params.orf2folder}/${params.species}_sizefiltered.fa | awk -F '>' '/^>/ { \$1 = ">SEQXXXX" ++i + 1000000 "_" } 1 ' OFS="" | \
    sed 's/SEQXXXX1/SEQ/' > ${params.orf2folder}/${params.species}_sizefiltered_named.fa   
    
    """
}

process CLUSTERING {
    cpus { 12 }
    time { '24h' }
    memory { '60.GB' }

    input:
    path orf2folder

    output:
    path orf2folder
    
    script:
    """
    # Check if the divergence folder exists
    if [ ! -d "${params.orf2folder}/${params.divergence}" ]; then
        # If it oesn't exist, create it
        mkdir -p "${params.orf2folder}/${params.divergence}"
    fi

    
    cd-hit-est -i ${params.orf2folder}/${params.species}_sizefiltered_named.fa \
    -o ${params.divergencefolder}/${params.species}_sizefiltered_named_subfams${params.divergence} -c 0.${params.divergence} -d 0
    """
}

process PARSING {
    
    input:
    path orf2folder

    output:
    path orf2folder


    script:
    """
    
    perl ${params.clstr2txt} ${params.divergencefolder}/${params.species}_sizefiltered_named_subfams${params.divergence}.clstr \
    | awk -F '\\t' '\$3 > 10 && \$5 == 1 {print \$1}' > ${params.divergencefolder}/${params.species}_tmp
    seqkit grep -r -f ${params.divergencefolder}/${params.species}_tmp ${params.divergencefolder}/${params.species}_sizefiltered_named_subfams${params.divergence} \
    | awk -F '>' '/^>/ { \$1 = ">${params.species}_TEXXXX" ++i + 100000 "_" } 1 ' OFS="" \
    | sed 's/TEXXXX1/TE/' | sed '/^>/ s/\$/#LINE/' \
    > ${params.divergencefolder}/${params.species}_sizefiltered_named_subfams${params.divergence}_parsed_FULLNAME.fa
    cat ${params.divergencefolder}/${params.species}_sizefiltered_named_subfams${params.divergence}_parsed_FULLNAME.fa | sed -E 's/_SEQ.*/#LINE/' \
    > ${params.divergencefolder}/${params.species}_sizefiltered_named_subfams${params.divergence}_parsed.fa
    """
}

process REPEATMASKER {

    tag "Run RepeatMasker on ${params.species} CR1 elements"

    input:
    path orf2folder

    output:
    path orf2folder

    script:
    """
    # Check if the folder exists
    if [ ! -d "${params.orf2folder}/RepeatMasker" ]; then
        # If it oesn't exist, create it
        mkdir -p "${params.orf2folder}/RepeatMasker"
    fi

    # Run RepeatMasker with specified parameters
    RepeatMasker -pa 4 -nolow -lib ${params.eukrepeats} -dir ${params.orf2folder}/RepeatMasker ${params.allelements}
    """
}

process CGEPARSING {
    tag "extract CGE elements"

    input:
    path orf2folder

    output:
    path orf2folder 

    script:
    """
    # Extract CGE elements from RepeatMasker output
    grep -E  "/CR1|/L2|/Jockey|/L2A|/L2B|/Rex1|/Daphne|/Crack" ${params.rmskout} | awk '{ print \$5 }' | sort | uniq > ${params.cgelist}
    
    # Use seqkit to grep CGE elements from the genomic sequences
    seqkit grep -r -f ${params.cgelist} ${params.allelements} > ${params.divergencefolder}/${params.species}_${params.divergence}_cgeelements.fa
    
    grep -r -f ${params.cgelist} ${params.rmskout} > ${params.divergencefolder}/${params.species}_${params.divergence}_cgeRMSKlist.txt
    """
}


process CHECK_FILE {
    input:
    path orf2folder

    script:
    """
    if [ -s ${params.cgelist} ]
    then
        echo "${params.species}" >> ${projectDir}/other/completed${params.divergence}_runs.txt
    else
        echo "Error: File is empty" >&2
        exit 1
    fi
    """
}

workflow {
    CONDA
    rename_ch = RENAMING(params.orf2folder)
    cluster_ch = CLUSTERING(rename_ch)
    parsing_ch = PARSING(cluster_ch)
    repeatmasker_ch = REPEATMASKER(parsing_ch)
    cge_ch = CGEPARSING(repeatmasker_ch)
    checkfile_ch = CHECK_FILE(cge_ch)
    }
