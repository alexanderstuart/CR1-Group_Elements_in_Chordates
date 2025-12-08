params.species = "setmeup"
params.genomefolder = "${projectDir}/genomes/${params.species}/ncbi_dataset/data"
params.genomestring = file("$projectDir/genomes/${params.species}/ncbi_dataset/data/*/*.fna")
params.genome = "${params.genomestring[0]}"
params.orf2folder = "$projectDir/analysis/${params.species}/ORF2"
params.repeats = "$projectDir/repeats/CGE_RVT_over200_80cluster.fa"
params.dbname = "${params.species}"
params.conda_env = '/hpcfs/users/a1749756/myconda/envs/bedtools'
params.clstr2txt = "${projectDir}/tools/clstr2txt.pl"
params.eukrepeats = "$projectDir/repeats/RepBase2022_RMSK_Eukaryota_withTacu.fa"
params.allelements = "${params.orf2folder}/${params.species}_sizefiltered_subfams80_parsed.fa"
params.rmskout = "${params.allelements}.out"


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


process REPEATMASKER {

    tag "Run RepeatMasker"

    input:
    path allelements
    path eukrepeats

    output:
    path orf2folder

    script:
    """
    RepeatMasker -pa 4 -nolow -lib ${params.eukrepeats} -dir ${params.orf2folder} ${params.allelements}
    """
}

process PREPAREBEDFILE {
    tag "extract CGE elements"

    input:
    path orf2folder

    output:
    path orf2folder

    script:
    """
    grep -r  "/CR1|/L2|/Jockey|/L2A|/L2B|/" *.out

    awk '{ if (\$2 > \$3) { t = \$2; \$2 = \$3; \$3 = t } print }' ${params.orf2folder}/${params.species}_raw.txt \
    | awk -F '\\t' '{ if ((\$3 - \$2) > 400) print }' \
    | awk -F '\\t' '{ if ((\$3 - \$2) < 850) print }' \
    | sort -k1,1 -k2,2n | bedtools merge \
    | awk -F '\\t' '{ if ((\$3 - \$2) < 850) print }' \
    > ${params.orf2folder}/${params.species}_sizefiltered.bed
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
