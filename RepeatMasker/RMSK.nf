params.species = "setmeup"
params.genomefolder = "${projectDir}/genomes/${params.species}/ncbi_dataset/data"
params.genomestring = file("$projectDir/genomes/${params.species}/ncbi_dataset/data/*/*.fna")
params.genome = "${params.genomestring[0]}"
params.rmskfolder = "$projectDir/analysis/${params.species}/Repeat_Masking"
params.repeats = "$projectDir/repeats/CGE_RVT_over200_80cluster.fa"
params.dbname = "${params.species}"
params.conda_env = '/hpcfs/users/a1749756/myconda/envs/bedtools'
params.clstr2txt = "${projectDir}/tools/clstr2txt.pl"
params.eukrepeats = "$projectDir/repeats/repeatmasker_libraries_for_CR1_project/new_combined_library_210825/tmp"


process CONDA{

    conda params.conda_env
    script:
    """
    source activate ${params.conda_env}
    """
}

process FOLDERSETUP {

    executor 'local'

    tag "set up ${params.species} Repeat_Masking folder"

    script:
    """
    # Check if the folder exists
    if [ ! -d "${projectDir}/analysis/${params.species}/Repeat_Masking2" ]; then
        # If it doesn't exist, create it
        mkdir -p "${projectDir}/analysis/${params.species}/Repeat_Masking2"
    fi
    """
}


process REPEATMASKER {
    cpus { 24 }
    time { '6h' }
    memory { '96.GB' }

    tag "Run RepeatMasker"

    input:
    path rmskfolder
    output:
    path rmskfolder

    script:
    """
    RepeatMasker -pa 24 \
      -lib ${params.eukrepeats} \
      -dir ${projectDir}/analysis/${params.species}/Repeat_Masking2 \
      ${params.genome} \
      > ${projectDir}/analysis/${params.species}/Repeat_Masking/repeatmasker.log 2>&1

    # Check if the expected .out file(s) exist
    if ls ${projectDir}/analysis/${params.species}/Repeat_Masking2/*.out 1> /dev/null 2>&1; then
        echo "RepeatMasker finished successfully, .out file created."
    else
        echo "ERROR: No .out file created in Repeat_Masking2/" >&2
        exit 1
    fi
    """
}

workflow {
    CONDA
    FOLDERSETUP()
    repeatmasker_ch = REPEATMASKER(params.genome)
}
