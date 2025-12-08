params.species = "setmeup"
params.divergence = "85"
params.genomefolder = "${projectDir}/genomes/${params.species}/ncbi_dataset/data"
params.genomestring = file("${projectDir}/genomes/${params.species}/ncbi_dataset/data/*/*.fna")
params.genome = "${params.genomestring[0]}"
params.orf2folder = "${projectDir}/analysis/${params.species}/ORF2"
params.repeatfolder = "${projectDir}/repeats"
params.repeatcuration = "${projectDir}/analysis/${params.species}/repeatcuration"
params.mchelperfolder = "${projectDir}/analysis/${params.species}/repeatcuration/MCHelper"
params.eukrepeats = "${projectDir}/repeats/RepBase2022_RMSK_Eukaryota_withTacu.fa"
params.hmm = "${projectDir}/repeats/dummyfile.hmm"
params.dbname = "${params.species}"
params.conda_env = '/hpcfs/users/a1749756/myconda/envs/MCHelper'


process CONDA{

    tag "Activate ${params.conda_env} environment"

    conda params.conda_env
    script:
    """
    source activate ${params.conda_env}
    """
}

process FOLDERSETUP {

    executor 'local'

    tag "makin folders :)"

    script:
    """
    # Check if the folder exists
    if [ ! -d "${projectDir}/analysis/${params.species}/repeatcuration/MCHelper" ]; then
        # If it doesn't exist, create it
        mkdir -p "${projectDir}/analysis/${params.species}/repeatcuration/MCHelper"
    fi
    """
}

process MCHELPER {
    tag "Run MCHelper on ${params.species} repeats"

    cpus { 30 }
    time { '72h' }
    memory { '120.GB' }


    input:
    path orf2folder
    path repeatfolder
    path mchelperfolder

    output:
    path mchelperfolder
    
    script:
    """    
    python3 /hpcfs/users/a1749756/tools/mchelper/MCHelper.py -l ${params.orf2folder}/${params.divergence}/${params.species}_${params.divergence}_CGEelements.fa \
    -o ${params.mchelperfolder} -g ${params.genome} --input_type fasta -a F -r A \
    -t 6 -m ${params.eukrepeats} -b ${params.hmm} -x 8 -e 300 -c 4
    """
}

process CHECK_FILE {
    tag "Check if file is empty"

    input:
    val mchelper_output
    path mchelperfolder

    script:
    """
    if [ -s ${params.mchelperfolder}/curated_sequences_NR.fa ]
    then
        echo "${params.species}" >> ${projectDir}/other/completed_MCHelper_runs.txt
    else
        echo "Error: File is empty" >&2
        exit 1
    fi
    """
}

workflow {
    CONDA
    foldersetup_ch = FOLDERSETUP()
    mchelper_ch = MCHELPER(foldersetup_ch, params.orf2folder, params.repeatfolder, params.mchelperfolder)
    CHECK_FILE(mchelper_ch, params.mchelperfolder)
}


