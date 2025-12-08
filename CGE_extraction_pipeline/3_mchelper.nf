params.species = "setmeup"
params.divergence = "85"
params.genomefolder = "/hpcfs/users/a1749756/CGE_project/genomes/${params.species}/ncbi_dataset/data"
params.genomestring = file("/hpcfs/users/a1749756/CGE_project/genomes/${params.species}/ncbi_dataset/data/*/*.fna")
params.genome = "${params.genomestring[0]}"
params.orf2folder = "/hpcfs/users/a1749756/CGE_project/analysis/${params.species}/ORF2"
params.repeatfolder = "/hpcfs/users/a1749756/CGE_project/repeats"
params.repeatcuration = "/hpcfs/users/a1749756/CGE_project/analysis/${params.species}/repeatcuration"
params.mchelperfolder = "/hpcfs/users/a1749756/CGE_project/analysis/${params.species}/repeatcuration/MCHelper3"
params.eukrepeats = "/hpcfs/users/a1749756/CGE_project/repeats/RepBase2022_RMSK_Eukaryota_withTacu.fa"
params.hmm = "/hpcfs/users/a1749756/CGE_project/repeats/dummyfile.hmm"
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
    if [ ! -d "/hpcfs/users/a1749756/CGE_project/analysis/${params.species}/repeatcuration/MCHelper3" ]; then
        # If it doesn't exist, create it
        mkdir -p "/hpcfs/users/a1749756/CGE_project/analysis/${params.species}/repeatcuration/MCHelper3"
    fi
    """
}

process MCHELPER3 {
    tag "Run MCHelper3 on ${params.species} repeats"

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
        echo "${params.species}" >> /hpcfs/users/a1749756/CGE_project/other/completed_MCHelper3_runs.txt
    else
        echo "Error: File is empty" >&2
        exit 1
    fi
    """
}

workflow {
    CONDA
    foldersetup_ch = FOLDERSETUP()
    mchelper_ch = MCHELPER3(foldersetup_ch, params.orf2folder, params.repeatfolder, params.mchelperfolder)
    CHECK_FILE(mchelper_ch, params.mchelperfolder)
}


