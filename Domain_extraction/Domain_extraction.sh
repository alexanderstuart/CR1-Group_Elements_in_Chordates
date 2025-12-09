#-----------------------------------------------------------------------------------------------------------------------------------
# define environment, fetch ORF
conda activate domain

data_source=test
mkdir ./data/parse/$data_source/
source_sequence="tmp/seq"
Database=./data/mmseqs2/Database/Protein/$data_source

mkdir -p $Database

getorf  -sequence tmp/seq \
        -min 200\
        -find 1 \
        -reverse N \
        -outseq $Database/$data_source.orf

cat $Database/$data_source.orf |  sed -e 's/ \[/_/g' -e 's/\]//g' -e 's/ - /_/g' > ./tmp/getorf

#-----------------------------------------------------------------------------------------------------------------------------------
# domain serach with hmsearch and Pfam

mkdir tmp/hmmer
rm tmp/hmmer/*
mkdir ./data/hmmer/$data_source/

orf=./tmp/getorf
db='../_general_data/hmm_datbase/Pfam_A_37/Pfam-A.hmm'

hmmsearch --cpu 20 \
          --noali \
          --domtblout ./data/hmmer/$data_source/$data_source.domtbl  \
          $db \
          $orf  

# parse hmmout
cat ./data/hmmer/$data_source/$data_source.domtbl  | sed 's/^#.*//g' \
                                                   | sed '/^[[:space:]]*$/d' \
                                                   | sed 's/ \+/ /g' \
                                                   | sed 's/ /\t/g' \
                                                   > ./data/hmmer/$data_source/$data_source.domtbl_parsed

                                                   
