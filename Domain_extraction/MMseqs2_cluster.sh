
# Multiple clustering with MMseqs2

data_source=test
species=test

seq=tmp/seq

start=80
stop=40
step=10
coverage=60
cov_mode=0
clus_mode=1

source_sequence=tmp/seq
Database=./data/mmseqs2/Database/Protein/$data_source/RVT/$species
Cluster=./data/mmseqs2/Cascade/Protein/$data_source/RVT/clusMode${clus_mode}_covMode${cov_mode}_Cov${coverage}/$species

mkdir -p $Database
mkdir -p $Cluster

rm $Database/*
rm $Cluster/*
rm -r ./tmp/mmseqs2/*

mmseqs createdb $source_sequence \
                $Database/$species

for (( i=$start; i>=$stop; i-=$step )); do
        # first clustering run
        mmseqs cluster $Database/$species \
                        $Cluster/$species'_'$i \
                        ./tmp/mmseqs2 \
                        --single-step-clustering 1 \
                        --min-seq-id 0.$i \
                        --cov-mode $cov_mode \
                        --cluster-mode $clus_mode \
                        --cluster-reassign 1 \
                        -c 0.$coverage \
                        -v 2
                        
        # create a tsv file with the clustering results
        mmseqs createtsv $Database/$species \
                        $Database/$species \
                        $Cluster/$species'_'$i \
                        $Cluster/$species'_'$i.tsv

        # count the number of sequences in each cluster
        cat $Cluster/$species'_'$i.tsv | cut -d$'\t' -f1 \
                                        | sort \
                                        | uniq -c \
                                        | sort -n \
                                        > $Cluster/$species'_'$i.count
        # cat $Cluster/$species'_'$i.count       
done                             
