conda activate domain

data_source=CR1_project
input_type=RVT

group=CR1_BEAST_aligned
# group=CR1_project

Database=./data/mmseqs2/Database/Protein/$data_source/$input_type/Group_$group
Search=./data/mmseqs2/Search/Protein/$data_source/$input_type/Group_$group

seq=tmp/seq

# Get pairwise distance

mkdir -p $Database
mkdir -p $Search

rm -r $Database/*
rm -r $Search/*

mmseqs createdb $seq \
                $Database/Group_$group

mmseqs search $Database/Group_$group \
              $Database/Group_$group \
              $Search/Group_$group \
              tmp/mmseqs2 \
              -a

mmseqs convertalis $Database/Group_$group \
                   $Database/Group_$group \
                   $Search/Group_$group \
                   $Search/Group_$group.m8

len_threshold=200
cat $Search/Group_$group.m8 | awk '$1 != $2' \
                            | awk -v th=$len_threshold '$4 >= th' \
                            | awk '{ if ($1 == $2) next; print $1, $2, $12, $3}' \
                            | awk '{print ($1 < $2 ? $1 FS $2 FS $3 FS $4: $2 FS $1 FS $3 FS $4)}' \
                            | sort -k1,1 -k4,4nr | awk '!seen[$1 FS $2]++' \
                            > tmp/all_unique_edges

# save top edges
cat tmp/all_unique_edges | awk '!seen[$1]++' > tmp/top_edges
