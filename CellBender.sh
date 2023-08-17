#!/bin/bash

bsub -R 'rusage[mem=40GB]' -Is bash

module load python/3.7.0 anaconda3/2019.07

#conda init bash
conda activate cellbender

#dir=$1
#cells=$2
#droplets=$3

droplets=5000

# MFCON007dcM
dir=/icgc/dkfzlsdf/analysis/B210/data/mf/cellranger201_count_24192-25608_4839STDY7131583_GRCh38/
cells=1828

# MFCON007dfM # bad one
dir=/icgc/dkfzlsdf/analysis/B210/data/mf/cellranger201_count_24192-25608_4839STDY7131582_GRCh38/
cells=2307

# MFCON020afM # bad one
dir=/icgc/dkfzlsdf/analysis/B210/data/mf/cellranger201_count_24192-25608_4839STDY7131581_GRCh38/
cells=1795

# MFCON020acM
dir=/icgc/dkfzlsdf/analysis/B210/data/mf/cellranger201_count_24192-25608_4839STDY7131584_GRCh38/
cells=162

# MFCON007efM
dir=/icgc/dkfzlsdf/analysis/B210/data/mf/cellranger201_count_23156_6_GRCh38/
cells=424

# MFCON010dfM
dir=/icgc/dkfzlsdf/analysis/B210/data/mf/cellranger201_count_23156_7_GRCh38/
cells=561

# MFCON018bfM
dir=/icgc/dkfzlsdf/analysis/B210/data/mf/cellranger201_count_23156_8_GRCh38
cells=111

cellbender remove-background \
--input $dir/outs/raw_gene_bc_matrices/GRCh38/ \
--output $dir/outs/raw_gene_bc_matrices/GRCh38/cellbender_matrix.h5 \
--expected-cells $cells \
--total-droplets-included $droplets


# test code:
#echo "${BASH_VERSION}"
#conda activate cellbender
#printf '#!/bin/bash\nconda init bash; conda activate cellbender; /bin/bash --version' > test.sh
#chmod +x test.sh; ./test.sh
