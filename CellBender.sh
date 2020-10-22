#!/bin/bash

module load python/3.7.0 anaconda3/2019.07
conda activate cellbender

dir=$1
cells=$2
droplets=$3

cellbender remove-background \
--input $dir/outs/raw_feature_bc_matrix/GRCh38/ \
--output $dir/outs/raw_feature_bc_matrix/GRCh38/cellbender_matrix.h5 \
--expected-cells $cells \
--total-droplets-included $droplets

