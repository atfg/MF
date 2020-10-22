#bsub -R "rusage[mem=40GB]" ""


module load python/3.7.0 anaconda3/2019.07
conda activate cellbender

cellbender remove-background \
--input $PATH/outs/raw_feature_bc_matrix/GRCh38/ \
--output $PATH/outs/raw_feature_bc_matrix/GRCh38/cellbender_matrix.h5 \
--expected-cells $CELLS \
--total-droplets-included $DROPLETS

