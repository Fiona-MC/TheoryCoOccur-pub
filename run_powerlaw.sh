#!/bin/bash
# this runs all of the power law analyses 

# Initialize conda for bash shell
source $(conda info --base)/etc/profile.d/conda.sh

SAVE_DIR="./degree_sequences_np"
# mkdir $SAVE_DIR
# mkdir $SAVE_DIR/cooc
# mkdir $SAVE_DIR/inter

# # run for 1000 species
echo 'running python process_data_degSeq.py 0.1 0.001 1000'
python process_data_degSeq.py 0.1 0.001 1000 $SAVE_DIR/
Rscript powerLaw.r 0.1 0.001 1000

echo 'running python process_data_degSeq.py 0.2 0.001 1000'
python process_data_degSeq.py 0.2 0.001 1000 $SAVE_DIR/
Rscript powerLaw.r 0.2 0.001 1000

echo 'running python process_data_degSeq.py 0.3333 0.001 1000'
python process_data_degSeq.py 0.3333 0.001 1000 $SAVE_DIR/
Rscript powerLaw.r 0.3333 0.001 1000

# # run for 100 species
echo 'running python process_data_degSeq.py 0.3333 0.1 100'
python process_data_degSeq.py 0.3333 0.1 100 $SAVE_DIR/
Rscript powerLaw.r 0.3333 0.1 100

echo 'running python process_data_degSeq.py 0.2 0.5 100'
# run for 100 species with 1/5 interactFactor and 0.5 pPositive
python process_data_degSeq.py 0.2 0.5 100 $SAVE_DIR/
Rscript powerLaw.r 0.2 0.5 100

echo 'running python process_data_degSeq.py 0.1 0.5 100'
python process_data_degSeq.py 0.1 0.5 100 $SAVE_DIR/
Rscript powerLaw.r 0.1 0.5 100
