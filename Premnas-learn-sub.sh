#!/bin/bash

if [ $# -ne 2 ]; then
	echo "Please run the shell script as: sh Premnas.sh single_cell_data.txt metadata.txt"
	exit
else
	echo "---------- Begin to run Premnas ----------"
fi

mkdir output_dir

echo "Processing single cell data:"
Rscript /src/preprocessing_sc.R $1 Project_name $2 > /dev/null 2 >&1

if [ ! -f "/output_dir/seurat_onlyfiltration.RDS" ]; then
	echo "ERROR: Failed processing"
	exit
fi

echo "Processing reduce and batch correction" 
Rscript /src/preRunACTIONet.R output_dir seurat_onlyfiltration.RDS > /dev/null 2 >&1

if [ ! -f "/output_dir/sc-after-reduce-and-batchcorrect.RDS" ]; then
	echo "ERROR: Failed reduce/batch correction"
	exit
fi

echo "run ACTIONet"
Rscript /src/runACTIONet.R output_dir sc-after-reduce-and-batchcorrect.RDS
if [ ! -f "/output_dir/ACTIONet-model" ]; then
	echo "ERROR: Failed to runACTIONet"
	exit
fi

echo "Write ACTIONet output"
Rscript /src/runActioNet_K.R output_dir > /dev/null 2 >&1

echo "Pruning out cells by archetypal explicit function (threshold = 0.6)"
python3 /src/Prune_cell.py

echo "Construct CIBERSORTx intput"
Rscript /src/prepareForCibersortx.R output_dir seurat_onlyfiltration.RDS pruned-assigned-subpopulation.txt > /dev/null 2 >&1


echo "---------- Done ----------"
