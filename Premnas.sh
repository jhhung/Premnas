#!/bin/bash


echo "------------------"
echo "Check docker image"
echo "------------------"

PREMNAS_DOCKER="toeric/premnas:1.0.6"
mode=false
susceptibility_threshold=0.9
consistency_threshold=0.8

if [[ "$(docker images $PREMNAS_DOCKER -q  2> /dev/null)" == "" ]]; then
	echo "Premnas images doesn't exit"
	echo "Downloading:"
	docker pull $PREMNAS_DOCKER
else 
	echo "Premnas images exist!"
fi


if [[ "$(docker images cibersortx/fractions:latest -q  2> /dev/null)" == "" ]]; then
	echo "CIBERSORTx images doesn't exit"
	echo "Downloading:"
	docker pull cibersortx/fractions
else 
	echo "CIBERSORTx images exist!"
fi


while getopts I:O:D:S:U:T:C:M:m:s:c: flag
do
	case "${flag}" in
		# required
		I) input_dir=${OPTARG};;
		O) output_dir=${OPTARG};;
		D) single_cell=${OPTARG};;
		S) sc_source=${OPTARG};;
		U) user_name=${OPTARG};;
		T) token=${OPTARG};;
		C) mixture=${OPTARG};;
		M) meta_data=${OPTARG};;

		# optional
		m) mode=${OPTARG};;
		s) susceptibility_threshold=${OPTARG};;
		c) consistency_threshold=${OPTARG};;
	esac
done

echo "--------------------------------------------"
echo "Learing ad hoc subpopulation characteristics"
echo "--------------------------------------------"


echo "Single cell data: $input_dir$single_cell"
echo "Single cell source: $input_dir$sc_source"
echo "Output folder: $output_dir"

docker run  -ti --rm \
	-v $input_dir:/input_dir \
	-v $output_dir:/output_dir \
	$PREMNAS_DOCKER \
	/Premnas.sh $single_cell $sc_source

if [[ "$mode" == "true" ]]; then
	echo "DONE!"
	exit
fi

echo "----------------------------"
echo "Performing digital cytometry"
echo "----------------------------"

exit
cp $output_dir/subpopulation-characteristic.txt $input_dir/

docker run -v $input_dir:/src/data \
	-v $output_dir:/src/outdir cibersortx/fractions \
	--username $user_name \
	--token $token \
	--single_cell TRUE  --fraction 0 --rmbatchSmode TRUE --perm 500 \
	--refsample subpopulation-characteristic.txt \
	--mixture $mixture

echo "------------------------------"
echo "Analyzing subpopulation change"
echo "------------------------------"

python3 src/Treatment-election.py $meta_data $output_dir/CIBERSORTx_Adjusted.txt $susceptibility_threshold $consistency_threshold

#remove tmp file
echo "Selected cocktail therapy is stored in Treatment-selection-output.csv"
echo "---------"
echo "ALL DONE!"
echo "---------"



