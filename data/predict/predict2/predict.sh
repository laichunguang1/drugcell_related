#!/bin/bash
index=2
inputdir="/home/cxl1214/drugcell/drugcell1.0/alldata/"
inputdir2="/home/cxl1214/drugcell/drugcell1.0/alldata/split52/"
gene2idfile=$inputdir"gene2ind.txt"
cell2idfile=$inputdir"cell2ind.txt"
drug2idfile=$inputdir"drug2ind.txt"
testdatafile=$inputdir2"drugcell_test"$index".txt"

mutationfile=$inputdir"cell2mutation.txt"
drugfile=$inputdir"drug2fingerprint.txt"

modelfile="/home/cxl1214/drugcell/test/model/"$index"split/Model_sample/model_292.pt"

resultdir="Result_sample"
hiddendir="Hidden_sample"

cudaid=$1

if [$cudaid = ""]; then
	cudaid=0
fi

mkdir $resultdir
mkdir $hiddendir

source activate pytorch3drugcell

python -u /home/cxl1214/drugcell/drugcell1.0/code/predict_drugcell.py -gene2id $gene2idfile -cell2id $cell2idfile -drug2id $drug2idfile -genotype $mutationfile -fingerprint $drugfile -hidden $hiddendir -result $resultdir -predict $testdatafile -load $modelfile -cuda $cudaid > test_sample.log
