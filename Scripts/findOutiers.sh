#!/bin/bash

declare -a model=(A R N D C Q E G H I L K M X F P S T W Y V X)

pathstring="/home/fatemeh/Leishmania_2019/Leishmania_2019/Results/Integrated_Genes/TriTrypgenemodel_Intersection/TriTryp_EditedCovea/"
pathtoOD="/home/fatemeh/Leishmania_2019/Leishmania_2019/Scripts/OD-Seq/OD-seq"
filenames=$(ls "${pathstring}"*.fasta)
mkdir -p "${pathstring}/outlier"
for ref in $filenames; do
    printf " %s\n" $ref3
    ref2=${ref#"${pathstring}"}
    ref3=${ref2%".fasta"}
    "${pathtoOD}" -i "${ref}" -o "${pathstring}/outlier/${ref3}.txt"
done


declare -A matrix
num_rows=2
num_columns=22

for ((i=1;i<=num_rows;i++)) do
    for ((j=1;j<=num_columns;j++)) do
	countoutlier=$(cat "${pathstring}outlier/${model[j-1]}.txt" | wc -l)
	totaloutlier=$(cat "${pathstring}${model[j-1]}.fasta" | wc -l)
        matrix[1,$j]=$(($countoutlier / 2)) # outliers
	matrix[2,$j]=$(($totaloutlier / 2))
    done
done

sum=0
printf "model %3s outlier#/total" > "${pathstring}stats.txt"
for ((j=1;j<=num_columns;j++)) do
    printf "\n ${model[j-1]} %9s ${matrix[1,$j]}/${matrix[2,$j]}" >> "${pathstring}stats.txt"
    sum=$((${matrix[1,$j]} + $sum))
done
printf "\n total number of outliers = ${sum}" >> "${pathstring}stats.txt"

