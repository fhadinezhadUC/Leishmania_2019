#!/bin/bash

array=(A R N D C Q E G H I L K M X F P Z S T W Y V)
echo "Array size: ${#array[*]}"

for item in ${array[*]}
do
    printf "   %s\n" $item
    cat tsfm_finalinput_HomoC_EditedCovea.fasta | fasgrep "_${item}$|_ara${item}$" > "/home/fatemeh/Leishmania_2019/Leishmania_2019/Results/tsfminput_final/inputsfiles/TriTryp/${item}.fasta"
    # for HOMOC genes
    cat tsfm_finalinput_HomoC_EditedCovea.fasta | fasgrep "Homo${item}" > "/home/fatemeh/Leishmania_2019/Leishmania_2019/Results/tsfminput_final/inputsfiles/HomoC/Homo${item}.fasta"

    cat "/home/fatemeh/Leishmania_2019/Leishmania_2019/Results/tsfminput_final/inputsfiles/TriTryp/${item}.fasta" | fasconvert -o clustalw > "/home/fatemeh/Leishmania_2019/Leishmania_2019/Results/tsfminput_final/inputsfiles/TriTryp/clustalW/TryTryp_${item}.aln"
    cat "/home/fatemeh/Leishmania_2019/Leishmania_2019/Results/tsfminput_final/inputsfiles/HomoC/Homo${item}.fasta" | fasconvert -o clustalw > "/home/fatemeh/Leishmania_2019/Leishmania_2019/Results/tsfminput_final/inputsfiles/HomoC/clustalW/Homo_${item}.aln"

done
