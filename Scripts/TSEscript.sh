#!/bin/bash
yourfilenames=$(ls /home/fatemeh/Leishmania_2019/GenomeData/remain/*.fasta)
for FILE in $yourfilenames; do
	printf "$FILE"
        output="${FILE%_Genome*}"
        output1="${output##*TriTrypDB-41_}.tse.out"
        output2="${output##*TriTrypDB-41_}.SS.tse.out"
        output3="${output##*TriTrypDB-41_}.tse.fasta"
        tRNAscan-SE "${FILE}" -f "/home/fatemeh/Leishmania_2019/Leishmania_2019/Results/tRNAScanGenes/${output2}" -o "/home/fatemeh/Leishmania_2019/Leishmania_2019/Results/tRNAScanGenes/${output1}" -a "/home/fatemeh/Leishmania_2019/Leishmania_2019/Results/tRNAScanGenes/${output3}"
done

# I ended up installing tRNAscan with this command: conda install -c bioconda trnascan-se 

