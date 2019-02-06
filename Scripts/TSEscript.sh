#!/bin/bash -l
#SBATCH -t 16:00:00 ##### check runtime
#SBATCH --mem=128G ##### check memory requirements
#SBATCH --mail-type=ALL
#SBATCH --cwd
#SBATCH --mail-user=fhadinezhad@ucmerced.edu
#SBATCH -c 1 # number of cores


yourfilenames=$(ls /home/fhadinezhad/Leishmania_2019/GenomeData/*.fasta)
for FILE in $yourfilenames; do
        output="${FILE%_Genome*}"
        output1="${output##*TriTrypDB-41_}.tse.out"
        output2="${output##*TriTrypDB-41_}.SS.tse.out"
        output3="${output##*TriTrypDB-41_}.tse.fasta"
        ./tRNAscan-SE "${FILE}" -f "/home/fhadinezhad/Leishmania_2019/Leishmania_2019/Results/tRNAScanGenes/${output2}" -o "/home/fhadinezhadLeishmania_2019/Leishmania_2019/Results/tRNAScanGenes/${output1}" -a "/home/fhadinezhad/Leishmania_2019/Leishmania_2019/Results/tRNAScanGenes/${output3}"
done

