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
        output="${output##*TriTrypDB-41_}.ara.out"
        ./aragorn -i116 -t -br -seq -w -e -l -d -o "${output}" "${FILE}"
done









