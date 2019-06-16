#!/bin/bash
# this script will run tsfm on all the genomes 
# for each genomes it will make a folder in each folder will make three folders: KLD,ID,Func, bubble_Table


folderpath="/home/fatemeh/Leishmania_2019/Leishmania_2019/Results/tsfminput_final/input3/"
tsfmpath="/home/fatemeh/Leishmania_Aug2018/phyloclassification/tsfm-master/tsfm"
folders=$(ls -l $folderpath | grep "^d" | awk -F" " '{print $9}')

# creating function logos


# creating ID, KLD, and tables for the bubble plots
for name in $folders; do 
if [ $name != "Logos" ]
then
printf "***************$name***************"
python3 "$tsfmpath/tsfm.py" -c "$folderpath/tsfm_finalinput_HomoC_structfile.txt" --idlogo "$folderpath/$name/clustalW/$name" "$folderpath/HOMO/clustalW/HOMO"
mv -- *.eps "$folderpath/Logos/$name/ID"

python3 "$tsfmpath/tsfm.py" -c "$folderpath/tsfm_finalinput_HomoC_structfile.txt" --kldlogo "$folderpath/$name/clustalW/$name" "$folderpath/HOMO/clustalW/HOMO"
mv -- *.eps "$folderpath/Logos/$name/KLD"

python3 "$tsfmpath/tsfm.py" -c "$folderpath/tsfm_finalinput_HomoC_structfile.txt" --bt "$folderpath/$name/clustalW/$name" "$folderpath/HOMO/clustalW/HOMO"
rm -- *.eps

mv *_Table.txt "$folderpath/Logos/$name/Bubble"
fi
done



