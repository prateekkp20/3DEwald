#!/usr/bin/bash

make clean > output.txt
make >> output.txt

cd run/
echo "File,Reciprocal,,Real" > parallel1.csv
# echo "File,Self,,Reciprocal,,Real" > original.csv
# i=5
# echo -n "POSCAR."$i >> outputfiles.csv
# echo ",hey" >>outputfiles.csv
# ./coulomb.x > files2.txt
# echo "," >> files2.txt
# echo $i >> outputfiles.csv 
# echo $i",Iteration,Self,,Real,,Reciprocal" >> outputfiles.csv
# echo "Iteration,Self,,Real,,Reciprocal" >> outputfiles.csv

for i in {1..14} 
do
    sed -i "s/Posfile = big\/POSCAR\.$((i-1))/Posfile = big\/POSCAR\.$i/g" input.in
    echo -n "POSCAR"$i >> parallel1.csv
    for j in {1..10}
    do
        ./coulomb.x >> parallel1.csv
        echo " " >> parallel1.csv
    done
    echo " " >> parallel1.csv
    echo " " >> parallel1.csv
done
