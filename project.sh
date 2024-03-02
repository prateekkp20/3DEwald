#!/usr/bin/bash

make clean > output.txt
make >> output.txt

cd run/
# i=5
echo "File,Self,,Real,,Reciprocal" > outputfiles1.csv
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
    echo -n "POSCAR"$i >> outputfiles1.csv
    for j in {1..5}
    do
        ./coulomb.x >> outputfiles1.csv
        echo " " >> outputfiles1.csv
    done
    echo " " >> outputfiles1.csv
    echo " " >> outputfiles1.csv
done


