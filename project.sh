#!/usr/bin/bash

# # echo "File,Reciprocal,,Real" > parallel1.csv
# i=5
# echo -n "POSCAR."$i >> outputfiles.csv
# echo ",hey" >>outputfiles.csv
# ./coulomb.x > files2.txt
# echo "," >> files2.txt
# echo $i >> outputfiles.csv 
# echo $i",Iteration,Self,,Real,,Reciprocal" >> outputfiles.csv
# echo "Iteration,Self,,Real,,Reciprocal" >> outputfiles.csv


## code for multiple cells and same threads

make clean > terminal.txt
make >> terminal.txt

cd run/
echo "File,Reciprocal,," > parallel1.csv
# echo "File,Self,,Reciprocal,,Real" > parallel1.csv

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



## code for multiple different threads experiment

# cd run/
# echo "# of Threads,Self,,Reciprocal,,Real" > original.csv
# cd ..
# for i in {2..15}
# do
#     sed -i "s/\#define NUM_THREADS $((i-1))/\#define NUM_THREADS $i/g" inc/const.h
#     make clean > output.txt
#     make >> output.txt

#     cd run/
#     echo -n $i >> original.csv
#     for j in {1..10}
#     do
#         ./coulomb.x >> original.csv
#         echo " " >> original.csv
#     done
#     echo " " >> original.csv
#     echo " " >> original.csv
#     cd ..
# done


# code for bench marking the calculations using lammps

# make clean > terminal.txt
# make >> terminal.txt

# cd run/
# echo "File,Ecoul,Elong,Total" > output.csv

# for i in {1..14}
# do 
#     sed -i "s/Posfile = big\/POSCAR\.$((i-1))/Posfile = big\/POSCAR\.$i/g" input.in
#     echo -n "POSCAR"$i >> output.csv
#     for j in {1..10}
#     do 
#         ./coulomb.x >> output.csv
#         echo " " >> output.csv
#     done
#     echo " " >> output.csv
#     echo " " >> output.csv
# done
