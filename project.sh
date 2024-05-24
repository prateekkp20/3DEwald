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

# make clean > terminal.txt
# make >> terminal.txt

# cd run/
# echo "File,Reciprocal,," > parallel1.csv
# # echo "File,Self,,Reciprocal,,Real" > parallel1.csv

# for i in {1..14} 
# do
#     sed -i "s/Posfile = big\/POSCAR\.$((i-1))/Posfile = big\/POSCAR\.$i/g" input.in
#     echo -n "POSCAR"$i >> parallel1.csv
#     for j in {1..10}
#     do
#         ./coulomb.x >> parallel1.csv
#         echo " " >> parallel1.csv
#     done
#     echo " " >> parallel1.csv
#     echo " " >> parallel1.csv
# done



## code for multiple threads experiment but same cell

# cd run/
# echo "# of Threads,Self,,Reciprocal,,Real" > original_poscar15_fifty.csv
# cd ..
# for i in {2..24}
# do
#     sed -i "s/\#define NUM_THREADS $((i-1))/\#define NUM_THREADS $i/g" inc/const.h
#     make clean > output.txt
#     make >> output.txt

#     cd run/
#     echo -n $i >> original_poscar15_fifty.csv
#     for j in {1..5}
#     do
#         ./coulomb.x >> original_poscar15_fifty.csv
#         echo " " >> original_poscar15_fifty.csv
#     done
#     echo " " >> original_poscar15_fifty.csv
#     echo " " >> original_poscar15_fifty.csv
#     cd ..
# done


# # code for bench marking the calculations using lammps
# make clean > terminal.txt
# make >> terminal.txt

# cd run/
# # echo "File,Ecoul,Elong,Total" > output.csv

# for i in {1..5}
# do 
#     sed -i "s/Posfile = bench\/benchfour$((i-1))\.POSCAR/Posfile = bench\/benchfour$i\.POSCAR/g" input.in
#     # echo -n "POSCAR"$i >> output.csv
#     # for j in {1..10}
#     # do 
#         ./coulomb.x >> output1.csv
#         echo " " >> output1.csv
#     # done
#     # echo " " >> output1.csv
#     # echo " " >> output1.csv
# done

## code for multiple threads for different settings with same size and same number of atoms

# for cell in {1..10}
# do
    # sed -i "s/Posfile = same_same\/POSCAR\.$((cell-1))/Posfile = same_same\/POSCAR\.$cell/g" run/input.in
#     cd run/
#     echo "# of Threads,,Reciprocal,,Real" > same_poscar$cell\file.csv
#     cd ..
#     for i in {2..24}
#     do
#         sed -i "s/\#define NUM_THREADS $((i-1)) /\#define NUM_THREADS $i /g" inc/const.h
#         make clean > output.txt
#         make >> output.txt

#         cd run/
#         echo -n $i >> same_poscar$cell\file.csv
#         for j in {1..5}
#         do
#             ./coulomb.x >> same_poscar$cell\file.csv
#             echo " " >> same_poscar$cell\file.csv
#         done
#         echo " " >> same_poscar$cell\file.csv
#         echo " " >> same_poscar$cell\file.csv
#         cd ..
#     done
#     sed -i "s/\#define NUM_THREADS 24 /\#define NUM_THREADS 1 /g" inc/const.h
# done

## Testing and Validating the bspline code

# for order in {2..7}
# do 
#     sed -i "s/int n\=$((order-1))\;/int n\=$order\; /g" Bspline_fftw.cpp
#     for grid in {2..7}
#     do 
#         sed -i "s/vector\<int\> K\=\{$((grid-1)*10)\,60\,60\}\;/"
#     done
# done

## code for different threads experiment but same box, same grid and order
echo "# of Threads,Time for FFTW" > run/may22_exp5.csv
for i in {1..24}
do
    sed -i "s/\#define NUM_THREADS $((i-1)) /\#define NUM_THREADS $i /g" inc/const.h
    echo -n $i >> run/may22_exp5.csv
    for run in {1..5}
    do
        make clean > output.txt
        make >> output.txt
        cd run/
        ./coulomb.x >> may22_exp5.csv
        echo " " >> may22_exp5.csv
        cd ..
    done
    echo " " >> run/may22_exp5.csv
    echo " " >> run/may22_exp5.csv
done


## code for grids and order experiment but same box, same threads
# echo "# of Grids,Order of Interpolation,Time for FFTW,Relative Error" > run/may22_exp3.csv
# for grids in {1..8}
# do
#     for order in {3..10}
#     do
#     g=$((40+10*grids))
#     sed -i "s/\#define GRID_SIZE $((40+10*grids-10)) /\#define GRID_SIZE $((40+10*grids)) /g" inc/const.h
#     sed -i "s/\#define BSPLINE_ORDER $((order-1)) /\#define BSPLINE_ORDER $order /g" inc/const.h
#     echo -n "${g},${order}" >> run/may22_exp3.csv
#     for run in {1..5}
#         do 
#             # echo $grids
#             make clean > output.txt
#             make >> output.txt
#             cd run/
#             ./coulomb.x >> may22_exp3.csv
#             echo " " >> may22_exp3.csv
#             cd ..
#         done
#     echo " " >> run/may22_exp3.csv
#     echo " " >> run/may22_exp3.csv
#     done
# done