#!/bin/sh
rmax_arr=(1e-3 1e-4 1e-5 1e-6 1e-7)
N_arr=(100 1000 10000 100000 1000000)

make

for N in "${N_arr[@]}"
do
    ./main --filename facebook --file_name ./graphs/facebook.txt --method abwalk --pairs vertex-pairs --N ${N} >> facebook-parameter-pair-N.txt
    ./main --filename facebook --file_name ./graphs/facebook.txt --method localtree --pairs vertex-pairs --N ${N} >> facebook-parameter-pair-N.txt
    ./main --filename facebook --file_name ./graphs/facebook.txt --method bipush --pairs vertex-pairs --N ${N} >> facebook-parameter-pair-N.txt
done

for rmax in "${rmax_arr[@]}"
do
    ./main --filename facebook --file_name ./graphs/facebook.txt --method push --pairs vertex-pairs --rmax ${rmax} --N 10000 >> facebook-parameter-pair-rmax.txt
    ./main --filename facebook --file_name ./graphs/facebook.txt --method bipush --pairs vertex-pairs --rmax ${rmax} --N 10000 >> facebook-parameter-pair-rmax.txt
done

for N in "${N_arr[@]}"
do
    ./main --filename facebook --file_name ./graphs/facebook.txt --method source_landmark --pairs node --N ${N} >> facebook-parameter-source-N.txt
    ./main --filename facebook --file_name ./graphs/facebook.txt --method source --pairs node --num_trials 50 --N ${N} >> facebook-parameter-source-N.txt
done

for rmax in "${rmax_arr[@]}"
do
    ./main --filename facebook --file_name ./graphs/facebook.txt --method source-rmax --pairs node --num_trials 50 --N 10000 --rmax ${rmax} >> facebook-parameter-source-rmax.txt
done