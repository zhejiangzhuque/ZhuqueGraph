#g++ -march=core2 -ffast-math -use_fast_math -pthread -std=c++11 -DSFMT_MEXP=607 -I SFMT-src-1.4.1/ -O3 -o TopPPR SFMT.c main.cpp 

rm youtube_u-result.txt

batch_size_arr=(1 2 3 4 5)

eps_arr=(0.5 0.4 0.3 0.2 0.1)

./TopPPR -d youtube_u -algo GEN_QUERY -n 50

./TopPPR -d youtube_u -algo GEN_GROUND_TRUTH -n 50 -a 0.2

./TopPPR -d youtube_u -algo COMPARE -a 0.2 -g undirected >> youtube_u-result.txt

for eps in "${eps_arr[@]}"
do
    ./TopPPR -d youtube_u -algo SS -n 50 -method PF -a 0.2 -epsilon ${eps} >> youtube_u-result.txt
    ./TopPPR -d youtube_u -algo SS -n 50 -method PW -a 0.2 -epsilon ${eps} >> youtube_u-result.txt
    ./TopPPR -d youtube_u -algo SS -n 50 -method PPF -a 0.2 -epsilon ${eps} -batch_size 2 >> youtube_u-result.txt
    ./TopPPR -d youtube_u -algo SS -n 50 -method PPW -a 0.2 -epsilon ${eps} -batch_size 2 >> youtube_u-result.txt
done

for eps in "${eps_arr[@]}"
do
    ./TopPPR -d youtube_u -algo SS -n 50 -method PF -a 0.2 -epsilon ${eps} -g undirected >> youtube_u-result.txt
    ./TopPPR -d youtube_u -algo SS -n 50 -method PPF -a 0.2 -epsilon ${eps} -batch_size 3 -g undirected >> youtube_u-result.txt
done

for eps in "${eps_arr[@]}"
do
    ./TopPPR -d youtube_u -algo PC -method MCW -a 0.2 -epsilon ${eps} >> youtube_u-result.txt
    ./TopPPR -d youtube_u -algo PC -method MCF -a 0.2 -epsilon ${eps} >> youtube_u-result.txt
    ./TopPPR -d youtube_u -algo PC -n 50 -method FORA -a 0.2 -epsilon ${eps} >> youtube_u-result.txt
    ./TopPPR -d youtube_u -algo PC -n 50 -method PF -a 0.2 -epsilon ${eps} >> youtube_u-result.txt
    ./TopPPR -d youtube_u -algo PC -n 50 -method PW -a 0.2 -epsilon ${eps} >> youtube_u-result.txt
    ./TopPPR -d youtube_u -algo PC -n 50 -method PPF -a 0.2 -epsilon ${eps} -batch_size 2 >> youtube_u-result.txt
    ./TopPPR -d youtube_u -algo PC -n 50 -method PPW -a 0.2 -epsilon ${eps} -batch_size 2 >> youtube_u-result.txt
done

for eps in "${eps_arr[@]}"
do
    ./TopPPR -d youtube_u -algo PC -method MCF -a 0.2 -epsilon ${eps} -g undirected >> youtube_u-result.txt
    ./TopPPR -d youtube_u -algo PC -n 50 -method PF -a 0.2 -epsilon ${eps} -g undirected >> youtube_u-result.txt
    ./TopPPR -d youtube_u -algo PC -n 50 -method PPF -a 0.2 -epsilon ${eps} -batch_size 2 -g undirected >> youtube_u-result.txt
done

./TopPPR -d youtube_u -algo GEN_GROUND_TRUTH -n 50 -a 0.01

./TopPPR -d youtube_u -algo COMPARE -a 0.01 -g undirected >> youtube_u-result.txt

for eps in "${eps_arr[@]}"
do
    ./TopPPR -d youtube_u -algo SS -n 50 -method PF -a 0.01 -epsilon ${eps} >> youtube_u-result.txt
    ./TopPPR -d youtube_u -algo SS -n 50 -method PW -a 0.01 -epsilon ${eps} >> youtube_u-result.txt
    ./TopPPR -d youtube_u -algo SS -n 50 -method PPF -a 0.01 -epsilon ${eps} -batch_size 5 >> youtube_u-result.txt
    ./TopPPR -d youtube_u -algo SS -n 50 -method PPW -a 0.01 -epsilon ${eps} -batch_size 5 >> youtube_u-result.txt
done

for eps in "${eps_arr[@]}"
do
    ./TopPPR -d youtube_u -algo SS -n 50 -method PF -a 0.01 -epsilon ${eps} -g undirected >> youtube_u-result.txt
    ./TopPPR -d youtube_u -algo SS -n 50 -method PPF -a 0.01 -epsilon ${eps} -batch_size 5 -g undirected >> youtube_u-result.txt
done

for eps in "${eps_arr[@]}"
do
    ./TopPPR -d youtube_u -algo PC -method MCW -a 0.01 -epsilon ${eps} >> youtube_u-result.txt
    ./TopPPR -d youtube_u -algo PC -method MCF -a 0.01 -epsilon ${eps} >> youtube_u-result.txt
    ./TopPPR -d youtube_u -algo PC -n 50 -method FORA -a 0.01 -epsilon ${eps} >> youtube_u-result.txt
    ./TopPPR -d youtube_u -algo PC -n 50 -method PF -a 0.01 -epsilon ${eps} >> youtube_u-result.txt
    ./TopPPR -d youtube_u -algo PC -n 50 -method PW -a 0.01 -epsilon ${eps} >> youtube_u-result.txt
    ./TopPPR -d youtube_u -algo PC -n 50 -method PPF -a 0.01 -epsilon ${eps} -batch_size 5 >> youtube_u-result.txt
    ./TopPPR -d youtube_u -algo PC -n 50 -method PPW -a 0.01 -epsilon ${eps} -batch_size 5 >> youtube_u-result.txt
done

for eps in "${eps_arr[@]}"
do
    ./TopPPR -d youtube_u -algo PC -n 50 -method PF -a 0.01 -epsilon ${eps} -g undirected >> youtube_u-result.txt
    ./TopPPR -d youtube_u -algo PC -n 50 -method PPF -a 0.01 -epsilon ${eps} -batch_size 5 -g undirected >> youtube_u-result.txt
done
