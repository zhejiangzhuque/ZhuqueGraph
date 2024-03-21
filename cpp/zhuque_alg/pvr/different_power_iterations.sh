g++ -march=core2 -ffast-math -use_fast_math -pthread -std=c++11 -DSFMT_MEXP=607 -I SFMT-src-1.4.1/ -O3 -o TopPPR SFMT.c main.cpp

l1_err_arr=(1e-3 1e-4 1e-5 1e-6 1e-7)

rm DPI.txt

for l1_error in "${l1_err_arr[@]}"
do
    ./TopPPR -d orkut -algo SS -n 1 -method DPI_naive -a 0.2 -l1_error ${l1_error} -g undirected >> DPI.txt
    ./TopPPR -d orkut -algo SS -n 1 -method DPI -a 0.2 -l1_error ${l1_error} -g undirected >> DPI.txt
    ./TopPPR -d orkut -algo SS -n 1 -method DPI_naive -a 0.01 -l1_error ${l1_error} -g undirected >> DPI.txt
    ./TopPPR -d orkut -algo SS -n 1 -method DPI -a 0.01 -l1_error ${l1_error} -g undirected >> DPI.txt
done