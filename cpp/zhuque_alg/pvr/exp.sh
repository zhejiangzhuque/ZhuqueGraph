#g++ -march=core2 -ffast-math -use_fast_math -pthread -std=c++11 -DSFMT_MEXP=607 -I SFMT-src-1.4.1/ -O3 -o ../build/TopPPR SFMT.c main.cpp 
rm log.log
rm youtube_u-result.txt
rm logtime.log

start=`date +%s`

eps_arr=(0.5 0.4 0.3 0.2 0.1)

for eps in "${eps_arr[@]}"
do
    ./TopPPR -d youtube_u -algo PC -method MCW -a 0.2 -epsilon ${eps} >> log.log
    ./TopPPR -d youtube_u -algo PC -method MCF -a 0.2 -epsilon ${eps} >> log.log
    ./TopPPR -d youtube_u -algo PC -n 50 -method FORA -a 0.2 -epsilon ${eps} >> log.log
    ./TopPPR -d youtube_u -algo PC -n 50 -method PF -a 0.2 -epsilon ${eps} >> log.log
    ./TopPPR -d youtube_u -algo PC -n 50 -method PW -a 0.2 -epsilon ${eps} >> log.log
    ./TopPPR -d youtube_u -algo PC -n 50 -method PPF -a 0.2 -epsilon ${eps} -batch_size 2 >> log.log
    ./TopPPR -d youtube_u -algo PC -n 50 -method PPW -a 0.2 -epsilon ${eps} -batch_size 2 >> log.log
done

end=`date +%s`
dif=$[ end - start ]
echo $dif

echo "**************************************************************************************py***********************************************************************"
rm youtube_u-result.txt

start=`date +%s`


python test.py >> log.log

end=`date +%s`
dif=$[ end - start ]
echo $dif