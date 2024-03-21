#!/bin/sh
landmark_arr=(degree core pagerank random)

make

for landmark in "${landmark_arr[@]}"
do
    ./main --filename facebook --file_name ./graphs/facebook.txt --method abwalk --pairs vertex-pairs --landmark ${landmark} >> facebook-landmark.txt
    ./main --filename facebook --file_name ./graphs/facebook.txt --method localtree --pairs vertex-pairs --landmark ${landmark} >> facebook-landmark.txt
    ./main --filename facebook --file_name ./graphs/facebook.txt --method bipush --pairs vertex-pairs --landmark ${landmark} >> facebook-landmark.txt
    ./main --filename facebook --file_name ./graphs/facebook.txt --method source_landmark --pairs node --landmark ${landmark} >> facebook-landmark.txt
    ./main --filename facebook --file_name ./graphs/facebook.txt --method source --pairs node --num_trials 50 --landmark ${landmark} >> facebook-landmark.txt
    ./main --filename facebook --file_name ./graphs/facebook.txt --method push --pairs vertex-pairs --landmark ${landmark} >> facebook-landmark.txt
done