#!/bin/sh

make

./main --filename facebook --file_name ./graphs/facebook.txt --method source_landmark --pairs node > facebook-result-sourcelandmark.txt
./main --filename facebook --file_name ./graphs/facebook.txt --method source --pairs node --num_trials 5 > facebook-result-source.txt

./main --filename facebook --file_name ./graphs/facebook.txt --method akp --pairs vertex-pairs > facebook-result-akp.txt
./main --filename facebook --file_name ./graphs/facebook.txt --method commute --pairs vertex-pairs > facebook-result-commute.txt
./main --filename facebook --file_name ./graphs/facebook.txt --method abwalk --pairs vertex-pairs > facebook-result-abwalk.txt
./main --filename facebook --file_name ./graphs/facebook.txt --method localtree --pairs vertex-pairs > facebook-result-localtree.txt
./main --filename facebook --file_name ./graphs/facebook.txt --method push --pairs vertex-pairs > facebook-result-push.txt
./main --filename facebook --file_name ./graphs/facebook.txt --method bipush --pairs vertex-pairs > facebook-result-bipush.txt