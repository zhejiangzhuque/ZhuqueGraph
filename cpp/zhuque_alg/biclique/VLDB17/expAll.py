import os
import sys
import commands 

s = "git  stackoverflow	mu	imdb	actor2	ama	dblp"
ss = s.split()

ans = ['' for i in range(len(ss))]
def walkFile2(file):
    global ans
    for root, dirs, files in os.walk(file):
        i = 0

        for f in files:
            if f.endswith(".txt") == False:
                continue
            
            j = 0
            while(j < len(ss) and f.count(ss[j]) == 0):
                j += 1
            if j == len(ss):
                continue

            totalF = os.path.join(root, f)
            ans[j] = totalF

walkFile2("../data/")


for i in range(len(ss)):
    run = "timeout 24h ./run"
    run += ' -f ' + ans[i]
    run += ' -e 500000'
    run += ' -all'
    run += ' > logs/'+ss[i]+'_5e5.txt'

    # print run 
    print commands.getstatusoutput(run)


