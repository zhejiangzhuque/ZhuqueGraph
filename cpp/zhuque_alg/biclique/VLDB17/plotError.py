import sys

f1 = sys.argv[1]
f2 = sys.argv[2]

print(f1)
print(f2)

ps = []
qs = []
vs = []

with open(f1, 'r') as f:
    lines = f.readlines()

    for line in lines:
        if line.count('-') > 0 and line.count(':') > 0:
            p1 = line.index('-')
            p2 = line.index(':')

            p = int(line[:p1])
            q = int(line[p1 + 1 : p2])

            if p == 1 or q == 1 or p >= 100 or q >= 100:
                continue

            v = int(line[p2 + 1:-1])

            # if p==2 and q==2:
            #     print line

            ps.append(p)
            qs.append(q)
            vs.append(v)

ps2 = []
qs2 = []
vs2 = []

tm = 0
sumE = 0.0

with open(f2, 'r') as f:
    lines = f.readlines()

    for line in lines:
        if line.startswith("totaltime:"):
            # print("there")
            i = line.index(':')
            j = line.index('ms')
            # print(i, j)

            # print(line[i + 1:j])

            tm = int(line[i + 1:j])

        if line.count('-') > 0 and line.count(':') > 0:
            p1 = line.index('-')
            p2 = line.index(':')

            p = int(line[:p1])
            q = int(line[p1 + 1 : p2])

            if p == 1 or q == 1:
                continue
                
            if line[p2 + 1] == ' ':
                p2 += 1
            
            v = int(float(line[p2 + 1:-1]))
            if line[-1].isdigit():
                v = int(float(line[p2 + 1:-1]))

            # ps2.append(p)
            # qs2.append(q)
            # vs2.append(v)
            for i in range(len(ps)):
                if ps[i] == p and qs[i] == q:
                    if p != 2:
                        sumE += 1.0 * abs(v - vs[i]) / vs[i]
                    print '{}-{}'.format(p,q), "{:.2f}%".format(100*(1.0 * abs(v - vs[i]) / vs[i]))

print "tm:", tm
print "avgErr:", sumE/ 6
