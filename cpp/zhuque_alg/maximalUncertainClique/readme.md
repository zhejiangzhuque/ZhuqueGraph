# Contents
This repository contains implementations of algorithms for maximal clique enumeration on uncertain graphs, appeared in:

*SIGMOD2022* - Fast Maximal Clique Enumeration on Uncertain Graphs: A Pivot-based Approach 

# Compile

```
make
```

# Usage
To execute the code, you need to run the following executable file:

```
./uclique [filepath] -a=[1,2] -k=[0,n] -e=[0,1]
```

This executable accepts the following optional parameters:
- "-a=": This is the executed algorithm, "-a=1" for pivot algorithm with topCore reduction and "-a=2" for pivot algorithm with topTriangle reduction.

- "-k=": This is the minsize constraint.

- "-e=": This is the clique probability constraint.
