# Contents
This repository contains implementations of algorithms for maximal s-defective clique enumeration.

# Compile

```
make
```

# Usage
To execute the code, you need to run the following executable file:

```
./sdclique branch|delay [filepath] -a=[1,7] -k=[1,5] -q=[2,n]
```
For example:
```
./sdclique branch datas/ca-GrQc.txt -k=1 -a=6 -q=10
```

This executable accepts the following optional parameters when using 'branch':
- "-a=": This is the executed algorithm, where "-a=1|2|3" for basic branching algorithms, "-a=4|5" for pivot-based branching algorithms enuemrating all maximal s-defective cliques, "-a=6|7" for pivot-based branching algorithms enuemrating relatively-large maximal s-defective cliques.

- "-k=": This is to set parameter 's'.

- "-q=": This is the constraint on the min-size of maximal s-defective cliques to be enumerated.
