# Code for MBC Algorithm

This repository contains a reference implementation of the algorithm for Mining Bursting Cores in Temporal Networks.

## Environment setup

Codes run on Python 2.7 or later. [PyPy](http://pypy.org/) compiler is recommended because it can make the computations quicker without change the codes.

You may use Git to clone the repository from
GitHub and run it manually like this:

    git clone https://github.com/qinhc/MBC.git
    cd MBC
    pip install click
    python run.py 

## Dateset description
We focus on mining the temporal networks so each edges are associated with a timestamp. Temporal edges are stored at the raw data in which each line is one temporal edge.
 
| from_id | \t  | to_id    | \t  |  timestamps  |
| :----:  |:----: | :----:   |:----:   | :----: |

## Running example
You can type in dataset number, parameters l, delta and method number to control the program:

    Type one number to chose the dataset: [1]Chess; [2]Lkml; [3]Enron. (int): 1
    l(int): 3
    Delta(int): 3
    Type one number to chose the algorithm: [1]MBC; [2]MBCPLUS; [3]POMBC. (int): 2
    max degree:263
    n:7301
    m:55899.0
    temporal edges:63689.0
    T:101
    Graph loaded!
    kCore Time:0:00:00.035273
    length(Vc):5472
    {9, 39, 40, 41, 81, 101, 114, 127, 128, 130, 152, 156, 223, 269, 668, 669, 712, 713, 2195, 2209, 2412, 2797, 3036, 3037, 3038, 3043, 3116, 3139, 3474, 3574, 3837, 3899, 3962, 3973, 3987, 4046, 4801, 4954, 4966}
    0:00:02.548877
    
## Tips

Due to the limit of space, we only upload some datasets of small size here. Welcome to e-mail me for more datasets of temporal networks.