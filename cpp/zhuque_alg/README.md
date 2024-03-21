#### 项目简介
cpp/zhuque_alg/下是10个图算法源码

* BiHSE biplex --- Hereditary Cohesive Subgraphs Enumeration on Bipartite Graphs: The Power of Pivot-based Approaches. SIGMOD, 2023.
* Biclique --- Efficient Biclique Counting in Large Bipartite Graphs. SIGMOD, 2023.
* DefectiveClique --- Maximal Defective Clique Enumeration. SIGMOD, 2023.
* dpcolor --- Lightning Fast and Space Efficient k-clique Counting. WWW, 2022.
* fairnessclique --- Fairness-aware Maximal Clique Enumeration. ICDE, 2022
* maximalUncertainClique --- Maximal Clique Enumeration on Uncertain Graphs: A Pivot-based Approach. SIGMOD, 2022.
* MBC --- Mining Bursting Core in Large Temporal Graph. Proc. VLDB Endow. 15(13): 3911-3923, 2022.
* pvr --- Efficient Personalized PageRank Computation: The Power of Variance-Reduced Monte Carlo Approaches. SIGMOD, 2023.
* Resistance-Landmark --- Efficient Resistance Distance Computation: The Power of Landmark-based Approaches. SIGMOD, 2023.
* RSFPPR --- Efficient Personalized PageRank Computation: A Spanning Forests Sampling Based Approach. SIGMOD, 2022.

#### 环境准备
##### pybind
```
pip install "pybind11[global]"
```
##### Eigen
```
sudo apt-get install libeigen3-dev
http://eigen.tuxfamily.org/index.php?title=Main_Page
```
##### 要确定你的Eigen3安装在/usr/local/include还是/usr/include
```
cd /usr/include
sudo ln -sf eigen3/Eigen Eigen
sudo ln -sf eigen3/unsupported unsupported
```
#### 编译整个项目
```
 mkdir build/
 cd build/
 rm -rf * && cmake .. && make
 ```
##### 编译结果
编译后build/下生成可以直接被Python调用的so库，调用过程参见目录zhuque_graph/analytics 
