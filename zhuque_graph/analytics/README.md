#### 项目简介
一一对应cpp/zhuque_alg/下的图算法
#### 运行
```
python3 xxxxx.py --help  #查看参数信息
python3 xxxxx.py #默认参数运行
```
#### 问题排查
调用时可能存在import找不到so库的问题ModuleNotFoundError

1. 首先确保so的路径在sys.path内，可以通过sys.path.append()添加
2. <so.name>.cpython-<python-version>-<cpu-platform>-linux-gnu.so 检查版本和so库名称是否正确
3. 确认当前Python支持的版本
```
import importlib.machinery
print(importlib.machinery.all_suffixes())
```
