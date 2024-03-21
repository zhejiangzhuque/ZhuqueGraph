### 文件依赖环境

* python=3.6.3

* conda 

  * ```txt
    _libgcc_mutex=0.1
    _openmp_mutex=4.5
    ca-certificates=2023.5.7
    cudatoolkit=10.0.130
    cudnn=7.6.5.32
    libgcc-ng=13.1.0
    libgomp=13.1.0
    libstdcxx-ng=13.1.0
    ncurses=5.9
    openssl=1.0.2u
    pip=21.3.1
    python=3.6.3
    python_abi=3.6
    readline=7.0
    setuptools=58.0.4
    sqlite=3.20.1
    tk=8.6.11
    wheel=0.37.1
    xz=5.2.6
    zlib=1.2.11
    ```

* pip 

  * ```txt
    absl-py==0.7.0
    astor==0.8.1
    cached-property==1.5.2
    certifi==2023.5.7
    charset-normalizer==2.0.12
    dataclasses==0.8
    decorator==4.4.2
    dm-sonnet==1.23
    gast==0.5.4
    grpcio==1.48.2
    h5py==3.1.0
    idna==3.4
    imageio==2.15.0
    importlib-metadata==4.8.3
    joblib==1.1.1
    Keras-Applications==1.0.8
    Keras-Preprocessing==1.1.2
    Markdown==3.3.7
    mock==5.0.2
    networkx==2.5.1
    numpy==1.19.5
    onnx==1.14.0
    pandas==1.1.5
    Pillow==8.4.0
    protobuf==4.21.0
    python-dateutil==2.8.2
    pytz==2023.3
    requests==2.27.1
    scikit-learn==0.24.2
    scipy==1.5.4
    six==1.16.0
    sklearn==0.0
    tensorboard==1.13.1
    tensorflow-estimator==1.13.0
    tensorflow-gpu==1.13.1
    tensorflow-probability==0.5.0
    termcolor==1.1.0
    tf2onnx==1.7.2
    threadpoolctl==3.1.0
    typing_extensions==4.1.1
    urllib3==1.26.16
    Werkzeug==2.0.3
    zipp==3.6.0
    ```

### 使用说明

1. 创建虚拟环境

```
conda create -n rl_env python=3.6.3
source activate rl_env
```

2. 安装cuda

```txt
conda install cudatoolkit=10.0.130
conda install cudnn=7.6.5.32
conda install libprotobuf
conda install -c anaconda cmake
```

3. 安装第三方库

```txt
pip install -r requirements.txt
```

4. run

```cmd
cd LSC/launch/LSCTrain
python LSCTrain.py
```

训练参数如下所示：

* --save_every 保存模型和可视化信息的频率
* --update_every 模型更新频率
* --n_round 训练轮数
* --render 是否保存可视化渲染信息
* --map_size 场景图大小
* --max_steps 每轮博弈的最大步长
* --len_nei 通信距离
* --usemsg 是否使用通信机制
* --idx 对手模型参数
* --seed 随机数种子
* --load_model_dir 导入敌方模型文件目录
* --save_dir 保存模型及可视化信息目录

训练完后，在LSC/rl/data目录下存放了运行代码时刻所保存的日志文件和模型文件。

点进某个时间节点目录，里面存放有models文件夹、render文件夹和tmp文件夹，如下所示。

![pCtNc8I.png](https://s1.ax1x.com/2023/06/24/pCtNc8I.png)

models文件夹里存放保存的模型文件，有onnx格式；render文件夹中存放可以可视化的文件；tmp存放了tensorboard的训练文件（在gil文件夹中）和训练日志logging.txt，日志内容如下所示。

```txt
Run: 1, {'ave_agent_reward': 0.047252174234469944, 'total_reward': 1134.1599444560707, 'kill': 0, 'oppo-ave-agent-reward': -0.04172265775054257, 'oppo-total-reward': -1068.1000384138897, 'oppo-kill': 6}
```

* Run 后面为训练轮次
* ave_agent_reward 为每一个智能体的平均奖励
* total_reward 为所有智能体的奖励和
* kill 指智能体击杀数量
* oppo-ave-agent-reward 为对方智能体的平均奖励
* oppo-total-reward 为对方智能体的奖励和
* oppo-kill 指对方智能体的击杀数量

查看可视化奖励曲线使用tensorboard，如下所示。

```cmd
tensorboard --logdir= 以events开头文件的路径（如/root/hzx/LSC/rl/data/2023.06.13/16_04_40/tmp/gilfixtiny3b6h3000-6/gil/events.out.tfevents.1686643487.node1）
```

5. 可视化

```cmd
cd /LSC/rl/algorithm
python render.py --load_model_iteration 1900
```

关键参数：

* --load_model_dir 需要可视化的模型文件目录，默认为"../data/2023.07.04/22_29_01/models/gilfixtiny256a6h3000-6"，如需更改，请根据默认路径的格式进行修改
* --load_model_iteration 需要可视化的相应训练迭代次数的模型文件目录
* --save_dir 保存可视化文件的根目录
* --render_num 需要渲染的局数
* --red_num 红方智能体的数量
* --blue_num 蓝方(敌方)智能体的数量
* --red_damage 红方智能体的攻击力
* --blue_damage 蓝方智能体的攻击力
* --red_hp 红方智能体的血量
* --blue_hp 蓝方智能体的血量
* --red_step_recover 红方智能体每步回复的血量
* --blue_step_recover 蓝方智能体每步回复的血量
* --max_steps 推演可视化时的最大步长，若当局博弈超过这个步长还未结束则会强制终止
* 其余参数无需改变，按照默认参数即可。

运行命令之后，会在save_dir的文件目录下生成一个日期/时间的多级目录，里面内容如下所示。其中config.json为生成的config文件，image.gif为生成的动图，其余内容为每一局游戏的可视化文本信息。在可视化动图中，红方为训练方，蓝方为敌方。

![pCgYxjf.png](https://s1.ax1x.com/2023/07/08/pCgYxjf.png)

### 文档目录结构

LSC 

* launch 存放执行代码
  * LSCTrain
    * LSC_utils.py
    * LSCTrain.py 为执行代码
* rl
  * data 存放生成的模型、可视化及日志文件，如上方所述。
  * algorithm 存放强化学习算法文件
    * base.py 
    * q_leanrning.py 存放MsgDQN、DQN、GDQN三种强化学习算法
    * tools.py
  * environment 存放强化学习环境代码文件
    * battle_model
      * build 用于build多智能体强化学习环境
      * magent 存放magent多智能体强化学习环境
      * src 强化学习环境必要的库
* zhuque_graph
  * dataloader 空
  * nn
    * tensorflow
      * core 空
      * layer 存放构建图网络等必要模块
      * model 存放多个图模型
  * sample 空
  * utils 存放编写的函数工具集，方便调用

### 数据效果

* 数据集

  * MAgent强化学习环境
  * MPE强化学习环境

* 指标及对比测试结果

  * 指标说明
  
    * $N_s$ : number of successive reaching of three agents to one landmark
    * $N_o$ : number of successive reaching of more than three agents to one landmark
    * Mean-reward : average per-episode reward of all agents
    *  $r_{kd}$ : kill to death ratio $N_k/N_d$

  * MAgent
  
  * | 方法        | LSC      | LSC-star | LSC-nbor | IDQN |
    | ----------- | -------- | -------- | -------- | ---- |
    | Mean-reward | **1.13** | 1.04     | 0.94     | 0.89 |
    | $N_k$       | **62.9** | 61.36    | 62.4     | 59.3 |
    | $N_d$       | **39.4** | 48.52    | 43.6     | 55.8 |
    | $r_{kd}$    | **1.60** | 1.21     | 1.43     | 1.06 |
  
  * MPE
  
  * | 方法        | LSC       | LSC-star | LSC-nbor | IDQN  |
    | ----------- | --------- | -------- | -------- | ----- |
    | $N_s$       | **156**   | 80       | 61       | 16    |
    | $N_o$       | 13        | 16       | 9        | **2** |
    | Mean-reward | **-54.6** | -68.2    | -69.1    | -78.1 |
  
    
