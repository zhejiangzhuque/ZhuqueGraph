"""Self Play
"""
import sys
import os

sys.path.append(os.path.abspath(os.path.join(os.path.abspath(__file__),"../../../rl/environment/battle_model/")))
sys.path.append(os.path.abspath(os.path.join(os.path.abspath(__file__),"../../../rl/environment")))
sys.path.append(os.path.abspath(os.path.join(os.path.abspath(__file__),"../../../rl/algorithm")))
sys.path.append(os.path.abspath(os.path.join(os.path.abspath(__file__),"../../../rl/")))


import argparse
import os
import magent
import warnings
warnings.filterwarnings('ignore')
warnings.filterwarnings("ignore", category=DeprecationWarning, module='tensorflow')


from LSC_utils import train,load_config

os.environ["TF_CPP_MIN_LOG_LEVEL"] = "3"

# os.environ['CUDA_VISIBLE_DEVICES'] = '/gpu:4'

def main(args):
    env=magent.GridWorld(load_config(size=args.map_size))
    # print(dir(env))
    train(env,args)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--save_every', type=int, default=100, help='decide the self-play update interval')
    parser.add_argument('--update_every', type=int, default=5, help='decide the udpate interval for q-learning, optional')
    parser.add_argument('--n_round', type=int, default=2000, help='set the trainning round')
    parser.add_argument('--render', default='True', action='store_true', help='render or not (if true, will render every save)')
    parser.add_argument('--map_size', type=int, default=40, help='set the size of map')  # then the amount of agents is 64
    parser.add_argument('--max_steps', type=int, default=400, help='set the max steps')
    parser.add_argument('--len_nei',type=int,default=6)
    parser.add_argument('--usemsg',type=str,default='True')
    parser.add_argument('--idx', type=str, default='6-1995selfnomnw')
    parser.add_argument('--seed',type=int,default=3000)
    parser.add_argument("--load_model_dir",type=str,default="../../rl/data/models")
    parser.add_argument("--save_dir",type=str,default="../../rl/data")

    args = parser.parse_args()
    args.algo = "gil"

    main(args)
