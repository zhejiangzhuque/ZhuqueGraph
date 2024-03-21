"""Self Play
"""
import sys
# sys.path.append("../rl/examples/battle_model/python")
# sys.path.append("../rl")


import argparse
import os

import tensorflow as tf
if type(tf.contrib) != type(tf): tf.contrib._warning = None
tf.compat.v1.logging.set_verbosity(tf.compat.v1.logging.ERROR)
import numpy as np
import magent
import pdb
from algorithm import spawn_ai
from algorithm import tools
from environment.battle_model.senario_battle import play
from tensorflow import set_random_seed

import time
date_running = time.strftime('%Y.%m.%d')
time_running = time.strftime('%H_%M_%S')


os.environ['CUDA_VISIBLE_DEVICES'] = '/gpu:4'
os.environ["TF_CPP_MIN_LOG_LEVEL"] = "3"

BASE_DIR = os.path.dirname(os.path.abspath(__file__))

def load_config(size):
    gw = magent.gridworld
    cfg = gw.Config()

    cfg.set({"map_width": size, "map_height": size})
    cfg.set({"minimap_mode": False})
    cfg.set({"embedding_size": 10})

    small = cfg.register_agent_type(
        "small",
        {'width': 1, 'length': 1, 'hp': 10, 'speed': 2,
         'view_range': gw.CircleRange(6), 'attack_range': gw.CircleRange(1.5),
         'damage': 2, 'step_recover': 0.1,

         'step_reward': -0.005,  'kill_reward': 5, 'dead_penalty': -0.1, 'attack_penalty': -0.1,
         })
    tiny = cfg.register_agent_type(
        "tiny",
        {'width': 1, 'length': 1, 'hp': 4, 'speed':1,
         'view_range': gw.CircleRange(6), 'attack_range': gw.CircleRange(1.5),
         'damage': 1, 'step_recover': 0.2,

         'step_reward': 0,  'kill_reward': 0, 'dead_penalty': -2, 'attack_penalty': -0.01,
         })

 
    g0 = cfg.add_group(tiny)
    g1 = cfg.add_group(small)

    a = gw.AgentSymbol(g0, index='any')
    b = gw.AgentSymbol(g1, index='any')

    # reward shaping to encourage attack
    cfg.add_reward_rule(gw.Event(a, 'attack', b), receiver=a, value=5)
    cfg.add_reward_rule(gw.Event(b, 'attack', a), receiver=b, value=0.2)

    return cfg

def train(env,args):
    np.random.seed(args.seed)
    set_random_seed(args.seed)
    
    parent_dir = f"{args.save_dir}/logs/{date_running}/{time_running}/"
    os.makedirs(parent_dir, exist_ok=True)
    # print(os.path.abspath(os.path.join(os.path.abspath(__file__),"../rl")))
    render_dir = os.path.join(parent_dir,"render")
    if not os.path.exists(render_dir):
        os.makedirs(render_dir)
    env.set_render_dir(render_dir)
    handles = env.get_handles()

    tf_config = tf.ConfigProto(allow_soft_placement=True, log_device_placement=False)
    tf_config.gpu_options.allow_growth = True

    log_dir = os.path.join(parent_dir,"tmp",'{}'.format(args.algo+'fixtiny3b6h'+str(args.seed)+'-'+str(args.len_nei)))
    model_dir = os.path.join(parent_dir,"models", '{}'.format(args.algo+'fixtiny256a6h'+str(args.seed)+'-'+str(args.len_nei)))
    
    
    # print('logdir',os.path.abspath(log_dir))
    sess = tf.Session(config=tf_config)
    main_model_dir = os.path.join(os.path.abspath(args.load_model_dir), '{}-0'.format(args.algo+'fixtiny64a6h'+str(args.seed)+'-'+str(args.len_nei)))
    oppo_model_dir = os.path.join(os.path.abspath(args.load_model_dir), '{}-1'.format('il'))
    main_msg_dir = os.path.join(os.path.abspath(args.load_model_dir), '{}-msg0'.format(args.algo+'fixtiny64a6h'+str(args.seed)+'-'+str(args.len_nei)))
    oppo_msg_dir = os.path.join(os.path.abspath(args.load_model_dir), '{}-msg1'.format(args.algo))

    # print(args.load_model_dir)

    models = [spawn_ai(args.algo, tf.train.AdamOptimizer, sess, env, handles[0], args.algo + '-me', args.max_steps,args.len_nei), spawn_ai('il',tf.train.AdamOptimizer,sess, env, handles[1],'il-oppo', args.max_steps,args.len_nei)]

    if args.usemsg!='None':
        MsgModels = [spawn_ai('msgdqn', tf.train.AdamOptimizer, sess, env, handles[0], 'msgdqn'+ '-me', args.max_steps), spawn_ai('msgdqn',tf.train.AdamOptimizer, sess, env, handles[1], 'msgdqn' + '-opponent', args.max_steps)]
    else:
        print('do not use msg models')
        MsgModels=[None,None]
    models.extend(MsgModels)
    sess.run(tf.global_variables_initializer())

    if args.idx!='None':
        print('load successfully'+str(args.idx))
        models[1].load(oppo_model_dir, step=args.idx)

    if args.algo in ['mfq', 'mfac']:
        use_mf = True
    else:
        use_mf = False 

    start_from = 0

    runner = tools.Runner(sess, env, handles, args.map_size, args.max_steps, models[:-2],models[-2:] ,play,
                            render_every=args.save_every if args.render else 0, save_every=args.save_every, tau=0.01, log_name=args.algo,
                            log_dir=log_dir, model_dir=model_dir,render_dir=render_dir, train=True ,len_nei=args.len_nei,rewardtype= 'self',crp=True,is_selfplay=False,is_fix=True)

    for k in range(start_from, start_from + args.n_round):
        eps = magent.utility.piecewise_decay(k, [0, 700,1400, 2000], [1, 0.3,0.01,0.01])
        # eps = 0
        runner.run(eps, k)
