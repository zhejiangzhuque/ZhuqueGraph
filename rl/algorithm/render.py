import warnings
import tensorflow as tf
import os
warnings.filterwarnings('ignore')
warnings.filterwarnings("ignore", category=DeprecationWarning, module='tensorflow')
if type(tf.contrib) != type(tf): tf.contrib._warning = None
tf.compat.v1.logging.set_verbosity(tf.compat.v1.logging.ERROR)

os.environ["TF_CPP_MIN_LOG_LEVEL"] = "3"
os.environ['CUDA_VISIBLE_DEVICES'] = '/gpu:4'

from PIL import Image, ImageDraw,ImageFont
import json
import os
import argparse
import time
import sys
sys.path.append(os.path.abspath(os.path.join(os.path.abspath(__file__),"../../environment/battle_model/")))
sys.path.append(os.path.abspath(os.path.join(os.path.abspath(__file__),"../../")))
sys.path.append(os.path.abspath(os.path.join(os.path.abspath(__file__),"../../../")))

# print(os.path.abspath(os.path.join(os.path.abspath(__file__),"../environment/battle_model/")))

import argparse
import os
import magent
from algorithm import spawn_ai
from algorithm import tools
from environment.battle_model.senario_battle import play

date_running = time.strftime('%Y.%m.%d')
time_running = time.strftime('%H_%M_%S')

def load_config(size, args):
    gw = magent.gridworld
    cfg = gw.Config()

    cfg.set({"map_width": size, "map_height": size})
    cfg.set({"minimap_mode": False})
    cfg.set({"embedding_size": 10})

    # small is blue
    small = cfg.register_agent_type(
        "small",
        {'width': 1, 'length': 1, 'hp': args.blue_hp, 'speed': 2,
         'view_range': gw.CircleRange(6), 'attack_range': gw.CircleRange(1.5),
         'damage': args.blue_damage, 'step_recover': args.blue_step_recover,

         'step_reward': -0.005,  'kill_reward': 5, 'dead_penalty': -0.1, 'attack_penalty': -0.1,
         })
    
    # tiny is red
    tiny = cfg.register_agent_type(
        "tiny",
        {'width': 1, 'length': 1, 'hp': args.red_hp, 'speed':1,
         'view_range': gw.CircleRange(6), 'attack_range': gw.CircleRange(1.5),
         'damage': args.red_damage, 'step_recover': args.red_step_recover,

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

def render(args):
    global save_dir

    env=magent.GridWorld(load_config(size=40, args=args))
    
    print("Initialize env successfully!")
    if not os.path.exists(save_dir):
        os.makedirs(save_dir)
    env.set_render_dir(save_dir)
    handles = env.get_handles()

    tf_config = tf.ConfigProto(allow_soft_placement=True, log_device_placement=False)
    tf_config.gpu_options.allow_growth = True

    sess = tf.Session(config=tf_config)
    main_model_dir = os.path.join(os.path.abspath(args.load_model_dir),f"{args.load_model_iteration}"+'-0')
    oppo_model_dir = os.path.join(os.path.abspath("../data/models"), '{}-1'.format('il'))
    main_msg_dir = os.path.join(os.path.abspath(args.load_model_dir),f"{args.load_model_iteration}"+'-msg0')

    models = [spawn_ai(args.algo, tf.train.AdamOptimizer, sess, env, handles[0], args.algo + '-me', args.max_steps,args.len_nei), spawn_ai('il',tf.train.AdamOptimizer,sess, env, handles[1],'il-oppo', args.max_steps,args.len_nei)]

    MsgModels = [spawn_ai('msgdqn', tf.train.AdamOptimizer, sess, env, handles[0], 'msgdqn'+ '-me', args.max_steps), spawn_ai('msgdqn',tf.train.AdamOptimizer, sess, env, handles[1], 'msgdqn' + '-opponent', args.max_steps)]
    models.extend(MsgModels)
    sess.run(tf.global_variables_initializer())
    print("Initialize model successfully!")

    if args.idx!='None':
        print('load successfully'+str(args.idx))
        models[1].load(oppo_model_dir, step=args.idx)

    if args.usemsg!='None':
       models[0].load(main_model_dir,step='6-0selfcrpw')
       MsgModels[0].load(main_msg_dir,step='6-0selfcrpw')
    else:
       models[0].load(main_model_dir,step='6-0selfnomnw')

    start_from = 0

    print("Running")
    runner = tools.Runner(sess, env, handles, args.map_size, args.max_steps, models[:-2],models[-2:] ,play,
                            render_every=1, save_every=1, tau=0.01, log_name=args.algo,
                            log_dir=None, model_dir=None,render_dir=save_dir, train=False,len_nei=args.len_nei,rewardtype= 'self',crp=True,is_selfplay=False,is_fix=True,agents_num=args.agents_num)

    for k in range(start_from, start_from + args.render_num):
        eps = 0
        # print(k)
        runner.run(eps, k)

class Wall:
    def __init__(self, images, width, height, coordinates) -> None:
        self.images = images
        self.width = width
        self.height = height
        self.coordinates = coordinates

    def draw(self):
        for coordinate in self.coordinates:
            x, y = coordinate
            for image in self.images:
                image.putpixel((x, y), (127, 127, 127,1))  # 灰色

class View:
    def __init__(self,args,config_path,info_path,save_path) -> None:
        # 图像大小
        self.args = args
        self.background_color = (255, 255, 255)
        self.scale_coef = 10
        self.save_path = save_path
        self.config_path = config_path
        self.info_path = info_path
        with open(os.path.join(config_path),"r") as f:
            config = json.load(f)
            self.width = config["width"]
            self.height = config["height"]
        self.coordinates = []
        self.wall_coordinates = []
        self.wall_flag = False
        self.agent_flag = False
        self.num_coordinates = []
        self.agent_counts = []

    def get_coordinates(self,):
        # 解析信息行和图形坐标
        for info_file in os.listdir(self.info_path):
            if info_file == "config.json":
                continue
            # 解析信息行和图形坐标
            with open(os.path.join(self.info_path,info_file), "r") as file:
                lines = file.readlines()
                for line in lines:
                    elements = line.strip().split(" ")
                    if elements[0] == "F" or elements[0] == "W":
                        if elements[0] == "F":
                            num_coordinate = int(elements[1])
                            self.num_coordinates.append(num_coordinate)
                            self.agent_counts.append([0,0])
                            self.agent_flag = True
                            self.wall_flag = False
                        elif elements[0] == "W":
                            self.agent_flag = False
                            self.wall_flag = True
                    else:
                        if elements[0] == "0" and len(elements) == 4:
                            continue
                        if self.agent_flag:
                            color = (255, 0, 0) if int(elements[0]) < self.args.agents_num[0] else (0, 0, 255)
                            x = int(elements[3])
                            y = int(elements[4])
                            self.coordinates.append((x, y, color))
                            if int(elements[0]) < self.args.agents_num[0]:
                                self.agent_counts[-1][0] += 1
                            else:
                                self.agent_counts[-1][1] += 1
                        if self.wall_flag:
                            x = int(elements[0])
                            y = int(elements[1])
                            self.wall_coordinates.append((x, y))

    def draw(self):
        # 创建图像对象
        self.get_coordinates()
        images = [Image.new("RGB", (self.width, self.height), self.background_color) for _ in range(len(self.num_coordinates))]
        wall = Wall(images, self.width, self.height, self.wall_coordinates)
        wall.draw()

        # 绘制坐标点
        idx = 0
        for image_idx, num_coordinate in enumerate(self.num_coordinates):
            for coordinate_idx in range(idx, num_coordinate + idx):
                x, y, color = self.coordinates[coordinate_idx]
                images[image_idx].putpixel((x, y), color)

            idx += num_coordinate

        # 在图像上显示智能体数量
        for i, image in enumerate(images):
            draw = ImageDraw.Draw(image)

        # 图形放大
        resized_images = [image.resize((self.width*self.scale_coef, self.height*self.scale_coef), resample=Image.NEAREST) for image in images]
        # 保存为GIF动画
        for i, image in enumerate(resized_images):
            draw = ImageDraw.Draw(image)
            count_text = f"Red Num:  {self.agent_counts[i][0]}\nBlue Num: {self.agent_counts[i][1]}"
            draw.text((20, 20), count_text, fill=(0, 0, 0))
        resized_images[0].save(self.save_path, save_all=True, append_images=resized_images[1:], optimize=False, duration=200, loop=0)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("--load_model_dir",type=str,default="../data/2023.07.04/22_29_01/models/gilfixtiny256a6h3000-6")
    parser.add_argument("--load_model_iteration",type=int,default=1900)
    parser.add_argument("--save_dir",type=str,default="../data")
    parser.add_argument("--render_num",type=int,default=5)
    
    parser.add_argument("--red_num",type=int,default=40)
    parser.add_argument("--blue_num",type=int,default=40)
    parser.add_argument("--red_damage",type=int,default=1)
    parser.add_argument("--blue_damage",type=int,default=2)
    parser.add_argument("--red_hp",type=int,default=4)
    parser.add_argument("--blue_hp",type=int,default=10)
    parser.add_argument("--red_step_recover",type=int,default=0.2)
    parser.add_argument("--blue_step_recover",type=int,default=0.1)

    parser.add_argument('--usemsg',type=str,default='True')
    parser.add_argument('--idx', type=str, default='6-1995selfnomnw')
    parser.add_argument('--len_nei',type=int,default=6)
    parser.add_argument('--max_steps', type=int, default=200, help='set the max steps')
    parser.add_argument('--map_size', type=int, default=40, help='set the size of map')

    args = parser.parse_args()
    save_dir = f"{args.save_dir}/logs/{date_running}/{time_running}/"
    args.agents_num = [args.red_num,args.blue_num]
    args.algo = "gil"

    render(args)

    config_path = os.path.join(save_dir,"config.json")
    save_path = os.path.join(save_dir,"image.gif")
    view = View(args,config_path,save_dir,save_path)
    view.draw()