import random
import math
import numpy as np
import time
import tensorflow as tf
import pdb 
import json
import pickle
import pandas as pd



def generate_map(env, map_size, handles, init_num=[0,0]):
    """ generate a map, which consists of two squares of agents"""
    width = height = map_size
    init_num = [map_size * map_size * 0.04,map_size * map_size * 0.04] if init_num == [0,0] else init_num
    gap = 3
    leftID=0
    rightID=1

    # print(init_num)

    # left
    n = init_num[0]
    side = int(math.sqrt(n)) * 2
    pos = []
    # for x in range(width//2 - gap - side, width//2 - gap - side + side, 2):
    #     for y in range((height - side)//2, (height - side)//2 + side, 2):
    #         pos.append([x, y, 0])
    i = 0
    while i < n:
        x = random.randint(width//2 - gap - side, width//2 - gap - side + side)
        y = random.randint((height - side)//2, (height - side)//2 + side)
        if [x,y,0] in pos:
            continue
        pos.append([x, y, 0])
        i += 1
    env.add_agents(handles[leftID], method="custom", pos=pos)
    
    # unique_data = list(set(map(tuple, pos)))
    # unique_length = len(unique_data)
    # print(unique_length)

    # right
    n = init_num[1]
    side = int(math.sqrt(n)) * 2
    gap=3
    pos = []
    # for x in range(width//2 + gap, width//2 + gap + side, 2):
    #     for y in range((height - side)//2, (height - side)//2 + side, 2):
    #         pos.append([x, y, 0])

    i = 0
    while i < n:
        x = random.randint(width//2 + gap, width//2 + gap + side)
        y = random.randint((height - side)//2, (height - side)//2 + side)
        if [x,y,0] in pos:
            continue
        pos.append([x, y, 0])
        i += 1
    # print(len(pos))
    env.add_agents(handles[rightID], method="custom", pos=pos)

    # nums = [env.get_num(handle) for handle in handles]
    # print(nums)
    # print('done')
    # print(len(pos))

def get_distance_matrix(env,handle,nei_len):
    poss=env.get_pos(handle)
    
    nei_space=nei_len**2
    poss_all=np.array([poss]*len(poss))
    poss_mine=np.array([ [t]*len(poss) for t in poss])
    indexs=((poss_all-poss_mine)**2).sum(axis=2)
    indexs=indexs-np.ones_like(indexs)*nei_space
    indexs[indexs<0]=0
    return indexs,poss

def check_isbetter(W,indexs,chs):
    chs=np.array(chs)
    for ch in chs:
        # pdb.set_trace()
        a=np.where(indexs[ch]==0)[0]
        tmp=[val for val in a if val in chs and val!=ch]
        if len(tmp)>0:
            for t in tmp:
                if W[ch]<=W[t]:
                    chs=np.delete(chs,np.where(chs==ch)[0][0])
                    break
    return chs

def form_chs(W,indexs,un):
    ch=[]
    uns=np.array(un)
    while len(uns)!=0:
        # 
        a=uns[np.argmax(W[uns])]
        ch.append(a)
        uns=np.array(list(set(uns)-set(np.where(indexs[a]==0)[0])))
        # pdb.set_trace()
        # p.delete(uns,)
    return list(set(ch))


def ids_to_ind(id,ids):
    ind=[]
    for i in id:
        a=np.where(ids==i)[0]
        if list(a) == []:
            continue
        else: 
            ind.append(a[0])
    return ind

def ind_to_id(ind,ids):
    id=[]
    for i in ind:
        id.append(ids[i])
    return id

def get_chs(indexs,hps,chs_id,ids,weight=None):
    chs=ids_to_ind(chs_id,ids)
    # pdb.set_trace()
    # s=time.time()
    if list(weight):
        W=weight
    else:
        W=calc_w(indexs,hps)
    # pdb.set_trace()
    # downgrade nearby chs based on their W
    chs=check_isbetter(W,indexs,chs)
    # find undecided nodes
    un=set([i for i in range(len(indexs))])
    for ch in chs:
        un=un.intersection(set(np.where(indexs[ch]>0)[0]))
    un=list(un)
    
    un=form_chs(W,indexs,un)
    # print('getchs',time.time()-s)
    chs=list(set(chs).union(set(un)))
    # pdb.set_trace()
    chs_id = ind_to_id(chs,ids)
    return chs,chs_id
        

def get_neis(senders,indexs,poss):
    # i=0
    neis=[]
    for r in indexs[senders]:
        for i in np.where(r<0)[0]:
            neis.append(i)
        # neis.append(i for i in np.where(r<0)[0])
    # neis=np.where(index[senders]<0)[0]
    return neis

  

def play(env, n_round, map_size, max_steps, handles, models, print_every, eps=1.0, render=False, train=False,len_nei=40,MsgModels=None,num_bits_msg=5,msg_space=10,rewardtype='self',use_concate_mean_global=False,crps=[False,False],crp=False,selfplay=False,is_fix=False,agents_num=[0,0]):
    """play a ground and train"""
    env.reset()
    # print(crp)

    generate_map(env, map_size, handles,init_num=agents_num)

    step_ct = 0
    done = False

    n_group = len(handles)
    state = [None for _ in range(n_group)]
    acts = [None for _ in range(n_group)]
    ids = [None for _ in range(n_group)]
    sender_ids = [None for _ in range(n_group)]
    ids_=[None for _ in range(n_group)]

    alives = [None for _ in range(n_group)]
    alives_ = [None for _ in range(n_group)]
    rewards = [None for _ in range(n_group)]
    sender_rewards = [None for _ in range(n_group)]
    nums = [env.get_num(handle) for handle in handles]
    # print(nums)
    max_nums = nums.copy()

    loss = [None for _ in range(n_group)]
    eval_q = [None for _ in range(n_group)]
    n_action = [env.get_action_space(handles[0])[0], env.get_action_space(handles[1])[0]]

    # print("\n\n[*] ROUND #{0}, EPS: {1:.2f} NUMBER: {2}".format(n_round, eps, nums))
    mean_rewards = [[] for _ in range(n_group)]
    total_rewards = [[] for _ in range(n_group)]


    for i in range(n_group):
        state[i] = list(env.get_observation(handles[i]))
    msg=[[],[]]
    l_msg=[[],[]]
    msgs=[[],[]]
    mean_msg=[[],[]]
    senders=[[],[]]
    feedbacks=[[],[]]
    sender_rewards_adv=[[],[]]
    sender_rewards_all=[[],[]]
    state_senders=[[],[]]
    chs=[[],[]]
    chs_id=[[],[]]
    hps=[[],[]]
    posss=[[],[]]
    indexss=[[],[]]
    images = []
    q_values=None
    q_values_pre=None
    # index0=[]
    first=True

    models[0].test_buffer()
    models[1].test_buffer()
    if MsgModels[0]:
        MsgModels[0].test_buffer()
    while not done and step_ct < max_steps:
        
        # take actions for every model
        # start=time.time() 
        for i in range(n_group):
            state[i] = list(env.get_observation(handles[i]))
            ids[i] = env.get_agent_id(handles[i])   
            hps[i]=env.get_hps(handles[i])   
            posss[i] = env.get_pos(handles[i])
            indexss[i],poss=get_distance_matrix(env,handles[i],len_nei)
            
            if MsgModels[i] :
                # pdb.set_trace()
                # print('msgmodel msgmodel.....')
                # if i==0:
                #     index0=indexs
                msg[i]=MsgModels[i].act(obs=state[i][0],eps=eps)

            chs[i],chs_id[i]=get_chs(indexss[i],hps[i],chs_id[i],ids[i],msg[i])
            # print(len(chs[0]))
            # pdb.set_trace()
            if crp:
                if i==1:
                    acts[i] = models[i].act(obs=state[i][0], eps=eps, msgs=msgs[i], mean_msg=mean_msg[i],
                                    crp=crp,poss=posss[i],chs=chs[i])
                else:
                    acts[i] = models[i].act(obs=state[i][0], eps=eps, msgs=msgs[i], mean_msg=mean_msg[i],
                                    crp=crp,poss=posss[i],chs=chs[i])
            else:
                if i==1:
                    acts[i] = models[i].act(obs=state[i][0], eps=eps, msgs=msgs[i], mean_msg=mean_msg[i],
                                        crp=crp,poss=posss[i]) 
                else:
                    acts[i] = models[i].act(obs=state[i][0], eps=eps, msgs=msgs[i], mean_msg=mean_msg[i],
                                        crp=crp,poss=posss[i]) 
#                e_res=models[i].get_e_res(state=state[i], eps=eps, msgs=msgs[i], mean_msg=mean_msg[i],
#                                        crp=crp,poss=posss[i])
                # cluster
                                    # poss=np.array([posss[i] for j in range(len(state[i][0]))]))
            # print('nn',time.time()-s0)
        # print('get acts',time.time()-start)

        #pdb.set_trace() 
        for i in range(n_group):
            env.set_action(handles[i], acts[i])

        # simulate one step
        done = env.step()
        sender_alives=[[],[]]
        
        for i in range(n_group):
            rewards[i] = env.get_reward(handles[i])
            alives[i] = env.get_alive(handles[i])
            # sender_rewards[i]=rewards[i][senders[i]]
            
        if not first:
            pre_reward=buffer['rewards']
           
        buf_chs=[None for i in range(len(posss[0]))]
        i=0
        for p in posss[0]:
            buf_chs[i]=chs[0][i%len(chs[0])]
            i+=1
        # print(chs[0])
        buffer = {
            'state': state[0], 'acts': acts[0],
            'alives': alives[0], 'ids': ids[0], 'poss': posss[0],'chs':chs[0],
        }
        buffer1 = {
            'state': state[1], 'acts': acts[1],
            'alives': alives[1], 'ids': ids[1], 'poss': posss[1],'chs':chs[1],
        }
        
        # print(buffer['alives'])
        buffer['msgs'] = msgs[0]
        buffer['mean_msg']=mean_msg[0]
        buffer['rewards']=rewards[0]
        
        buffer1['msgs'] = msgs[1]
        buffer1['mean_msg']=mean_msg[1]
        buffer1['rewards']=rewards[1]
        
            
        # print(buffer)
         
        if not first:
            buffer_=buffer.copy()
            buffer_['obs']=buffer_['state'][0]
            buffer_['feature']=buffer_['state'][1]
            r_indexs=[]
            for j in range(len(ids[0])):
                if ids[0][j] in ids_[0]:
                    r_indexs.append(j)
            buffer_['rewards'] = pre_reward[r_indexs]
            buffer_['dones']=np.array([not done for done in buffer_['alives']], dtype=np.bool)

            # print('----------------------------')
            # print(len(buffer_['dones']))
            # print(len(buffer_['rewards']))
            # print(len(buffer_['obs']))

        #     q_values_pre=models[0].calc_target_q(**buffer_)
        # # if not first:
        #     feedbacks[0]=q_values[r_indexs]-q_values_pre       
            
            sender_rewards_adv[0]=np.array([sum(feedbacks[0])/len(rewards[0]) for k in range(len(senders[0]))])
        # q_values=models[0].get_q(**buffer)
        # sender_rewards_all[0]=np.array([sum(rewards[0])/len(rewards[0]) for k in range(len(senders[0]))])
        sender_alives[0]=alives[0][senders[0]]
        
        # if use_concate_mean_global:
        #     msg[0]=l_msg[0][senders[0]]
        if rewardtype == 'self':
            buffer_msg = {
                    'state': state[0], 'acts': msg[0], 'rewards': rewards[0],
                    'alives': alives[0], 'ids': ids[0]
            }
        
        for i in range(n_group):
            ids_[i]=ids[i]
            # act_onehot=np.eye(n_action[i])[np.array(acts[i])]
            indexs,poss=get_distance_matrix(env,handles[i],len_nei)

        # pdb.set_trace()
        
        if first:
            first=False
        # stat info
        
        if render:
            # obs = env.get_observation()  # 获取观察值
            # obs_view = obs["view"]  # 获取视图观察值
            # image = np.array(obs_view)  # 将观察值转换为NumPy数组
            # images.append(image)
            env.render()

        # clear dead agents
        env.clear_dead()
        nums = [env.get_num(handle) for handle in handles]
        
#        if done:
            
#            if nums[1]==0:
#                buffer['alives']=[False for i in range(len(alives[0]))]
#                buffer['rewards'] = [10 for i in range(len(ids[0]))]
#            elif nums[0]==0:
#                buffer['rewards']=[-10 for i in range(len(ids[0]))]

        if train:
            models[0].flush_buffer(**buffer)
            models[1].flush_buffer(**buffer1)
            if not first and MsgModels[0]:
                MsgModels[0].flush_buffer(**buffer_msg)

        for i in range(n_group):
            sum_reward = sum(rewards[i])
            if nums[i]!=0:
                rewards[i] = sum_reward/nums[i]
            else:
                rewards[i] = sum_reward 
            mean_rewards[i].append(rewards[i])
            total_rewards[i].append(sum_reward)
        info = {"Ave-Reward": np.round(rewards, decimals=6), "NUM": nums}
        step_ct += 1

        # if step_ct % print_every == 0:
        #     print("> step #{}, info: {}".format(step_ct, info))
        #     print("Run: 10, Epoch: 500, Loss: 7255.9017, Learning Rate: 0.0010, ")

    if train:
        models[0].train()
        if not  selfplay and not is_fix: 
            models[1].train()
        if MsgModels[0]:
            MsgModels[0].train()
    for i in range(n_group):
        mean_rewards[i] = sum(mean_rewards[i]) / len(mean_rewards[i])
        total_rewards[i] = sum(total_rewards[i])
    return max_nums, nums, mean_rewards, total_rewards, images

