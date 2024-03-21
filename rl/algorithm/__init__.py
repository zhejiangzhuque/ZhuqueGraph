from . import q_learning

IL = q_learning.DQN
MSGDQN = q_learning.MsgDQN
GIL=q_learning.GDQN

def spawn_ai(algo_name, optimizer, sess, env, handle, human_name, max_steps,len_nei=6):
    if algo_name == 'il':
        model = IL(optimizer,sess, human_name, handle, env, max_steps, len_nei,memory_size=80)
    elif algo_name == 'gil':
        model = GIL(optimizer,sess, human_name, handle, env, max_steps,len_nei, memory_size=8,isHie=False)
    elif algo_name == 'msgdqn':     
        model=MSGDQN(optimizer,sess,human_name,handle,env,max_steps,len_nei,msg_bits=3,memory_size=80000)

    return model
