import visdom


def open_win(env_name):
    global vis
    global e_name
    e_name = env_name
    vis = visdom.Visdom(env=e_name)


def init_line(win_name='train_loss', plot_name='train_loss', start_x=0., start_y=0.):
    vis.line([start_y],  # Y的第一个点的坐标
             [start_x],  # X的第一个点的坐标
             win=win_name,  # 画图窗口的名称
             opts=dict(title=plot_name)  # 图像的标例
             )
    vis.save([e_name])


def append_line(x, y, win_name='train_loss'):
    vis.line([y], [x], win=win_name, update='append')
    vis.save([e_name])
