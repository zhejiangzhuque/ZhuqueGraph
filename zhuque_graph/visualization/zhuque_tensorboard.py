from torch.utils.tensorboard import SummaryWriter


def open_win(filename_suffix):
    global writer
    writer = SummaryWriter('/OUTPUT/tensorboard/' + filename_suffix + '/', filename_suffix='.' + filename_suffix)


def draw_line(x, y, tag_name='loss'):
    writer.add_scalar(tag_name, y, global_step=x)
    writer.flush()
