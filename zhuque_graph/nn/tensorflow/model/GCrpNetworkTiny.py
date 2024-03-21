from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.abspath(__file__),"../../layer/")))
sys.path.append(os.path.abspath(os.path.join(os.path.abspath(__file__),"../../../utils/")))
import modules
import rl_tf_utils
import sonnet as snt
from rl_sonnet_utils import *

class GCrpNetworkTiny(snt.AbstractModule):
  """GraphNetwork with MLP edge, node, and global models."""

  def __init__(self, name="GCrpNetworkTiny"):
    super(GCrpNetworkTiny, self).__init__(name=name)
    with self._enter_variable_scope():
      self._obsEncoder=modules.obsEncoder(encoder_fn=make_conv_model)
      self._network = modules.CommNet(
          edge_model_fn=make_edge_model,
          node_model_fn=make_node_model)

      self._hnetwork=modules.HCommNet(
         edge_model_fn=make_edge_model,
          node_model_fn=make_Hnode_model)
  
      self._Lnetwork = modules.LCommNet(
          edge_model_fn=make_edge_model,
        node_model_fn=make_Lnode_model)
      self._qnet=modules.qEncoder(mlp_fn=get_q_model2)


  def _build(self, inputs):
    
    #return self._qnet(self._obsEncoder(inputs))
    return self._qnet(self._Lnetwork(self._hnetwork(self._network(self._obsEncoder(inputs)))))
    # return