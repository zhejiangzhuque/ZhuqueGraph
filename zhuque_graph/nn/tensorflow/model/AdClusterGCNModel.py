# coding=utf-8
# Copyright 2019 The Google Research Authors.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""Collections of different Models."""

from zhuque_graph.nn.tensorflow.layer.AdClusterGCNLayer import *
import tensorflow.compat.v1 as tf


FLAGS = tf.flags.FLAGS


def masked_softmax_cross_entropy(preds, labels, mask):
    """Softmax cross-entropy loss with masking."""
    loss = tf.nn.softmax_cross_entropy_with_logits(logits=preds, labels=labels)
    mask = tf.cast(mask, dtype=tf.float32)
    mask /= tf.reduce_mean(mask)
    loss *= mask
    return tf.reduce_mean(loss)


def masked_sigmoid_cross_entropy(preds, labels, mask):
    """Sigmoid cross-entropy loss with masking."""
    loss_all = tf.nn.sigmoid_cross_entropy_with_logits(
        logits=preds, labels=labels)
    mask = tf.cast(mask, dtype=tf.float32)
    mask /= tf.reduce_mean(mask)
    loss = tf.multiply(loss_all, mask[:, tf.newaxis])
    return tf.reduce_mean(loss)


def masked_accuracy(preds, labels, mask):
    """Accuracy with masking."""
    correct_prediction = tf.equal(tf.argmax(preds, 1), tf.argmax(labels, 1))
    accuracy_all = tf.cast(correct_prediction, tf.float32)
    mask = tf.cast(mask, dtype=tf.float32)
    mask /= tf.reduce_mean(mask)
    accuracy_all *= mask
    return tf.reduce_mean(accuracy_all)


def masked_accuracy_multilabel(preds, labels, mask):
    """Multilabel accuracy with masking."""
    preds = preds > 0
    labels = labels > 0.5
    correct_prediction = tf.equal(preds, labels)
    accuracy_all = tf.cast(correct_prediction, tf.float32)
    mask = tf.cast(mask, dtype=tf.float32)
    mask /= tf.reduce_mean(mask)
    accuracy_all = tf.multiply(accuracy_all, mask[:, tf.newaxis])
    return tf.reduce_mean(accuracy_all)


def masked_softmax_cross_entropy_w(preds, labels, mask, weight):
    """Softmax cross-entropy loss with masking."""
    loss = tf.nn.softmax_cross_entropy_with_logits(
        logits=preds, labels=labels) * weight
    mask = tf.cast(mask, dtype=tf.float32)
    mask /= tf.reduce_mean(mask)
    loss *= mask
    return tf.reduce_mean(loss)


def masked_sigmoid_cross_entropy_w(preds, labels, mask, weight):
    """Sigmoid cross-entropy loss with masking."""
    loss_all = tf.nn.sigmoid_cross_entropy_with_logits(
        logits=preds, labels=labels)
    mask = tf.cast(mask, dtype=tf.float32)
    mask /= tf.reduce_mean(mask)
    loss = tf.multiply(tf.multiply(loss_all, mask[:, tf.newaxis]),weight[:,tf.newaxis])
    return tf.reduce_mean(loss)


def masked_new_weight(preds, labels, mask, weight,acc):
    weight = tf.cast(weight, tf.float32)
    mask = tf.cast(mask, tf.float32)
    correct_prediction = tf.cast(tf.equal(tf.argmax(preds, 1), tf.argmax(labels, 1)), tf.float32)
    #i_weight =tf.cast(tf.clip_by_value(weight*mask*tf.exp(tf.cast((1-2*correct_prediction)*acc/(1.0-acc+1e-6),tf.float32)) +1.0-mask,0.0,1000.0) , tf.float32)
    i_weight =tf.cast(tf.clip_by_value(weight*mask*tf.exp(tf.cast((1-2*correct_prediction)*tf.log(tf.sqrt(acc/(1.0-acc+1e-6))),tf.float32)) +1.0-mask,0.0,1000.0) , tf.float32)
    return i_weight/ tf.reduce_mean(i_weight)
    # makes the average weight of nodes is 1

def masked_new_weight_multilabel(preds,labels,mask,weight,acc):
    weight =tf.cast(weight,tf.float32)
    mask = tf.cast(mask,tf.float32)
    preds = preds > 0
    labels = labels > 0.5
    correct_prediction = tf.equal(preds, labels)
    accuracy_all = tf.reduce_mean(tf.cast(correct_prediction, tf.float32),axis=-1)
    i_weight =tf.cast(weight*mask*tf.cast((1-2*accuracy_all)*acc/(1.0-acc+1e-6),tf.float32) +1.0-mask , tf.float32)
    i_weight /=tf.reduce_mean(i_weight)
    return i_weight


class Model(object):
    """Model class to be inherited."""

    def __init__(self, **kwargs):
        allowed_kwargs = {
            'name', 'logging', 'multilabel', 'norm', 'precalc', 'num_layers'
        }
        for kwarg, _ in kwargs.items():
            assert kwarg in allowed_kwargs, 'Invalid keyword argument: ' + kwarg
        name = kwargs.get('name')
        if not name:
            name = self.__class__.__name__.lower()
        self.name = name

        logging = kwargs.get('logging', False)
        self.logging = logging

        self.vars = {}
        self.placeholders = {}

        self.layers = []
        self.activations = []

        self.inputs = None
        self.outputs = None

        self.loss = 0
        self.accuracy = 0
        self.pred = 0
        self.optimizer = None
        self.opt_op = None
        self.multilabel = kwargs.get('multilabel', False)
        self.norm = kwargs.get('norm', False)
        self.precalc = kwargs.get('precalc', False)
        self.num_layers = kwargs.get('num_layers', 2)

    def _build(self):
        raise NotImplementedError

    def build(self):
        """Wrapper for _build()."""
        with tf.variable_scope(self.name):
            self._build()

        # Build sequential layer model
        self.activations.append(self.inputs)
        for layer in self.layers:
            hidden = layer(self.activations[-1])
            if isinstance(hidden, tuple):
                tf.logging.info('{} shape = {}'.format(layer.name,
                                                       hidden[0].get_shape()))
            else:
                tf.logging.info(
                    '{} shape = {}'.format(
                        layer.name,
                        hidden.get_shape()))
            self.activations.append(hidden)
        self.outputs = self.activations[-1]

        # Store model variables for easy access
        variables = tf.get_collection(
            tf.GraphKeys.GLOBAL_VARIABLES, scope=self.name)
        self.vars = variables
        for k in self.vars:
            tf.logging.info((k.name, k.get_shape()))

        # Build metrics
        self._loss()
        self._accuracy()
        self._predict()

        self.opt_op = self.optimizer.minimize(self.loss)

    def _loss(self):
        """Construct the loss function."""
        # Weight decay loss
        if FLAGS.weight_decay > 0.0:
            for var in self.layers[0].vars.values():
                self.loss += FLAGS.weight_decay * tf.nn.l2_loss(var)

        # Cross entropy error
        if self.multilabel:
            self.loss += masked_sigmoid_cross_entropy(
                self.outputs, self.placeholders['labels'],
                self.placeholders['labels_mask'])
        else:
            self.loss += masked_softmax_cross_entropy(
                self.outputs, self.placeholders['labels'],
                self.placeholders['labels_mask'])

    def _accuracy(self):
        if self.multilabel:
            self.accuracy = masked_accuracy_multilabel(
                self.outputs, self.placeholders['labels'],
                self.placeholders['labels_mask'])
        else:
            self.accuracy = masked_accuracy(
                self.outputs,
                self.placeholders['labels'],
                self.placeholders['labels_mask'])

    def _predict(self):
        if self.multilabel:
            self.pred = tf.nn.sigmoid(self.outputs)
        else:
            self.pred = tf.nn.softmax(self.outputs)

    def save(self, sess=None):
        if not sess:
            raise AttributeError('TensorFlow session not provided.')
        saver = tf.train.Saver(self.vars)
        save_path = saver.save(sess, 'tmp/%s.ckpt' % self.name)
        tf.logging.info('Model saved in file:', save_path)

    def load(self, sess=None):
        if not sess:
            raise AttributeError('TensorFlow session not provided.')
        saver = tf.train.Saver(self.vars)
        save_path = 'tmp/%s.ckpt' % self.name
        saver.restore(sess, save_path)
        tf.logging.info('Model restored from file:', save_path)


class GCN(Model):
    """Implementation of GCN model."""

    def __init__(self, placeholders, input_dim, **kwargs):
        super(GCN, self).__init__(**kwargs)

        self.inputs = placeholders['features']
        self.input_dim = input_dim
        self.output_dim = placeholders['labels'].get_shape().as_list()[1]
        self.placeholders = placeholders

        self.optimizer = tf.train.AdamOptimizer(
            learning_rate=FLAGS.learning_rate)

        self.build()

    def _build(self):
        # note that input_dim will not be doubled if precalc,input is Nx2d , Weights id 2d*d
        # if not precalc, the input is Nxd,in call AX+X makes it to Nx2d, and
        # Weights is 2d*d
        self.layers.append(
            GraphConvolution(
                input_dim=self.input_dim if self.precalc else self.input_dim * 2,
                output_dim=FLAGS.hidden1,
                placeholders=self.placeholders,
                act=tf.nn.relu,
                dropout=True,
                sparse_inputs=False,
                logging=self.logging,
                norm=self.norm,
                precalc=self.precalc))

        for _ in range(self.num_layers - 2):
            self.layers.append(
                GraphConvolution(
                    input_dim=FLAGS.hidden1 * 2,
                    output_dim=FLAGS.hidden1,
                    placeholders=self.placeholders,
                    act=tf.nn.relu,
                    dropout=True,
                    sparse_inputs=False,
                    logging=self.logging,
                    norm=self.norm,
                    precalc=False))

        self.layers.append(
            GraphConvolution(
                input_dim=FLAGS.hidden1 * 2,
                output_dim=self.output_dim,
                placeholders=self.placeholders,
                act=lambda x: x,
                dropout=True,
                logging=self.logging,
                norm=False,
                precalc=False))


class GCN_hc(Model):
    """Implementation of GCN model."""

    def __init__(
            self,
            placeholders,
            placeholders_hc,
            input_dim,
            num_hier_layer,
            **kwargs):
        super(GCN_hc, self).__init__(**kwargs)

        self.inputs = placeholders['features']
        self.input_dim = input_dim
        self.output_dim = placeholders['labels'].get_shape().as_list()[1]
        self.placeholders = placeholders
        self.placeholders_hc = placeholders_hc
        self.num_hier_layer = num_hier_layer
        self.hier_layers = [[] for _ in range(self.num_hier_layer)]
        self.hier_activations = [[] for _ in range(self.num_hier_layer)]
        self.hier_outputs = []
        self.optimizer = tf.train.AdamOptimizer(
            learning_rate=FLAGS.learning_rate)

        self.build()

    def _build(self):
        # We do not precalc
        # basic layers
        self.layers.append(
            GraphConvolution(
                input_dim=self.input_dim * 2 * 2,
                output_dim=FLAGS.hidden1,
                placeholders=self.placeholders,
                act=tf.nn.relu,
                dropout=True,
                sparse_inputs=False,
                logging=self.logging,
                norm=self.norm,
                precalc=False))

        for _ in range(self.num_layers - 2):
            self.layers.append(
                GraphConvolution(
                    input_dim=FLAGS.hidden1 * 2,
                    output_dim=FLAGS.hidden1,
                    placeholders=self.placeholders,
                    act=tf.nn.relu,
                    dropout=True,
                    sparse_inputs=False,
                    logging=self.logging,
                    norm=self.norm,
                    precalc=False))

        self.layers.append(
            GraphConvolution(
                input_dim=FLAGS.hidden1 * 2,
                output_dim=self.output_dim,
                placeholders=self.placeholders,
                act=lambda x: x,
                dropout=True,
                logging=self.logging,
                norm=False,
                precalc=False))

        # hier layers

        for ii in range(self.num_hier_layer):
            self.hier_layers[ii].append(
                GraphConvolution(
                    input_dim=self.input_dim * 3**(self.num_hier_layer - ii) if ii == 0 else self.input_dim * (1 + 3**(self.num_hier_layer - ii)),
                    output_dim=FLAGS.hidden1,
                    placeholders=self.placeholders_hc[ii],
                    act=tf.nn.relu,
                    dropout=True,
                    sparse_inputs=False,
                    logging=self.logging,
                    norm=self.norm,
                    precalc=True
                )
            )

            self.hier_layers[ii].append(
                GraphConvolution(
                    input_dim=FLAGS.hidden1,
                    output_dim=self.input_dim,
                    placeholders=self.placeholders_hc[ii],
                    act=tf.nn.relu,
                    dropout=True,
                    logging=self.logging,
                    norm=False,
                    precalc=True)
            )

    def build(self):
        """Wrapper for _build()."""
        with tf.variable_scope(self.name):
            self._build()

        # Build hier layer model
        for ii in range(self.num_hier_layer):
            inputs = self.placeholders_hc[ii]['features']
            if ii > 0:
                inputs = tf.concat((inputs, self.hier_outputs[-1]), axis=1)
            self.hier_activations[ii].append(inputs)
            for layer in self.hier_layers[ii]:
                hidden = layer(self.hier_activations[ii][-1])
                if isinstance(hidden, tuple):
                    tf.logging.info(
                        '{} shape = {}'.format(
                            layer.name, hidden[0].get_shape()))
                else:
                    tf.logging.info(
                        '{} shape = {}'.format(
                            layer.name, hidden.get_shape()))
                self.hier_activations[ii].append(hidden)
            self.hier_outputs.append(tf.nn.embedding_lookup(
                self.hier_activations[ii][-1], self.placeholders_hc[ii]['lookup']))

        # Build sequential layer model
        inputs = tf.concat((self.inputs, self.hier_outputs[-1]), axis=1)
        self.activations.append(inputs)
        for layer in self.layers:
            hidden = layer(self.activations[-1])
            if isinstance(hidden, tuple):
                tf.logging.info('{} shape = {}'.format(layer.name,
                                                       hidden[0].get_shape()))
            else:
                tf.logging.info(
                    '{} shape = {}'.format(
                        layer.name,
                        hidden.get_shape()))
            self.activations.append(hidden)
        self.outputs = self.activations[-1]

        # Store model variables for easy access
        variables = tf.get_collection(
            tf.GraphKeys.GLOBAL_VARIABLES, scope=self.name)
        self.vars = variables
        for k in self.vars:
            tf.logging.info((k.name, k.get_shape()))

        # Build metrics
        self._loss()
        self._accuracy()
        self._predict()

        self.opt_op = self.optimizer.minimize(self.loss)


class GCN_w(Model):
    """Implementation of GCN model."""

    def __init__(self, placeholders, input_dim, **kwargs):
        super(GCN_w, self).__init__(**kwargs)

        self.inputs = placeholders['features']
        self.input_dim = input_dim
        self.output_dim = placeholders['labels'].get_shape().as_list()[1]
        self.placeholders = placeholders

        self.optimizer = tf.train.AdamOptimizer(
            learning_rate=FLAGS.learning_rate)
        #self.loss_w = 0.0
        self.nodes_weight_new = None

        self.build()
        # self._loss_w()
        self.opt_op_w = self.optimizer.minimize(self.loss)
        self._nodes_weight_new()

    def _build(self):
        # note that input_dim will not be doubled if precalc,input is Nx2d , Weights id 2d*d
        # if not precalc, the input is Nxd,in call AX+X makes it to Nx2d, and
        # Weights is 2d*d
        self.layers.append(
            GraphConvolution(
                input_dim=self.input_dim if self.precalc else self.input_dim * 2,
                output_dim=FLAGS.hidden1,
                placeholders=self.placeholders,
                act=tf.nn.relu,
                dropout=True,
                sparse_inputs=False,
                logging=self.logging,
                norm=self.norm,
                precalc=self.precalc))

        for _ in range(self.num_layers - 2):
            self.layers.append(
                GraphConvolution(
                    input_dim=FLAGS.hidden1 * 2,
                    output_dim=FLAGS.hidden1,
                    placeholders=self.placeholders,
                    act=tf.nn.relu,
                    dropout=True,
                    sparse_inputs=False,
                    logging=self.logging,
                    norm=self.norm,
                    precalc=False))

        self.layers.append(
            GraphConvolution(
                input_dim=FLAGS.hidden1 * 2,
                output_dim=self.output_dim,
                placeholders=self.placeholders,
                act=lambda x: x,
                dropout=True,
                logging=self.logging,
                norm=False,
                precalc=False))

    def _loss_w(self):
        """Construct the loss function."""
        # Weight decay loss
        if FLAGS.weight_decay > 0.0:
            for var in self.layers[0].vars.values():
                self.loss += FLAGS.weight_decay * tf.nn.l2_loss(var)

        # Cross entropy error
        if self.multilabel:
            self.loss += masked_sigmoid_cross_entropy_w(
                self.outputs, self.placeholders['labels'],
                self.placeholders['labels_mask'],
                self.placeholders['nodes_weight'])
        else:
            self.loss += masked_softmax_cross_entropy_w(
                self.outputs, self.placeholders['labels'],
                self.placeholders['labels_mask'],
                self.placeholders['nodes_weight'])

    def _loss(self):
        """Construct the loss function."""
        # Weight decay loss
        if FLAGS.weight_decay > 0.0:
            for var in self.layers[0].vars.values():
                self.loss += FLAGS.weight_decay * tf.nn.l2_loss(var)

        # Cross entropy error
        if self.multilabel:
            self.loss += masked_sigmoid_cross_entropy_w(
                self.outputs, self.placeholders['labels'],
                self.placeholders['labels_mask'],
                self.placeholders['nodes_weight'])
        else:
            self.loss += masked_softmax_cross_entropy_w(
                self.outputs, self.placeholders['labels'],
                self.placeholders['labels_mask'],
                self.placeholders['nodes_weight'])

    def _nodes_weight_new(self):
        if self.multilabel:
            i_weight = masked_new_weight_multilabel(
                self.outputs,
                self.placeholders['labels'],
                self.placeholders['labels_mask'],
                self.placeholders['nodes_weight'],
                self.accuracy)
        else:
            i_weight = masked_new_weight(
                self.outputs,
                self.placeholders['labels'],
                self.placeholders['labels_mask'],
                self.placeholders['nodes_weight'],
                self.accuracy)

        for i in range(FLAGS.PageRank_round):
            i_weight = tf.reshape(i_weight, [-1, 1])
            i_weight = (1 - FLAGS.PageRank_alpha) * dot(self.placeholders['support'],i_weight, True) + FLAGS.PageRank_alpha * i_weight

        # i_weight /= tf.reduce_mean(i_weight)
        #i_weight = tf.reshape(i_weight, [-1, 1])

        self.nodes_weight_new = i_weight
        self.nodes_weight_new = tf.reshape(i_weight/tf.reduce_mean(i_weight),[-1,1])

        #self.nodes_weight_new = tf.reshape(self.placeholders['nodes_weight'],[-1,1])

class GCN_adclick(Model):
    """Implementation of GCN model."""

    def __init__(self, placeholders, input_dim, total_user,total_item,id_dim,**kwargs):
        super(GCN_adclick, self).__init__(**kwargs)

        self.user_raw_features = placeholders['user_raw_features']
        self.item_raw_features = placeholders['item_raw_features']
        self.user_id_lists = placeholders['user_id_lists']
        self.item_id_lists = placeholders['item_id_lists']

        self.total_user_id_features = [glorot(shape=[total_user,id_dim]) for _ in range(len(self.user_id_lists))]
        self.total_item_id_features = [glorot(shape=[total_item,id_dim]) for _ in range(len(self.item_id_lists))]

        self.user_features = [tf.nn.embedding_lookup(self.total_user_id_features[ii],self.user_id_lists[ii]) for ii in range(len(self.user_id_lists))]
        self.user_features.append(self.user_raw_features)
        self.item_features = [tf.nn.embedding_lookup(self.total_item_id_features[ii],self.item_id_lists[ii]) for ii in range(len(self.item_id_lists))]
        self.item_features.append(self.item_raw_features)

        self.user_transform = glorot(shape=[id_dim*len(self.user_id_lists)+self.user_raw_features.shape[1],input_dim])
        self.item_transform = glorot(shape=[id_dim*len(self.item_id_lists)+self.item_raw_features.shape[1],input_dim])

        self.inputs =tf.concat([tf.matmul(tf.concat(self.user_features,axis=1),self.user_transform),tf.matmul(tf.concat(self.item_features,axis=1),self.item_transform)],axis=0)

        self.input_dim = input_dim
        self.output_dim = placeholders['labels'].get_shape().as_list()[1]
        self.placeholders = placeholders

        self.optimizer = tf.train.AdamOptimizer(learning_rate=FLAGS.learning_rate)
        #self.loss_w = 0.0
        self.nodes_weight_new = None

        self.build()
        # self._loss_w()
        self.opt_op_w = self.optimizer.minimize(self.loss)
        self._nodes_weight_new()

    def _build(self):
        # note that input_dim will not be doubled if precalc,input is Nx2d , Weights id 2d*d
        # if not precalc, the input is Nxd,in call AX+X makes it to Nx2d, and
        # Weights is 2d*d
        self.layers.append(
            GraphConvolution(
                input_dim=self.input_dim if self.precalc else self.input_dim * 2,
                output_dim=FLAGS.hidden1,
                placeholders=self.placeholders,
                act=tf.nn.relu,
                dropout=True,
                sparse_inputs=False,
                logging=self.logging,
                norm=self.norm,
                precalc=self.precalc))

        for _ in range(self.num_layers - 2):
            self.layers.append(
                GraphConvolution(
                    input_dim=FLAGS.hidden1 * 2,
                    output_dim=FLAGS.hidden1,
                    placeholders=self.placeholders,
                    act=tf.nn.relu,
                    dropout=True,
                    sparse_inputs=False,
                    logging=self.logging,
                    norm=self.norm,
                    precalc=False))

        self.layers.append(
            GraphConvolution(
                input_dim=FLAGS.hidden1 * 2,
                output_dim=self.output_dim,
                placeholders=self.placeholders,
                act=lambda x: x,
                dropout=True,
                logging=self.logging,
                norm=False,
                precalc=False))

    def _loss_w(self):
        """Construct the loss function."""
        # Weight decay loss
        if FLAGS.weight_decay > 0.0:
            for var in self.layers[0].vars.values():
                self.loss += FLAGS.weight_decay * tf.nn.l2_loss(var)

        # Cross entropy error
        if self.multilabel:
            self.loss += masked_sigmoid_cross_entropy_w(
                self.outputs, self.placeholders['labels'],
                self.placeholders['labels_mask'],
                self.placeholders['nodes_weight'])
        else:
            self.loss += masked_softmax_cross_entropy_w(
                self.outputs, self.placeholders['labels'],
                self.placeholders['labels_mask'],
                self.placeholders['nodes_weight'])

    def _loss(self):
        """Construct the loss function."""
        # Weight decay loss
        if FLAGS.weight_decay > 0.0:
            for var in self.layers[0].vars.values():
                self.loss += FLAGS.weight_decay * tf.nn.l2_loss(var)

        # Cross entropy error
        if self.multilabel:
            self.loss += masked_sigmoid_cross_entropy_w(
                self.outputs, self.placeholders['labels'],
                self.placeholders['labels_mask'],
                self.placeholders['nodes_weight'])
        else:
            self.loss += masked_softmax_cross_entropy_w(
                self.outputs, self.placeholders['labels'],
                self.placeholders['labels_mask'],
                self.placeholders['nodes_weight'])

    def _nodes_weight_new(self):
        if self.multilabel:
            i_weight = masked_new_weight_multilabel(
                self.outputs,
                self.placeholders['labels'],
                self.placeholders['labels_mask'],
                self.placeholders['nodes_weight'],
                self.accuracy)
        else:
            i_weight = masked_new_weight(
                self.outputs,
                self.placeholders['labels'],
                self.placeholders['labels_mask'],
                self.placeholders['nodes_weight'],
                self.accuracy)
        for i in range(FLAGS.PageRank_round):
            i_weight = tf.reshape(i_weight, [-1, 1])
            i_weight = (1 - FLAGS.PageRank_alpha) * dot(self.placeholders['support'],i_weight, True) + FLAGS.PageRank_alpha * i_weight
        i_weight /= tf.reduce_mean(i_weight)
        i_weight = tf.reshape(i_weight, [-1, 1])
        self.nodes_weight_new = i_weight

        self.nodes_weight_new = tf.reshape(i_weight/tf.reduce_mean(i_weight),[-1,1])
