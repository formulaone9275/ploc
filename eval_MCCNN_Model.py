from __future__ import print_function
import tensorflow as tf
import numpy as np
import pickle
from constant import *
from CNN_Model_Yifan import CNNModel

from build_tfrecord_data_aimed import load_sentence_matrix, load_tagged, create_tensor

def evaluate_cont_model(eval_data,model_path):
    #map_index = index_2_word()
    with tf.Graph().as_default():
        model = CNNModel()
        model.build_graph()
        saver = tf.train.Saver(tf.global_variables())
        with tf.Session() as sess:
            saver.restore(sess, model_path)
            sent_mx, head_mx, dep, entity, label = eval_data
            feed_dict={
                model.ph_sent: sent_mx,
                model.ph_head: head_mx,
                model.label_placeholder: labels,
                model.drop_rate: 0,
                model.drop_rate_dense: 0,
                model.is_training: False,
            }
            p, r, f, prob, l = sess.run(
                [model.precision, model.recall, model.fscore, model.prob, model.loss],
                feed_dict=feed_dict)
            #print(prob)
            print(p, r, f, l)


if __name__ == '__main__':

    #for folder_name in ['baseline','CP','CP_TW','filtered']:
    for folder_name in ['baseline','CP','CP_TW','CP_HP','filtered']: #'baseline','CP','CP_TW','filtered'
        for model_path in ['model_cnn_new/'+folder_name+'/pcnn_model_10']:
            print(model_path)

            test_data, entity, labels = load_sentence_matrix('data/aimed.txt')
            sent,head_sent, dep = test_data
            sent_mx, sent_len = create_tensor(*sent)
            head_mx, head_sent_len = create_tensor(*head_sent)

            dep_mx, dep_len = create_tensor(*dep)
            eval_data = [sent_mx, head_mx, dep_mx, entity, labels]
            evaluate_cont_model(eval_data,model_path)





































