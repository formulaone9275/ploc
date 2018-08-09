from __future__ import print_function
import tensorflow as tf
import numpy as np
import pickle
from collections import defaultdict
from constant import *
from CNN_Model_Yifan import CNNModel
from glob import glob
from sklearn.metrics import precision_recall_curve
from sklearn.metrics import average_precision_score
from build_tfrecord_data_aimed import load_context_matrix, load_tagged, create_tensor,load_sentence_matrix_v1
from build_tfrecord_data_new import load_sentence_matrix_v1_new
import random
from word_embedding import VOCAB_SIZE,EMBEDDING

def finetune_cont_model(train_data,dev_data,test_data,folder_name,cv_i):
    with tf.Graph().as_default():
        model = CNNModel()
        model.build_graph(True)
        saver = tf.train.Saver(tf.global_variables(),max_to_keep=None)

        with tf.Session() as sess:
            init = tf.global_variables_initializer()
            init_l = tf.local_variables_initializer()
            sess.run([init, init_l])
            if tf.flags.FLAGS.open_transfer is True:
                saver.restore(sess, 'model_cnn_new/'+folder_name+'/pcnn_model_10')
            sess.run(tf.assign(model.embedding_weights, EMBEDDING))
            sess.run(tf.assign(model.entity_type_weights, [[0, 0, 0, 0],
                                                           [1, 0, 0, 0],
                                                           [0, 1, 0, 0],
                                                           [0, 0, 1, 0],
                                                           [0, 0, 0, 1],
                                                           [0, 0, 1, 1]]))
            sess.run(tf.assign(model.global_step, 0))
            lr = sess.run(model.learning_rate)
            #print(lr)
            sent, sent_len, head, head_len, dep, dep_len, entity, labels = train_data

            feed_dict={
                model.ph_sent: sent,
                model.ph_head: head,
                model.label_placeholder: labels,
                model.drop_rate: tf.flags.FLAGS.drop_rate,
                model.drop_rate_dense: tf.flags.FLAGS.drop_rate_dense,
                model.is_training: True,
            }

            p, r, f, pred = sess.run(
                [model.precision, model.recall, model.fscore, model.pred],
                feed_dict=feed_dict)
            print(p,r,f)

            data_size = len(sent)
            batch_size = 128
            batch_num = int(data_size / batch_size)
            dev_f_score_step=[]
            dev_precision_step=[]
            dev_recall_step=[]
            dev_pred_step=[]
            test_f_score_step=[]
            test_precision_step=[]
            test_recall_step=[]
            for epoch in range(250):
                shuf = np.random.permutation(np.arange(data_size))
                sent = sent[shuf]
                head = head[shuf]

                labels = labels[shuf]

                for batch in range(batch_num):
                    batch_start = batch * batch_size
                    batch_end = (batch + 1) * batch_size
                    if batch_end > data_size:
                        batch_end = data_size

                    b_sent = sent[batch_start:batch_end]
                    b_head = head[batch_start:batch_end]

                    b_label = labels[batch_start:batch_end]

                    feed_dict={
                        model.ph_sent: b_sent,
                        model.ph_head: b_head,
                        model.label_placeholder: b_label,
                        model.drop_rate: tf.flags.FLAGS.drop_rate,
                        model.drop_rate_dense: tf.flags.FLAGS.drop_rate_dense,
                        model.is_training: True,
                    }

                    _, mini_loss = sess.run([model.train_op, model.loss],
                                            feed_dict=feed_dict)
                #see the performance on development set
                if (epoch+1)%25==0:
                    sent_dev, sent_len_dev, head_dev, head_len_dev, dep_dev, dep_len_dev, entity_dev, labels_dev = dev_data
                    feed_dict={
                        model.ph_sent: sent_dev,
                        model.ph_head: head_dev,
                        model.label_placeholder: labels_dev,
                        model.drop_rate: tf.flags.FLAGS.drop_rate,
                        model.drop_rate_dense: tf.flags.FLAGS.drop_rate_dense,
                        model.is_training: True,
                    }

                    p, r, f, pred = sess.run(
                        [model.precision, model.recall, model.fscore, model.pred],
                        feed_dict=feed_dict)
                    print(epoch+1,'performance on dev set',p,r,f)
                    if tf.flags.FLAGS.save_model:
                        global_step = sess.run(model.global_step)
                        path = saver.save(sess, 'model_new/'+folder_name+'/'+tf.flags.FLAGS.name+'cv'+str(cv_i)+'_transfer_step'+str(epoch+1))
                    dev_f_score_step.append(f)
                    dev_precision_step.append(p)
                    dev_recall_step.append(r)
                    dev_pred_step.append(pred)


                    #see the test set performance

                    sent_test, sent_len_test, head_test, head_len_test, dep_test, dep_len_test, entity_test, labels_test = test_data
                    feed_dict={
                        model.ph_sent: sent_test,
                        model.ph_head: head_test,
                        model.label_placeholder: labels_test,
                        model.drop_rate: tf.flags.FLAGS.drop_rate,
                        model.drop_rate_dense: tf.flags.FLAGS.drop_rate_dense,
                        model.is_training: False,
                    }

                    p, r, f, pred = sess.run(
                        [model.precision, model.recall, model.fscore, model.pred],
                        feed_dict=feed_dict)
                    print('Performance on test set',p,r,f)
                    test_f_score_step.append(f)
                    test_precision_step.append(p)
                    test_recall_step.append(r)

            del model
            print('epoch {}, batch {}, loss {}'.format(epoch, batch, mini_loss))
            return dev_f_score_step,test_f_score_step,dev_precision_step,test_precision_step,dev_recall_step,test_recall_step,dev_pred_step

def get_train_test_file(filename,cv,test=True):

    files_all = glob(filename)
    if test:
        files = glob(files_all[cv])

    else:
        del files_all[cv]
        files=files_all
    return files

if __name__ == '__main__':

    cv=10

    file_name='./data/model_document/fold*.txt'
    model_folder='model_cnn_new/'
    hyperparameter='mc_2'
    html_dev_file='data/html/html_dev.txt'
    for folder_name in ['baseline','CP','CP_TW','CP_HP']: #,'CP','CP_TW','CP_HP','filtered'

        if tf.flags.FLAGS.model == 'pcnn':

            cv_prec, cv_recall, cv_fscore = [], [], []
            cv_prec_dev, cv_recall_dev,cv_fscore_dev=[],[],[]
            cv_pred_dev=[]
            for ii in range(cv):
                file_to_train=get_train_test_file(file_name,ii,False)
                random.shuffle(file_to_train)

                #the first 8 files to train, the last file for development
                sent_mx_train, head_mx_train, dep_mx_train, entity_mx_train,label_train=[],[],[],[],[]
                for jj in range(len(file_to_train)):
                    sent_mx, head_mx,dep_mx, entity_mx, labels=[],[],[],[],[]
                    sent_mx, head_mx,dep_mx, entity_mx, labels = load_sentence_matrix_v1(file_to_train[jj])

                    sent_mx_train+=sent_mx
                    head_mx_train+=head_mx
                    dep_mx_train+=dep_mx
                    entity_mx_train+=entity_mx

                    label_train+=list(labels)

                sent, sent_len = create_tensor(sent_mx_train,160, [0] * 6)
                head, head_len = create_tensor(head_mx_train,160, [0] * 6)

                dep, dep_len = create_tensor(dep_mx_train,20, [0] * 6)
                entity = np.array(entity_mx_train)
                train_data = [sent, sent_len, head, head_len, dep, dep_len, entity,np.array(label_train)]

                #get development data
                sent_mx_dev, head_mx_dev, dep_mx_dev, entity_mx_dev,labels_dev=[],[],[],[],[]
                sent_mx_dev, head_mx_dev, dep_mx_dev, entity_mx_dev,labels_dev = load_sentence_matrix_v1(file_to_train[-1])

                sent_dev, sent_len_dev = create_tensor(sent_mx_dev,160, [0] * 6)
                head_dev, head_len_dev = create_tensor(head_mx_dev,160, [0] * 6)

                dep_dev, dep_len_dev = create_tensor(dep_mx_dev,20, [0] * 6)
                entity_dev = np.array(entity_mx_dev)
                dev_data = [sent_dev, sent_len_dev, head_dev, head_len_dev, dep_dev, dep_len_dev, entity_dev, np.array(labels_dev)]

                #get the test set
                file_to_test=get_train_test_file(file_name,ii,True)
                sent_mx_test, head_mx_test, dep_mx_test, entity_mx_test,labels_test=[],[],[],[],[]
                sent_mx_test, head_mx_test, dep_mx_test, entity_mx_test,labels_test = load_sentence_matrix_v1(file_to_test[-1])

                sent_test, sent_len_test = create_tensor(sent_mx_test,160, [0] * 6)
                head_test, head_len_test = create_tensor(head_mx_test,160, [0] * 6)

                dep_test, dep_len_test = create_tensor(dep_mx_test,20, [0] * 6)
                entity_test = np.array(entity_mx_test)
                test_data = [sent_test, sent_len_test, head_test, head_len_test, dep_test, dep_len_test, entity_test, np.array(labels_test)]

                #train the model and
                dev_f,test_f,dev_p,test_p,dev_r,test_r,dev_pred = finetune_cont_model(train_data, dev_data,test_data,folder_name,ii)
                #cv_prec.append(p)
                #cv_recall.append(r)
                cv_fscore.append(test_f)
                cv_fscore_dev.append(dev_f)
                cv_prec.append(test_p)
                cv_prec_dev.append(dev_p)
                cv_recall.append(test_r)
                cv_recall_dev.append(dev_r)
                cv_pred_dev.append(dev_pred)
            #print(sum(cv_prec)/len(cv_prec), sum(cv_recall)/len(cv_recall), sum(cv_fscore)/len(cv_fscore))
            with open(model_folder+folder_name+'/test_fscore_doc_new_'+hyperparameter+'.pickle', 'wb') as f:
                # Pickle the 'data' dictionary using the highest protocol available.
                pickle.dump(cv_fscore, f, pickle.HIGHEST_PROTOCOL)
            with open(model_folder+folder_name+'/dev_fscore_doc_new_'+hyperparameter+'.pickle', 'wb') as f:
                # Pickle the 'data' dictionary using the highest protocol available.
                pickle.dump(cv_fscore_dev, f, pickle.HIGHEST_PROTOCOL)


            with open(model_folder+folder_name+'/test_precision_doc_new_'+hyperparameter+'.pickle', 'wb') as f:
                # Pickle the 'data' dictionary using the highest protocol available.
                pickle.dump(cv_prec, f, pickle.HIGHEST_PROTOCOL)
            with open(model_folder+folder_name+'/dev_precision_doc_new_'+hyperparameter+'.pickle', 'wb') as f:
                # Pickle the 'data' dictionary using the highest protocol available.
                pickle.dump(cv_prec_dev, f, pickle.HIGHEST_PROTOCOL)


            with open(model_folder+folder_name+'/test_recall_doc_new_'+hyperparameter+'.pickle', 'wb') as f:
                # Pickle the 'data' dictionary using the highest protocol available.
                pickle.dump(cv_recall, f, pickle.HIGHEST_PROTOCOL)
            with open(model_folder+folder_name+'/dev_recall_doc_new_'+hyperparameter+'.pickle', 'wb') as f:
                # Pickle the 'data' dictionary using the highest protocol available.
                pickle.dump(cv_recall_dev, f, pickle.HIGHEST_PROTOCOL)

            with open(model_folder+folder_name+'/dev_pred_doc_new_'+hyperparameter+'.pickle', 'wb') as f:
                # Pickle the 'data' dictionary using the highest protocol available.
                pickle.dump(cv_pred_dev, f, pickle.HIGHEST_PROTOCOL)




