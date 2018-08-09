from __future__ import print_function
from pymongo import MongoClient
import json
import codecs
import os
from pprint import pprint


json_path='/usa/psu/Documents/ploc/data/LocText/'
filepath='/usa/psu/Documents/ploc/data/'
outputfile='loctext'
file_list=os.listdir(json_path)
count=0
with codecs.open(filepath+outputfile,'w+', encoding='utf8') as fo:
    for fi in  file_list:
        with open(json_path+fi) as f:
            data=json.load(f)
        #seperate the sentencess
        sen_text=data['text']
        #did not remove the title
        sen_text_title_index=sen_text.find('.\n')
        sen_text=sen_text.replace('\n',' ')
        #sen_text=sen_text_title_split[1]
        sen_text_index=[i for i in range(len(sen_text)) if sen_text.startswith('. ',i)]
        sen_text_index=[i+2 for i in sen_text_index]
        sen_text_index=[0]+sen_text_index+[len(sen_text)]

        #build for a dictionary for the sentence index
        sen_index_dic={}
        for si in range(len(sen_text_index)-1):
            sen_index_dic[si+1]=[sen_text_index[si],sen_text_index[si+1]]
        '''
        for si in range(len(sen_text_index)-1):
            print(si)
            print(sen_text[sen_text_index[si]:sen_text_index[si+1]]+'\"')
        '''
        #build the relation set for the searching purpose
        all_relation=[]
        if 'relations' in data.keys():
            for ri in data['relations']:
                e1=ri['subj']
                e2=ri['obj']
                relation_pair=[e1,e2]
                relation_pair.sort()
                all_relation.append(relation_pair)

        #build the protein and location list
        protein_list=[]
        location_list=[]
        if 'denotations' in data.keys():
            for di in data['denotations']:
                #get all the entities we waned

                if di['obj'].startswith('uniprot'):
                    protein_list.append(di['id'])
                if di['obj'].startswith('go'):
                    location_list.append(di['id'])

        #consider all the possible combinations
        for pi in protein_list:
            for li in location_list:
                cur_pair=[pi,li]

                #get the entity loaction
                for di in data['denotations']:
                    if di['id']==cur_pair[0]:
                        e1_start=di['span']['begin']
                        e1_end=di['span']['end']
                    if di['id']==cur_pair[1]:
                        e2_start=di['span']['begin']
                        e2_end=di['span']['end']
                two_entity_in_one_sentence=False
                sen_index=-1
                for ei in sen_index_dic.keys():
                    if e1_start>=sen_index_dic[ei][0] and e1_end<sen_index_dic[ei][1]\
                        and e2_start>=sen_index_dic[ei][0] and e2_end<sen_index_dic[ei][1]:
                        two_entity_in_one_sentence=True
                        sen_index=ei
                #print(sen_text[e1_start:e1_end],'-',sen_text[e2_start:e2_end])
                if two_entity_in_one_sentence:
                    this_sent=sen_text[sen_index_dic[sen_index][0]:sen_index_dic[sen_index][1]]
                    e1_start_new=e1_start-sen_index_dic[sen_index][0]
                    e1_end_new=e1_end-sen_index_dic[sen_index][0]
                    e2_start_new=e2_start-sen_index_dic[sen_index][0]
                    e2_end_new=e2_end-sen_index_dic[sen_index][0]
                    if e2_end>e1_end:
                        this_sent=this_sent[0:e2_end_new]+'</Subcellular_location>'+this_sent[e2_end_new:]
                        this_sent=this_sent[0:e2_start_new]+'<Subcellular_location>'+this_sent[e2_start_new:]
                        this_sent=this_sent[0:e1_end_new]+'</Gene>'+this_sent[e1_end_new:]
                        this_sent=this_sent[0:e1_start_new]+'<Gene>'+this_sent[e1_start_new:]
                    else:
                        this_sent=this_sent[0:e1_end_new]+'</Gene>'+this_sent[e1_end_new:]
                        this_sent=this_sent[0:e1_start_new]+'<Gene>'+this_sent[e1_start_new:]
                        this_sent=this_sent[0:e2_end_new]+'</Subcellular_location>'+this_sent[e2_end_new:]
                        this_sent=this_sent[0:e2_start_new]+'<Subcellular_location>'+this_sent[e2_start_new:]

                    cur_pair.sort()
                    if cur_pair in all_relation:
                        #this is a positive instance

                        fo.write('R_PLOC ')
                        fo.write(data['sourceid'])
                        fo.write(' ')
                        fo.write(str(sen_index))
                        fo.write('\t')
                        fo.write(this_sent)
                        fo.write('\n')
                    else:
                        #this is a negative instance
                        fo.write('NONE ')
                        fo.write(data['sourceid'])
                        fo.write(' ')
                        fo.write(str(sen_index))
                        fo.write('\t')
                        fo.write(this_sent)
                        fo.write('\n')
