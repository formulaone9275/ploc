import codecs
import os
import pickle


def count_plot_loc():
    filename='./data/test/ploc_original.txt'
    loc_dic={}
    fre_dic={}
    with codecs.open(filename, encoding='utf8') as f:
        for line in f:
            if not line.startswith(('NONE', 'R_PLOC')):
                continue

            line = line.strip()
            info, tagged = line.split('\t')
            gene_start=tagged.find('<Gene>')
            gene_end=tagged.find('</Gene>')
            location_start=tagged.find('<Subcellular_location>')
            location_end=tagged.find('</Subcellular_location>')
            loc_name=tagged[location_start+22:location_end]
            gene_name=tagged[gene_start+6:gene_end]
            if loc_name not in loc_dic.keys():
                loc_dic[loc_name]=[gene_name]
            else:
                loc_dic[loc_name].append(gene_name)
    for li in loc_dic.keys():
        fre_dic[li]=len(loc_dic[li])
    print(loc_dic.keys())
    print(len(fre_dic))
    s_fre_dic=sorted(fre_dic.items(),key=lambda x:x[1],reverse=True)
    print(s_fre_dic)
    return loc_dic

def generate_pos_ploc_from_craven(inputfolder,outputfile,picklefile):

    file_list=os.listdir(inputfolder)
    loc_dic={}

    with codecs.open(outputfile,'w+', encoding='utf8') as fo:
        for fi in file_list:
            line_count=1
            pair_list=[]
            #deal with the case that location name and protein name appear multiple times
            loc_index_dic={}
            protein_index_dic={}
            with codecs.open(inputfolder+fi, encoding='utf8') as f:
                for line in f:
                    if line_count==1:
                        line=line.replace('\"','')
                        instance_text=line.strip()

                    elif line_count==3:

                        pair_index1=[i for i in range(len(line)) if line.startswith('[',i)]
                        pair_index2=[i for i in range(len(line)) if line.startswith(']',i)]
                        pair_index1.sort()
                        pair_index2.sort()

                        for pi in range(len(pair_index1)):
                            pair_str=line[pair_index1[pi]+1:pair_index2[pi]]
                            pair_num_str=pair_str.split(',')
                            pair_list.append(pair_num_str)

                    line_count+=1
            #generate instance for each protein-location pair
            first_pair=True
            first_pair_loc=True
            for ii in pair_list:

                with codecs.open(inputfolder+fi, encoding='utf8') as f:
                    for line in f:
                        protein_name=''
                        location_name=''
                        line=line.strip()
                        if line.startswith(ii[0]+' ') and line.startswith(ii[1]+' '):
                            line_seg=line.split(' ')
                            w_ind=1
                            pre_ind=0
                            w_ind_loc=1
                            pre_ind_loc=0
                            for wi in line_seg:
                                if ':PROTEIN}' in wi or ':PROTEIN:LOCATION}' in wi:
                                    symbol_index=wi.find('{')
                                    if pre_ind==0:
                                        pre_ind=w_ind
                                        protein_name_part=wi[0:symbol_index]
                                        protein_name_part=protein_name_part.replace('_',' ')
                                        protein_name+=(protein_name_part+' ')
                                    elif (w_ind-pre_ind)==1:
                                        pre_ind=w_ind
                                        protein_name_part=wi[0:symbol_index]
                                        protein_name_part=protein_name_part.replace('_',' ')
                                        protein_name+=(protein_name_part+' ')
                                    else:
                                        if first_pair is False:
                                            protein_name_part=wi[0:symbol_index]
                                            protein_name_part=protein_name_part.replace('_',' ')
                                            protein_name=protein_name_part
                                            pre_ind=w_ind
                                        first_pair=False
                                else:
                                    if protein_name!='':
                                        protein_name=protein_name.strip()
                                        if protein_name not in protein_index_dic.keys():
                                            protein_index_dic[protein_name]=[ii[0]]
                                        else:
                                            protein_index_dic[protein_name].append(ii[0])
                                if ':LOCATION}' in wi:
                                    symbol_index=wi.find('{')
                                    if pre_ind_loc==0:
                                        pre_ind_loc=w_ind_loc
                                        location_name_part=wi[0:symbol_index]
                                        location_name_part=location_name_part.replace('_',' ')
                                        location_name+=(location_name_part+' ')
                                    elif (w_ind_loc-pre_ind_loc)==1:
                                        pre_ind_loc=w_ind_loc
                                        location_name_part=wi[0:symbol_index]
                                        location_name_part=location_name_part.replace('_',' ')
                                        location_name+=(location_name_part+' ')
                                    else:
                                        if first_pair_loc is False:
                                            location_name_part=wi[0:symbol_index]
                                            location_name_part=location_name_part.replace('_',' ')
                                            location_name=location_name_part
                                            pre_ind_loc=w_ind_loc
                                        first_pair_loc=False
                                else:
                                    if location_name!='':
                                        location_name=location_name.strip()
                                        if location_name not in loc_index_dic.keys():
                                            loc_index_dic[location_name]=[ii[1]]
                                        else:
                                            loc_index_dic[location_name].append(ii[1])
                                w_ind+=1
                                w_ind_loc+=1
                            protein_name=protein_name.strip()
                            if protein_name not in protein_index_dic.keys():
                                protein_index_dic[protein_name]=[ii[0]]
                            else:
                                protein_index_dic[protein_name].append(ii[0])
                            location_name=location_name.strip()
                            if location_name not in loc_index_dic.keys():
                                loc_index_dic[location_name]=[ii[1]]
                            else:
                                loc_index_dic[location_name].append(ii[1])
                        elif line.startswith(ii[0]+' '):
                            line_seg=line.split(' ')
                            w_ind=1
                            pre_ind=0
                            for wi in line_seg:

                                if ':PROTEIN}' in wi or ':PROTEIN:LOCATION}' in wi:
                                    symbol_index=wi.find('{')
                                    if pre_ind==0:
                                        pre_ind=w_ind
                                        protein_name_part=wi[0:symbol_index]
                                        protein_name_part=protein_name_part.replace('_',' ')
                                        protein_name+=(protein_name_part+' ')
                                    elif (w_ind-pre_ind)==1:
                                        pre_ind=w_ind
                                        protein_name_part=wi[0:symbol_index]
                                        protein_name_part=protein_name_part.replace('_',' ')
                                        protein_name+=(protein_name_part+' ')
                                    else:
                                        if first_pair is False:
                                            protein_name_part=wi[0:symbol_index]
                                            protein_name_part=protein_name_part.replace('_',' ')
                                            protein_name=protein_name_part
                                            pre_ind=w_ind
                                        first_pair=False
                                else:
                                    if protein_name!='':
                                        protein_name=protein_name.strip()
                                        if protein_name not in protein_index_dic.keys():
                                            protein_index_dic[protein_name]=[ii[0]]
                                        else:
                                            protein_index_dic[protein_name].append(ii[0])
                                w_ind+=1

                            protein_name=protein_name.strip()
                            if protein_name not in protein_index_dic.keys():
                                protein_index_dic[protein_name]=[ii[0]]
                            else:
                                protein_index_dic[protein_name].append(ii[0])
                        elif line.startswith(ii[1]+' '):
                            line_seg=line.split(' ')
                            w_ind_loc=1
                            pre_ind_loc=0
                            for wi in line_seg:
                                if ':LOCATION}' in wi:
                                    symbol_index=wi.find('{')

                                    if pre_ind_loc==0:
                                        pre_ind_loc=w_ind_loc
                                        location_name_part=wi[0:symbol_index]
                                        location_name_part=location_name_part.replace('_',' ')
                                        location_name+=(location_name_part+' ')
                                    elif (w_ind_loc-pre_ind_loc)==1:
                                        pre_ind_loc=w_ind_loc
                                        location_name_part=wi[0:symbol_index]
                                        location_name_part=location_name_part.replace('_',' ')
                                        location_name+=(location_name_part+' ')
                                    else:
                                        if first_pair_loc is False:
                                            location_name_part=wi[0:symbol_index]
                                            location_name_part=location_name_part.replace('_',' ')
                                            location_name=location_name_part
                                            pre_ind_loc=w_ind_loc
                                        first_pair_loc=False
                                else:
                                    if location_name!='':
                                        location_name=location_name.strip()
                                        if location_name not in loc_index_dic.keys():
                                            loc_index_dic[location_name]=[ii[1]]
                                        else:
                                            loc_index_dic[location_name].append(ii[1])
                                w_ind_loc+=1
                            location_name=location_name.strip()
                            if location_name not in loc_index_dic.keys():
                                loc_index_dic[location_name]=[ii[1]]
                            else:
                                loc_index_dic[location_name].append(ii[1])

            #add something here

            for kli in loc_index_dic.keys():
                loc_index_dic[kli]=list(set(loc_index_dic[kli]))
            for kpi in protein_index_dic.keys():
                protein_index_dic[kpi]=list(set(protein_index_dic[kpi]))
            #print(loc_index_dic)
            #print(protein_index_dic)
            #add the order index for the proteins and location names if they appear more than one time
            loc_appear_index_dic={}
            protein_appear_index_dic={}
            for kli in loc_index_dic.keys():
                if len(loc_index_dic[kli])==1:
                    loc_appear_index_dic[kli]=[1]
                else:
                    #generate the number list
                    num_list=[int(ni) for ni in loc_index_dic[kli]]
                    num_list.sort()
                    temp_index_list=[num_list.index(int(ni))+1 for ni in loc_index_dic[kli]]
                    loc_appear_index_dic[kli]=temp_index_list

            for kpi in protein_index_dic.keys():
                if len(protein_index_dic[kpi])==1:
                    protein_appear_index_dic[kpi]=[1]
                else:
                    #generate the number list
                    num_list=[int(ni) for ni in protein_index_dic[kpi]]
                    num_list.sort()
                    temp_index_list=[num_list.index(int(ni))+1 for ni in protein_index_dic[kpi]]
                    protein_appear_index_dic[kpi]=temp_index_list

            #print(loc_appear_index_dic)
            #print(protein_appear_index_dic)
            first_appear_in_one_line=True
            #first_location_in_one_line=True
            for idi in range(len(pair_list)):
                ii=pair_list[idi]
                text=instance_text

                if idi !=pair_list.index(ii):
                    first_appear_in_one_line=False
                #get the protein and location name
                #in some cases, there are more than one protein/location name in one line
                for kpi in protein_index_dic.keys():
                    if ii[0] in protein_index_dic[kpi]:
                        if first_appear_in_one_line is True:
                            protein_name=kpi
                            protein_appear_index=protein_appear_index_dic[kpi][protein_index_dic[kpi].index(ii[0])]
                            break
                        else:
                            first_appear_in_one_line=True
                            continue
                for kli in loc_index_dic.keys():
                    if ii[1] in loc_index_dic[kli]:
                        location_name=kli
                        location_appear_index=loc_appear_index_dic[kli][loc_index_dic[kli].index(ii[1])]
                        break

                #print('protein appear index:',protein_appear_index)
                text_low=text.lower()
                ind=0
                cur_protein_appear_index=1
                cur_loc_appear_index=1
                break_count=0
                while True:

                    #print(fi)
                    #print(protein_name)
                    #print(text)
                    gene_index=text_low.find(protein_name,ind)
                    if gene_index==-1:
                        protein_name_sp=protein_name.split(' ')
                        protein_name=protein_name_sp[0]
                    if text_low[gene_index-1] not in [' ','-','/',',','.','>'] or text_low[gene_index+len(protein_name)] not in [' ','-','/',',','.','<']:
                        ind=gene_index+len(protein_name)
                    else:
                        #print('Current protein index:',cur_protein_appear_index)
                        if protein_appear_index>cur_protein_appear_index:
                            ind=gene_index+len(protein_name)
                            cur_protein_appear_index+=1
                        else:
                            break
                    if break_count>5:
                        break
                    break_count+=1
                if gene_index==-1:
                    print(fi)
                    print(protein_name)
                    print(text)
                text=text[0:gene_index+len(protein_name)]+'</Gene>'+text[gene_index+len(protein_name):]
                text=text[0:gene_index]+'<Gene>'+text[gene_index:]

                text_low=text.lower()
                indd=0
                break_count=0
                while True:

                    #print(fi)
                    #print(location_name)
                    #print(text)
                    loc_index=text_low.find(location_name,indd)

                    if text_low[loc_index-1] not in [' ','-','/',',','.','(','>'] or text_low[loc_index+len(location_name)] not in [' ','-','/',',','.',';','<']:

                        indd=loc_index+len(location_name)
                    else:
                        if location_appear_index>cur_loc_appear_index:
                            indd=loc_index+len(location_name)
                            cur_loc_appear_index+=1
                        else:
                            break
                    if break_count>5:
                        break
                    break_count+=1
                text=text[0:loc_index+len(location_name)]+'</Subcellular_location>'+text[loc_index+len(location_name):]
                text=text[0:loc_index]+'<Subcellular_location>'+text[loc_index:]
                fo.write('R_PLOC ')
                fo.write(fi)
                fo.write('\t')
                fo.write(text)
                fo.write('\n')
                fo.write('\n')
                fo.write('\n')
                fo.write('\n')
                #print(location_name,'-',protein_name)
                if location_name not in loc_dic.keys():
                    loc_dic[location_name]=[protein_name]
                else:
                    loc_dic[location_name].append(protein_name)
    print(file_list)
    fre_dic={}
    print(loc_dic)
    print(len(loc_dic))
    for li in loc_dic.keys():
        fre_dic[li]=len(loc_dic[li])

    s_fre_dic=sorted(fre_dic.items(),key=lambda x:x[1],reverse=True)
    with open(picklefile, 'wb') as fp:
        # Pickle the 'data' dictionary using the highest protocol available.
        pickle.dump(s_fre_dic, fp, pickle.HIGHEST_PROTOCOL)

#generate the negative instance of ploc
def generate_neg_ploc_from_craven(inputfolder,outputfile,picklefile):

    file_list=os.listdir(inputfolder)
    with open(picklefile, 'rb') as f:
        loc_name_dic = pickle.load(f,encoding="latin1")

    loc_name_list=list(loc_name_dic.keys())
    protein_name_list_=[]
    for ki in loc_name_dic.keys():
        protein_name_list_+=loc_name_dic[ki]
    protein_name_list_=list(set(protein_name_list_))
    #conver the protein names to lower case
    protein_name_list=[]
    for pi in protein_name_list_:
        pi=pi.lower()
        if pi not in protein_name_list:
            protein_name_list.append(pi)


    with codecs.open(outputfile,'w+', encoding='utf8') as fo:
        for fi in file_list:
            line_count=1
            pos_ind=False
            with codecs.open(inputfolder+fi, encoding='utf8') as f:
                for line in f:
                    if line_count==1:
                        line=line.replace('\"','')
                        instance_text=line.strip()

                    elif line_count==2:
                        if line.startswith('['):
                            pos_ind=True
                            break

                    line_count+=1
            if pos_ind is True:
                continue
            text=instance_text
            text_low=text.lower()
            found_protein=False
            for protein_name in protein_name_list:
                protein_name=protein_name.lower()
                ind=0
                protein_search_count=0
                while True:

                    protein_search_count+=1
                    gene_index=text_low.find(protein_name,ind)
                    if gene_index==-1:
                        protein_name_sp=protein_name.split(' ')
                        protein_name=protein_name_sp[0]
                    if text_low[gene_index-1] not in [' ','-','/',',','.','>'] or text_low[gene_index+len(protein_name)] not in [' ','-','/',',','.','<']:
                        ind=gene_index+len(protein_name)
                    else:
                        break
                    if protein_search_count>10:
                        break
                if gene_index>-1 and protein_search_count<=10 and found_protein is False:
                    text=text[0:gene_index+len(protein_name)]+'</Gene>'+text[gene_index+len(protein_name):]
                    text=text[0:gene_index]+'<Gene>'+text[gene_index:]
                    found_protein=True

            if found_protein is False:
                print(fi)
                print(text)
            text_low=text.lower()
            found_loc=False
            for location_name in loc_name_list:
                location_name=location_name.lower()
                loc_search_count=0
                indd=0
                while True:

                    loc_search_count+=1
                    loc_index=text_low.find(location_name,indd)

                    if text_low[loc_index-1] not in [' ','-','/',',','.','(','>'] or text_low[loc_index+len(location_name)] not in [' ','-','/',',','.',';','<']:

                        indd=loc_index+len(location_name)
                    else:
                        break
                    if loc_search_count>10:
                        break
                if loc_index>-1 and loc_search_count<=10 and found_loc is False:
                    text=text[0:loc_index+len(location_name)]+'</Subcellular_location>'+text[loc_index+len(location_name):]
                    text=text[0:loc_index]+'<Subcellular_location>'+text[loc_index:]
                    found_loc=True

            if found_loc is False:
                print(fi)
                print(text)
            if found_protein is True and found_loc is True:
                fo.write('NONE ')
                fo.write(fi)
                fo.write('\t')
                fo.write(text)
                fo.write('\n')



if __name__ == '__main__':

    inputfolder='./data/YPD/pos/'
    outputfile='./data/pos_new.txt'
    negfolder='./data/YPD/all/'
    negfile='./data/neg.txt'
    picklefile='./data/location_dic.pickle'
    generate_pos_ploc_from_craven(inputfolder,outputfile,picklefile)
    #generate_neg_ploc_from_craven(negfolder,negfile,picklefile)
    #test of reading the pickle file
    with open(picklefile, 'rb') as f:
        data = pickle.load(f,encoding="latin1")
    print(data)
    '''
    #count top data
    topnum=10
    top_count=0
    for ti in range(topnum):
        top_count+=data[ti][1]
    print(top_count)

    '''