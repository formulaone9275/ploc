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

                text=instance_text
                protein_name=''
                location_name=''
                with codecs.open(inputfolder+fi, encoding='utf8') as f:
                    for line in f:
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

                                w_ind+=1
                                w_ind_loc+=1
                            location_name=location_name.strip()
                            protein_name=protein_name.strip()
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

                                w_ind+=1

                            protein_name=protein_name.strip()
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

                                w_ind_loc+=1
                            location_name=location_name.strip()
                if location_name not in loc_dic.keys():
                    loc_dic[location_name]=[protein_name]
                else:
                    loc_dic[location_name].append(protein_name)

                text_low=text.lower()
                ind=0
                while True:

                    print(fi)
                    print(protein_name)
                    print(text)
                    gene_index=text_low.find(protein_name,ind)
                    if gene_index==-1:
                        protein_name_sp=protein_name.split(' ')
                        protein_name=protein_name_sp[0]
                    if text_low[gene_index-1] not in [' ','-','/',',','.','>'] or text_low[gene_index+len(protein_name)] not in [' ','-','/',',','.','<']:
                        ind=gene_index+len(protein_name)
                    else:
                        break
                text=text[0:gene_index+len(protein_name)]+'</Gene>'+text[gene_index+len(protein_name):]
                text=text[0:gene_index]+'<Gene>'+text[gene_index:]

                text_low=text.lower()
                indd=0
                while True:

                    print(fi)
                    print(location_name)
                    print(text)
                    loc_index=text_low.find(location_name,indd)

                    if text_low[loc_index-1] not in [' ','-','/',',','.','(','>'] or text_low[loc_index+len(location_name)] not in [' ','-','/',',','.',';','<']:

                        indd=loc_index+len(location_name)
                    else:
                        break
                text=text[0:loc_index+len(location_name)]+'</Subcellular_location>'+text[loc_index+len(location_name):]
                text=text[0:loc_index]+'<Subcellular_location>'+text[loc_index:]
                fo.write('R_PLOC ')
                fo.write(fi)
                fo.write('\t')
                fo.write(text)
                fo.write('\n')

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
    outputfile='./data/pos.txt'
    negfolder='./data/YPD/all/'
    negfile='./data/neg.txt'
    picklefile='./data/location_dic.pickle'
    #generate_pos_ploc_from_craven(inputfolder,outputfile,picklefile)
    generate_neg_ploc_from_craven(negfolder,negfile,picklefile)
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