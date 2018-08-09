import pickle
import codecs


def filter_ds_data(filepath,filename,topnum):
    picklefile='./data/location_frequency.pickle'

    #test of reading the pickle file
    with open(picklefile, 'rb') as f:
        data = pickle.load(f,encoding="latin1")
    loc_name_list=[]
    for ti in range(topnum):
        loc_name_list.append(data[ti][0])
    #generate the file based on the input filename
    filename_split=filename.split('.')
    print(filename_split)
    outputfile=filename_split[0]+'_top'+str(topnum)+'.txt'
    with codecs.open(filepath+outputfile,'w+', encoding='utf8') as fo:
        with codecs.open(filepath+filename, encoding='utf8') as f:
            for line in f:
                if not line.startswith(('NONE', 'R_PLOC')):
                    continue

                line = line.strip()
                info, tagged = line.split('\t')
                loc_index_start=tagged.find('<9>')
                loc_index_end=tagged.find('</9>')
                loc_name=tagged[loc_index_start+3:loc_index_end]
                loc_name=loc_name.lower()
                if loc_name in loc_name_list:

                    fo.write(line)
                    fo.write('\n')

def filter_test_data(filepath,filename,topnum):
    picklefile='./data/location_frequency.pickle'

    #test of reading the pickle file
    with open(picklefile, 'rb') as f:
        data = pickle.load(f,encoding="latin1")
    loc_name_list=[]
    for ti in range(topnum):
        loc_name_list.append(data[ti][0])

    positive_count=0
    negative_count=0
    outputfile='ploc_top.txt'
    with codecs.open(filepath+outputfile,'w+', encoding='utf8') as fo:
        with codecs.open(filepath+filename, encoding='utf8') as f:
            for line in f:
                if not line.startswith(('NONE', 'R_PLOC')):
                    continue

                line = line.strip()
                info, tagged = line.split('\t')
                loc_index_start=tagged.find('<Subcellular_location>')
                loc_index_end=tagged.find('</Subcellular_location>')
                loc_name=tagged[loc_index_start+len('<Subcellular_location>'):loc_index_end]
                loc_name=loc_name.lower()
                if loc_name in loc_name_list:
                    if line.startswith('NONE'):
                        negative_count+=1
                    elif line.startswith('R_PLOC'):
                        positive_count+=1

                    fo.write(line)
                    fo.write('\n')

if __name__ == '__main__':
    filepath='./data/test/'
    '''
    for filename in ['ploc.txt','CP.txt','CP_TW.txt','CP_HP.txt','filtered.txt']:

        topnum=10
        filter_ds_data(filepath,filename,topnum)
    '''
    topnum=10
    filename='ploc_original.txt'
    filter_test_data(filepath,filename,topnum)