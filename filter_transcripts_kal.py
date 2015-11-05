import sys


def do_all(reconstr_per,fp_file,kal_file):
    # def write_cutoff(cutoff, weights, transcripts, false_positives):
    #     weights = weights[int(0.1*cutoff*len(weights)):]
    #     good_names = set([x[0] for x in weights])
    #     transcripts = [(name, seq) for name, seq in transcripts if name in good_names]
    #     hits = 0
    #     for name, seq in transcripts:
    #         if name in false_positives:
    #             hits += 1
    #     print('{} false positives for {}'.format(hits, cutoff))
    #     util.seqs_to_fasta('drop{}/transcripts.fa'.format(cutoff), transcripts)
    len_threshold = 200
    cutoff_range = range(10)



    ab_dict  = {}
    weights = [] 
    with open(kal_file) as f:
        f.readline()
        for line in f:
            name, tr_len, _, _, weight = line.split()
            if float(tr_len) < 200:
                continue
            weights.append((name, float(weight)))
            ab_dict[name] = [float(weight)]
    weights.sort(key = lambda x: x[1])
    
    
    weight_cutoff = {}
    tr_rec_dict = {}

    for cutoff in cutoff_range:
        weight_cutoff[cutoff] = weights[int(0.1*cutoff*len(weights))]
        tr_rec_dict[cutoff] = {}
        for tr in ab_dict:
            tr_rec_dict[cutoff][tr] = 0


    f=open(reconstr_per,'r')
    #lines = f.readlines()
    best_rec_dict = {}
    not_first_time = 0
    for line in f:
        '''if i<5:
            continue'''
        tokens = line.split()
        org = tokens[9]; #need to select for the tokens[9][29:]  when having larger prefix
        rec = tokens[13]; rec_len = int(tokens[0]); tr_len = int(tokens[10])
        if rec_len > 0.9*tr_len: 
            for cutoff in cutoff_range:
                if ab_dict.get(rec,0) > weight_cutoff[cutoff]:
                    #z = best_rec_dict[cutoff][0].get(org,[None,0]);                
                    #if rec_len>z[1]:
                    #    best_rec_dict[cutoff][org] = [rec,rec_len,tr_len,tr_ab/tr_len];
                    tr_rec_dict[cutoff][org] = 1

    recall = {}
    for cutoff in cutoff_range:
        recall[cutoff] = sum(tr_rec_dict[cutoff].values())

    
    tot_rec = {}
    fp_rec = {}
    for cutoff in cutoff_range:
        tot_rec[cutoff] = 0
        fp_rec[cutoff] = 0


    with open(fp_file) as f:
        for line in f:
            name, code, length = line.split()
            if int(length) < len_threshold: continue
            for cutoff in cutoff_range:
                if ab_dict.get(name,0) > weight_cutoff[cutoff]:
                    tot_rec[cutoff] +=1
                    if code == '0':
                        fp_rec[cutoff] =1

def main():
    if len(sys.argv) == 1:
        arguments = ['', 'reconstr_per.txt', 'reconstructed_fp.log', 'abundance.kal']
    else:
        arguments = sys.argv
    reconstr_per,fp_file,kal_file = arguments[1], arguments[2],arguments[3] 
    do_all(reconstr_per,fp_file,kal_file)

if __name__ == '__main__':
    main()
