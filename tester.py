import operator, pdb, os
#from Bio.Seq import Seq
#from Bio.Alphabet import generic_dna

def reverse_complement(s):
    '''Find the reverse complement of a DNA string using biopython interface'''
    replacements = [('A', 't'), ('T', 'a'), ('C', 'g'), ('G', 'c')]
    for ch1, ch2 in replacements:
        bases = re.sub(ch1, ch2, s.upper())
        return bases[::-1].upper()
    '''d = Seq(s,generic_dna)
    return str(d.reverse_complement())'''

def check_palindromes(directory,filename):
    f = open(directory+'temp.dict','w')
    for line in open(directory+filename):
        fields=line.strip().split();
        if (reverse_complement(fields[0]) == fields[0]):
            f.write(line)


def double_strandify(infile,outfile):
    ''' Input is a directory and filename. Create a temporary file named temp.dict and swap the output to the input file'''
    f = open(outfile,'w')
    for line in open(infile):
        fields=line.strip().split();
        
        if (reverse_complement(fields[0]) != fields[0]):
            f.write(line)
            f.write(reverse_complement(fields[0])+'\t'+fields[1]+'\n')
        else:
            f.write(fields[0]+'\t '+str(2*int(fields[1]))+'\n')
    #f.close()

def analyzer_blat(Repeats_File,Dest_File,exp_file,N):
    tr_lengths = {}
    tr_abundance = {}
    sum_exp = 0
    sum_len = 0
    for lines in open(exp_file):
        len_field = -1  #original: len_field = 1
        ab_field = 1 #original: ab_field = 7
        fields=lines.strip().split();
        if lines[0]=='#' or fields[len_field]==0:
            continue;
        if len_field != -1:
            tr_lengths[fields[0]] = int(fields[len_field])
            sum_len += int(fields[len_field])
        else:
            tr_lengths[fields[0]] = 1
            sum_len += 1
        tr_abundance[fields[0]] = float(fields[ab_field])
        sum_exp += float(fields[ab_field])
        

  
    f=open(Repeats_File,'r')
    lines = f.readlines()
    best_rec_dict = {}
    not_first_time = 0
    for (i,line) in enumerate(lines):
        if i<5:
            continue
        #pdb.set_trace() 
        #print('lenth',len(lines))
        tokens = line.split()
        org = tokens[9]; #need to select for the tokens[9][29:]  when having larger prefix
        rec = tokens[13]; rec_len = int(tokens[0]); tr_len = int(tokens[10])
        tr_ab = tr_abundance.get(org,0)
        z = best_rec_dict.get(org,[None,0]);                
        if rec_len>z[1]:
            best_rec_dict[org] = [rec,rec_len,tr_len,tr_ab/tr_len];
    perf = 0
    sum_ret = 0
    with open(Dest_File,'w') as write_file:
        for (org,z) in best_rec_dict.iteritems():
            write_file.write(org+'\t'+str(z[0])+'\t'+str(z[1])+'\t'+str(z[2]) + '\t' + str(z[3])+'\n')
            perf += float(z[1])/float(z[2])
            sum_ret += z[1]
        if len(best_rec_dict):    
            write_file.write("#Average fractional contig of transcripts retrived:" + '\t' + str(perf/len(best_rec_dict)) + '\n')
        else:
            write_file.write("#Average fractional contig of transcripts retrived:" + '\t' + str(perf) + '  best_rec=0' + '\n')
        if sum_len > 0:
            write_file.write("#Sum of positions retreived / sum of lengths:" + '\t' + str(float(sum_ret)/float(sum_len)) + '\n')

def analyzer_blat_rsem(Repeats_File,Dest_File,exp_file,N):
    tr_lengths = {}
    tr_abundance = {}
    sum_exp = 0
    sum_len = 0
    for lines in open(exp_file):
        len_field = 2  #original: len_field = 1
        ab_field = 5 #original: ab_field = 7
        fields=lines.strip().split();
        if lines[0]=='#' or fields[len_field]==0:
            continue;
        tr_lengths[fields[0]] = int(fields[len_field])
        tr_abundance[fields[0]] = float(fields[ab_field])
        sum_exp += float(fields[ab_field])
        sum_len += int(fields[len_field])

  
    f=open(Repeats_File,'r')
    lines = f.readlines()
    best_rec_dict = {}
    not_first_time = 0
    for (i,line) in enumerate(lines):
        if i<5:
            continue
        #pdb.set_trace() 
        #print('lenth',len(lines))
        tokens = line.split()
        org = tokens[9]; #need to select for the tokens[9][29:]  when having larger prefix
        rec = tokens[13]; rec_len = int(tokens[0]); tr_len = int(tokens[10])
        tr_ab = tr_abundance.get(org,0)
        z = best_rec_dict.get(org,[None,0]);                
        if rec_len>z[1]:
            best_rec_dict[org] = [rec,rec_len,tr_len,tr_ab/tr_len];
    perf = 0
    sum_ret = 0
    with open(Dest_File,'w') as write_file:
        for (org,z) in best_rec_dict.iteritems():
            write_file.write(org+'\t'+str(z[0])+'\t'+str(z[1])+'\t'+str(z[2]) + '\t' + str(z[3])+'\n')
            perf += float(z[1])/float(z[2])
            sum_ret += z[1]
        if len(best_rec_dict):    
            write_file.write("#Average fractional contig of transcripts retrived:" + '\t' + str(perf/len(best_rec_dict)) + '\n')
        else:
            write_file.write("#Average fractional contig of transcripts retrived:" + '\t' + str(perf) + '  best_rec=0' + '\n')
        if sum_len > 0:
            write_file.write("#Sum of positions retreived / sum of lengths:" + '\t' + str(float(sum_ret)/float(sum_len)) + '\n')


def analyzer_blat_noExp(Repeats_File,Dest_File,exp_file,N):
    f=open(Repeats_File,'r')
    lines = f.readlines()
    best_rec_dict = {}
    not_first_time = 0
    for (i,line) in enumerate(lines):
        '''if i<5:
            continue'''
        #pdb.set_trace() 
        #print('lenth',len(lines))
        tokens = line.split()
        org = tokens[9]; #need to select for the tokens[9][29:]  when having larger prefix
        rec = tokens[13]; rec_len = int(tokens[0]); tr_len = int(tokens[10])
        rec_tr_len = int(tokens[14]); #tr_abundance.get(org,0)
        z = best_rec_dict.get(org,[None,0]);                
        if rec_len>z[1]:
            best_rec_dict[org] = [rec,rec_len,tr_len,rec_tr_len];
    perf = 0
    sum_ret = 0
    no_tr_rec = 0 
    with open(Dest_File,'w') as write_file:
        for (org,z) in best_rec_dict.iteritems():
            write_file.write(org+'\t'+str(z[0])+'\t'+str(z[1])+'\t'+str(z[2]) + '\t' + str(z[3])+'\n')
            if z[1] >= 0.9*z[2]:
                no_tr_rec +=1
            perf += float(z[1])/float(z[2])
            sum_ret += z[1]
        if len(best_rec_dict):    
            write_file.write("#Average fractional contig of transcripts retrived:" + '\t' + str(perf/len(best_rec_dict)) + '\n')
        else:
            write_file.write("#Average fractional contig of transcripts retrived:" + '\t' + str(perf) + '  best_rec=0' + '\n')
        write_file.write("# of transcripts at greater than 90%:" + '\t' + str(no_tr_rec))


def analyzer(Repeats_File,Dest_File,exp_file,N):
    tr_lengths = {}
    tr_abundance = {}
    sum_exp = 0
    sum_len = 0
    for lines in open(exp_file):
        fields=lines.strip().split();
        if lines[0]=='#' or fields[7]==0:
            continue;
        tr_lengths[fields[0]] = int(fields[1])
        tr_abundance[fields[0]] = float(fields[7])
        sum_exp += float(fields[7])
        sum_len += int(fields[1])

  
    f=open(Repeats_File,'r')
    lines = f.readlines()
    best_rec = []
    not_first_time = 0
    for (i,line) in enumerate(lines):
        #pdb.set_trace() 
        #print('lenth',len(lines))
        tokens = line.split()
        if tokens[0] == ">" or i == len(lines)-1:
            if not_first_time:
                sorted_list = sorted(tempdict.iteritems(), key=operator.itemgetter(1),reverse=True)
                tr_len = tr_lengths.get(identity,-1)
                tr_ab = tr_abundance.get(identity,0)
                #print(sorted_list)
                #raw_input()
                if len(sorted_list) > 0:
                    best_rec.append([identity, sorted_list[0][0],int(sorted_list[0][1]),tr_len,tr_ab/tr_len])
                else:
                    best_rec.append([identity, 'None',0, tr_len,tr_ab/tr_len])
            not_first_time = 1 
            identity = tokens[1][29:]
            tempdict  = {}
        else:
            tempdict[tokens[0]] = tempdict.get(tokens[0],0)+int(tokens[3])
    perf = 0
    sum_ret = 0
    with open(Dest_File,'w') as write_file:
        for rec in best_rec:
            write_file.write(str(rec[0])+'\t'+str(rec[1])+'\t'+str(rec[2])+'\t'+str(rec[3]) + '\t' + str(rec[4])+'\n')
            perf += float(rec[2])/float(rec[3])
            sum_ret += rec[2]
        if len(best_rec):    
            write_file.write("#Average fractional contig of transcripts retrived:" + '\t' + str(perf/len(best_rec)) + '\n')
        else:
            write_file.write("#Average fractional contig of transcripts retrived:" + '\t' + str(perf) + '  best_rec=0' + '\n')
        if sum_len > 0:
            write_file.write("#Sum of positions retreived / sum of lengths:" + '\t' + str(float(sum_ret)/float(sum_len)) + '\n')


def ab_performance(log_file,exp_file):
    '''Log is the performance log of a given reconstruction with respect to the true transcripts. Exp_file is a three column file with the third column comprising values proportional to abundance'''
    tr_lengths = {}
    tr_abundance = {}
    sum_exp = 0
    sum_len = 0
    for lines in open(exp_file):
        len_field = 1  #original: len_field = 1
        ab_field = 2 #original: ab_field = 7
        fields=lines.strip().split();
        if lines[0]=='#' or fields[len_field]==0:
            continue;
        if len_field != -1:
            tr_lengths[fields[0]] = int(fields[len_field])
            sum_len += int(fields[len_field])
        else:
            tr_lengths[fields[0]] = 1
            sum_len += 1
        tr_abundance[fields[0]] = float(fields[ab_field])
        sum_exp += float(fields[ab_field])

    print(sum_exp)

    rec_ab=0
    f=open(log_file,'r')
    lines = f.readlines()
    best_rec_dict = {}
    not_first_time = 0
    for (i,line) in enumerate(lines):
        if i<5:
            continue
        #pdb.set_trace() 
        #print('lenth',len(lines))
        tokens = line.split()
        org = tokens[9]; #need to select for the tokens[9][29:]  when having larger prefix
        rec = tokens[13]; rec_len = int(tokens[0]); tr_len = int(tokens[10])
        tr_ab = tr_abundance.get(org,0) / sum_exp
        z = best_rec_dict.get(org,[None,0]);                
        if rec_len>z[1]:
            best_rec_dict[org] = [rec,rec_len,tr_len,tr_ab/tr_len];
            if rec_len > 0.9*tr_len:
                rec_ab += tr_ab
    print('Fraction of abundance reconstuctedd: ' + str(rec_ab))



def false_positive(Rec_fasta,LongReads_rec_per,Dest_File):
    tr_dict = {}
    tr_matches = {}; tr_attributes = {}
    tr_ratio = {}; tr_att2 = {}
    curr_name = ''
    tot = 0
    for lines in open(Rec_fasta):
        tokens=lines.strip().split();
        #pdb.set_trace()
	#print(lines)
        if tokens[0][0]!='>':
	    clen = tr_dict.get(curr_name,[0,0]); clen = clen[1]
            tr_dict[curr_name] = [0,clen+len(tokens[0])]  #Code,Length
            tr_matches[curr_name] = 0
            tr_attributes[curr_name] = [0,0,clen+len(tokens[0]),'']  #matchSize, qSize, tSize, qName
            tr_att2[curr_name] = [0,0,clen+len(tokens[0]),'']
            tr_ratio[curr_name] = 0
            continue
        curr_name = tokens[0][1:]
        tot += 1

    rec = 0
    for lines in open(LongReads_rec_per):
        tokens=lines.strip().split();
        qName = tokens[9]; qSize = int(tokens[10]) #need to select for the tokens[9][29:]  when having larger prefix
        tName = tokens[13]; tSize = int(tokens[14])
        matchSize = int(tokens[0])
        if matchSize >= tr_matches.get(tName,0):
            tr_matches[tName] = matchSize
            tr_attributes[tName] = [matchSize, qSize, tSize, qName]

        if float(matchSize)/qSize >= tr_ratio.get(tName,0):
            tr_att2[tName] = [matchSize, qSize, tSize, qName]            

        if matchSize >= 0.9 * min(qSize,tSize):
            if tr_dict[tName][0] == 0:
                rec+=1
                tr_dict[tName][0] = 1

        if matchSize >= 0.9 * tSize:
            if tr_dict[tName][0] == 0:
                rec+=1
            tr_dict[tName][0] = 2
    print(str(rec)+','+str(tot))
    with open(Dest_File,'w') as write_file:
        for (tName,val) in tr_dict.iteritems():
            #write_file.write(tName+'\t'+str(val[0])+'\t'+str(val[1])+'\n')
            write_file.write(tName+'\t'+str(tr_attributes[tName][0])+'\t'+str(tr_attributes[tName][1])+ '\t' + str(tr_attributes[tName][2]) + '\t' + str(tr_attributes[tName][3])+'\n'+'\t'+str(tr_att2[tName][0])+'\t'+str(tr_att2[tName][1])+ '\t' + str(tr_att2[tName][2]) + '\t' + str(tr_att2[tName][3])+'\n')


def performance_plot(reconstr_log,trinity_log,plot_file,L,S,N):
    norm = 0
    for lines in open(reconstr_log):
        fields=lines.strip().split();
        if fields[0][0] != '#' and float(fields[3])>S:  #fragment size
            norm += float(fields[4])*float(fields[2])

    ab_list = [1,10,25,50,100,1e10]
    tot = [0] * len(ab_list)
    frac = [0] * len(ab_list)
    no = [0] * len(ab_list)
    length = [0] * len(ab_list)

    for lines in open(reconstr_log):
        fields=lines.strip().split();
        if fields[0][0] != '#' and float(fields[3])>S:  #fragment size
            #pdb.set_trace()
            tr_rec = float(fields[2]); tr_len = float(fields[3]); tr_ab = float(fields[4])
            for (i,ab) in enumerate(ab_list):
                if i >= len(ab_list)-1:
                    continue;
                tr_cov = tr_ab*L / norm * N 
                if tr_cov>ab_list[i] and tr_cov<=ab_list[i+1]:
                    no[i] +=1
                    length[i] += tr_len
                    tot[i] += min(tr_rec,tr_len)
                    if tr_rec > 0.95*tr_len:
                        frac[i]+=1
    trinity_tot = [0] * len(ab_list)
    trinity_frac = [0] * len(ab_list)
    trinity_no = [0] * len(ab_list)
    trinity_length = [0] * len(ab_list)

    for lines in open(trinity_log):
        fields=lines.strip().split();
        if fields[0][0] != '#' and float(fields[3])>S: #fragment size
            tr_rec = float(fields[2]); tr_len = float(fields[3]); tr_ab = float(fields[4])
            for (i,ab) in enumerate(ab_list):
                tr_cov = tr_ab*L / norm * N 
                if tr_cov>ab_list[i] and tr_cov<=ab_list[i+1]:
                    trinity_no[i] +=1
                    trinity_length[i] += tr_len
                    trinity_tot[i] += tr_rec 
                    if tr_rec > 0.95*tr_len:
                        trinity_frac[i]+=1

    with open(plot_file,'w') as plotFile:
        plotFile.write('Abundance \t No. Trans \t Frac. reconstructed (Ours) \t Frac. reconstructed (Trinity) \t Bases reconstructed (Ours) \t Bases reconstr (Trinity) \t No bases\n')
        for (i,ab) in enumerate(ab_list):
            plotFile.write(str(ab)+'\t'+str(no[i])+'\t'+str(frac[i])+'\t'+str(trinity_frac[i])+'\t'+str(tot[i])+'\t'+str(trinity_tot[i])+'\t' + str(length[i])+'\n')

def extract_isoforms(fasta_file,mummer_file,out_file):
    tr_dict = {}
    ctr = 1
    for lines in open(fasta_file):
        fields=lines.strip().split()
        tr_id = fields[0][1:]
        if fields[0][0]=='>':
            tr_dict[tr_id] = ctr;  ctr+=1

    gn= 0

    query_id = ''
    for lines in open(mummer_file):
        #pdb.set_trace()
        fields=lines.strip().split()
        if fields[0][0]=='>':
            target_id = fields[1][:]
            continue
            #pdb.set_trace()

        else:
            query_id = fields[0][:]

        if tr_dict.get(target_id) and tr_dict.get(query_id):
            gn +=1
            #pdb.set_trace()
            tr_dict[target_id]=min(tr_dict[target_id],tr_dict[query_id])
            tr_dict[query_id]=min(tr_dict[target_id],tr_dict[query_id])

    list_of_vals = tr_dict.values()
    val_dict = {}

    for (tr,val) in tr_dict.items():
        val_dict[val] = val_dict.get(val,0)+1

    
    with open(out_file,'w') as plotFile:
        plotFile.write('Name\t\t\t\t\tGroupId\tNIsoform\n')
        for (tr,val) in tr_dict.items():
            plotFile.write(str(tr)+'\t0\t0\t0\t0\t' + str(val) + '\t' + str(val_dict[val]) + '\n')





def performance_plot_express(oracle_fasta,num_isoform_file,reconstr_log,trinity_log,cufflinks_log,plot_file,exp_file,exp_format,L,S,N,iso_low,iso_high,use_oracle):
    #exp_file has to be a express output
    #N = 76885022
    tr_dict_full = {}
    #pdb.set_trace()
    
    if exp_format=='EXP':
        tot_ab = 0
        for lines in open(exp_file):
            fields=lines.strip().split();
            tr_id = fields[1].upper()
            if not fields[2].isdigit():
                continue;
            try:
                tr_len = float(fields[2]);  tr_ab = float(fields[6]) / N; tr_cov = float(fields[6]) / max(1,tr_len-L) * L # * 20 / 135
                tot_ab = tot_ab + tr_ab
                tr_dict_full[tr_id] = [tr_len, tr_cov, tr_ab, 0, 0, 0, 0] #last 4 entries are num_iso, ours,trinity,cufflinks
            except ValueError:
                continue
    elif exp_format == 'SIM':  #Type is simulation format
        #pdb.set_trace()
	tot_ab = 0
        for lines in open(exp_file):
            fields=lines.strip().split();
            tr_id = 'hg19_wgEncodeGencodeBasicV17_'+fields[0].upper(); tr_id = tr_id.upper()
            if not fields[1].isdigit():
                continue;
            try:
                tr_len = float(fields[1]);  tr_ab = float(fields[7]) ; tr_cov = float(fields[7])  # * 20 / 135
                tot_ab = tot_ab + tr_ab
                tr_dict_full[tr_id] = [tr_len, tr_cov, tr_ab, 0, 0, 0, 0] #last 4 entries are num_iso, ours,trinity,cufflinks
            except ValueError:
                continue
	if tot_ab == 0:
		tot_ab=1; # max(tot_ab,1);
        for (tr_id,vec) in tr_dict_full.iteritems():
            tr_dict_full[tr_id] = [vec[0],vec[1]/tot_ab*N*L/max(1,vec[0]),vec[2]/tot_ab,vec[3],vec[4],vec[5],vec[6]];
    elif exp_format=='KAL':
        tot_ab = 0
        for lines in open(exp_file):
            fields=lines.strip().split();
            tr_id = fields[0].upper()
            if not fields[1].isdigit():
                continue;
	    if float(fields[1]) < S:
		continue;
            try:
                tr_len = float(fields[1]);  tr_ab = float(fields[4])/1e6; tr_cov = float(fields[3]) / max(1,(tr_len-L)) * L # * 20 / 135
                tot_ab = tot_ab + tr_ab
                tr_dict_full[tr_id] = [tr_len, tr_cov, tr_ab, 0, 0, 0, 0] #last 4 entries are num_iso, ours,trinity,cufflinks
            except ValueError:
                continue



    #print('total abundnace:' + str(tot_ab))
    pdb.set_trace()
    '''Filter to files in oracle_fasta'''
    tr_os = {}
    for lines in open(oracle_fasta):
        fields=lines.strip().split()
        tr_id = fields[0][1:].upper()
        if fields[0][0]=='>':
            if tr_dict_full.get(tr_id):
                tr_os[tr_id] = tr_dict_full[tr_id]
    if not use_oracle:
       tr_os = tr_dict_full; #Use to bypass filtering

    pdb.set_trace()
    '''Filter based on number of isoforms'''
    #No of isoforms is specified 
    iso_lower = iso_low; iso_upper = iso_high; #consider only transcripts that have iso_lower<=num_iso <= iso_upper
    tr_dict = {}
    tot_ab = 0
    tot_no = 0
    ab_list = [0,1,10,25,50,100,1e10]
    no = [0] * len(ab_list)

    for lines in open(num_isoform_file):
        fields=lines.strip().split()
        tr_id = fields[0].upper()
	if exp_format=='SIM':
		tr_id = 'hg19_wgEncodeGencodeBasicV17_'+fields[0]
		tr_id = tr_id.upper()
        if tr_os.get(tr_id):
            num_iso = float(fields[6])
            if iso_lower<=num_iso and num_iso <= iso_upper:
                tr_dict[tr_id] = tr_os[tr_id]
                tot_no += 1
		tot_ab += tr_dict[tr_id][2]; tr_cov = tr_dict[tr_id][1]
                for (i,ab) in enumerate(ab_list):
                	if i >= len(ab_list)-1:
                    		continue;
                #tr_cov = tr_ab*L / norm * N 
                	if tr_cov>=ab_list[i] and tr_cov<ab_list[i+1]:
                    		no[i] +=1

    pdb.set_trace()
    #tr_dict = tr_os

    '''Run through our reconstructed file'''
    tot = [0] * len(ab_list)
    frac = [0] * len(ab_list)
    length = [0] * len(ab_list)
    ab_rec = 0
    no_rec = 0
    #tot_ab = 0

    for lines in open(reconstr_log):
        fields=lines.strip().split();
        tr_id = fields[0].upper()
        if (not tr_dict.get(tr_id)):
            continue
            #pdb.set_trace()
        tr_rec = float(fields[2]); tr_len = tr_dict[tr_id][0]; tr_cov = tr_dict[tr_id][1]; tr_ab = tr_dict[tr_id][2]
        #tot_ab += tr_ab
        if tr_id == 'ENST00000537877.1':
                        pass #pdb.set_trace()

        if 1:
            for (i,ab) in enumerate(ab_list):
                if i >= len(ab_list)-1:
                    continue;
                #tr_cov = tr_ab*L / norm * N 
                if tr_id == 'hg19_wgEncodeGencodeBasicV17_ENST00000560268.1':
			pdb.set_trace()
		if tr_cov>=ab_list[i] and tr_cov<ab_list[i+1]:
                    #no[i] +=1
                    length[i] += tr_len
                    tot[i] += min(tr_rec,tr_len)
                    if tr_rec > 0.9*tr_len:
                        tr_dict[tr_id][4]=1
                        frac[i]+=1
                        ab_rec += tr_ab
			no_rec += 1

    #pdb.set_trace()
    '''Run through trinity'''
    trinity_tot = [0] * len(ab_list)
    trinity_frac = [0] * len(ab_list)
    trinity_no = [0] * len(ab_list)
    trinity_length = [0] * len(ab_list)
    ab_trinity = 0
    no_trinity = 0

    for lines in open(trinity_log):
        fields=lines.strip().split();
        tr_id = fields[0].upper()
        if (not tr_dict.get(tr_id)):
            continue
            #pdb.set_trace()
        tr_rec = float(fields[2]); tr_len = tr_dict[tr_id][0]; tr_cov = tr_dict[tr_id][1]; tr_ab = tr_dict[tr_id][2]
        if 1:
            for (i,ab) in enumerate(ab_list):
                #tr_cov = tr_ab*L / norm * N 
                if tr_cov>=ab_list[i] and tr_cov<ab_list[i+1]:
                    trinity_no[i] +=1
                    trinity_length[i] += tr_len
                    trinity_tot[i] += tr_rec 
                    if tr_rec > 0.9*tr_len:
                        tr_dict[tr_id][5]=1
                        trinity_frac[i]+=1
                        ab_trinity += tr_ab
			no_trinity += 1
                        
    cufflinks_tot = [0] * len(ab_list)
    cufflinks_frac = [0] * len(ab_list)
    cufflinks_no = [0] * len(ab_list)
    cufflinks_length = [0] * len(ab_list)
    ab_cufflinks = 0
    no_cufflinks = 0

    for lines in open(cufflinks_log):
        fields=lines.strip().split();
        tr_id = fields[0].upper()
        if (not tr_dict.get(tr_id)):
            continue
            #pdb.set_trace()
        tr_rec = float(fields[2]); tr_len = tr_dict[tr_id][0]; tr_cov = tr_dict[tr_id][1]; tr_ab = tr_dict[tr_id][2]
        if 1:
            for (i,ab) in enumerate(ab_list):
                #tr_cov = tr_ab*L / norm * N 
                if tr_cov>=ab_list[i] and tr_cov<ab_list[i+1]:
                    cufflinks_no[i] +=1
                    cufflinks_length[i] += tr_len
                    cufflinks_tot[i] += tr_rec 
                    if tr_rec > 0.9*tr_len:
                        tr_dict[tr_id][6]=1
                        cufflinks_frac[i]+=1
                        ab_cufflinks += tr_ab
			no_cufflinks += 1



    with open(plot_file,'w') as plotFile:
        plotFile.write('Abundance \t No. Trans \t Frac. reconstructed (Ours) \t Frac. reconstructed (Trinity) \t Frac. reco (Cufflinks) \t Bases reconstructed (Ours) \t Bases reconstr (Trinity) \t No bases\n')
        for (i,ab) in enumerate(ab_list):
            plotFile.write(str(ab)+'\t'+str(no[i])+'\t'+str(frac[i])+'\t'+str(trinity_frac[i])+'\t' + str(cufflinks_frac[i])+'\t'+str(tot[i])+'\t'+str(trinity_tot[i])+'\t' + str(length[i])+'\n')

    with open(plot_file+'_all','w') as plotFile:
        for (tr_id, tr_info) in tr_dict.items():
            plotFile.write(str(tr_id)+'\t'+str(tr_info[0])+'\t'+str(tr_info[1])+'\t'+str(tr_info[2])+'\t'+str(tr_info[3])+'\t'+str(tr_info[4])+'\t' + str(tr_info[5])+'\t' + str(tr_info[6])+'\n')

    print('Total,Ours,Trinity,Cufflinks='+str(tot_ab)+','+str(ab_rec)+','+str(ab_trinity)+','+str(ab_cufflinks))
    print('Total,Ours,Trinity,Cufflinks='+str(tot_no)+','+str(no_rec)+','+str(no_trinity)+','+str(no_cufflinks))

def performance_plot_correct(reconstr_log,trinity_log,exp_file,plot_file,L,S,N):
    tr_abundance = {}
    sum_exp = 0
    sum_len = 0
    norm = 0

    for lines in open(exp_file):
        len_field = -1  #original: len_field = 1
        ab_field = 1 #original: ab_field = 7
        fields=lines.strip().split();
        if lines[0]=='#' or fields[len_field]==0:
            continue;
        if len_field != -1:
            tr_lengths[fields[0]] = int(fields[len_field])
            sum_len += int(fields[len_field])
        else:
            tr_lengths[fields[0]] = 1
            sum_len += 1
        tr_abundance[fields[0]] = float(fields[ab_field])
        sum_exp += float(fields[ab_field])
        

    ab_list = [1,10,25,50,100,1e10]
    tot = [0] * len(ab_list)
    frac = [0] * len(ab_list)
    no = [0] * len(ab_list)
    length = [0] * len(ab_list)



    norm = 0
    for lines in open(reconstr_log):
        fields=lines.strip().split();
        if fields[0][0] != '#' and float(fields[3])>S:  #fragment size
            norm += float(fields[4])*float(fields[2])


    for lines in open(reconstr_log):
        fields=lines.strip().split();
        if fields[0][0] != '#' and float(fields[3])>S:  #fragment size
            #pdb.set_trace()
            tr_rec = float(fields[2]); tr_len = float(fields[3]); tr_ab = float(fields[4])
            for (i,ab) in enumerate(ab_list):
                if i >= len(ab_list)-1:
                    continue;
                tr_cov = tr_ab*L / norm * N 
                if tr_cov>ab_list[i] and tr_cov<=ab_list[i+1]:
                    no[i] +=1
                    length[i] += tr_len
                    tot[i] += min(tr_rec,tr_len)
                    if tr_rec > 0.95*tr_len:
                        frac[i]+=1
    trinity_tot = [0] * len(ab_list)
    trinity_frac = [0] * len(ab_list)
    trinity_no = [0] * len(ab_list)
    trinity_length = [0] * len(ab_list)

    for lines in open(trinity_log):
        fields=lines.strip().split();
        if fields[0][0] != '#' and float(fields[3])>S: #fragment size
            tr_rec = float(fields[2]); tr_len = float(fields[3]); tr_ab = float(fields[4])
            for (i,ab) in enumerate(ab_list):
                tr_cov = tr_ab*L / norm * N 
                if tr_cov>ab_list[i] and tr_cov<=ab_list[i+1]:
                    trinity_no[i] +=1
                    trinity_length[i] += tr_len
                    trinity_tot[i] += tr_rec 
                    if tr_rec > 0.95*tr_len:
                        trinity_frac[i]+=1

    with open(plot_file,'w') as plotFile:
        plotFile.write('Abundance \t No. Trans \t Frac. reconstructed (Ours) \t Frac. reconstructed (Trinity) \t Bases reconstructed (Ours) \t Bases reconstr (Trinity) \t No bases\n')
        for (i,ab) in enumerate(ab_list):
            plotFile.write(str(ab)+'\t'+str(no[i])+'\t'+str(frac[i])+'\t'+str(trinity_frac[i])+'\t'+str(tot[i])+'\t'+str(trinity_tot[i])+'\t' + str(length[i])+'\n')


def reverse_performance_plot(reconstr_rev_log,trinity_rev_log,plot_file,L,S,N):
    no = 0; rec = 0;
    bases=0; rec_bases = 0;
    for lines in open(reconstr_rev_log):
        fields=lines.strip().split();
        if fields[0][0] != '#' and float(fields[3])>S:  #fragment size
            no +=1
            bases += float(fields[3])
            rec_bases += float(fields[2])
            #pdb.set_trace()
            if float(fields[2]) > 0.95*float(fields[3]):
                rec+=1
    trinity_no = 0; trinity_rec = 0;
    trinity_bases=0; trinity_rec_bases = 0;
    for lines in open(trinity_rev_log):
        fields=lines.strip().split();
        if fields[0][0] != '#' and float(fields[3])>S:  #fragment size
            trinity_no +=1
            trinity_bases += float(fields[3])
            trinity_rec_bases += float(fields[2])
            #pdb.set_trace()
            if float(fields[2]) > 0.95*float(fields[3]):
                trinity_rec+=1

    with open(plot_file,'a') as plotFile:
        plotFile.write('\n False positives \n')
        
        plotFile.write('iTrans reconstructed ' +  str(no) + ' transcripts, and out of this, ' + str(rec) +' of them have 0.95 overlap with real transcripts. \n') #' Fraction = ' + str(rec/no) + '\n')
        plotFile.write('Trinity reconstructed ' +  str(trinity_no) + ' transcripts, and out of this, ' + str(trinity_rec) +' of them have 0.95 overlap with real transcripts. \n') # ' Fraction = ' + str(trinity_rec/trinity_no))
        plotFile.write('iTrans reconstructed ' +  str(rec_bases) + ' bases correctly out of ' + str(bases) +' total bases in reconstructed transcripts. \n') #' Fraction = ' + str(rec/no) + '\n')
        plotFile.write('Trinity reconstructed ' +  str(trinity_rec_bases) + ' bases correctly out of ' + str(trinity_bases) +' total bases in reconstructed transcripts. \n') #' Fraction = ' + str(rec/no) + '\n')
        

def reverse_analyzer_blat(Repeats_File,Dest_File,fasta_file,N):
    f=open(Repeats_File,'r')
    lines = f.readlines()
    best_rec_dict = {}
    not_first_time = 0
    for (i,line) in enumerate(lines):
        if i<5:
            continue
        #pdb.set_trace() 
        #print('lenth',len(lines))
        tokens = line.split()
        org = tokens[9]; rec = tokens[13]; # use the threshold
        rec_len = int(tokens[0]); tr_len = int(tokens[10])
        z = best_rec_dict.get(org,[None,0]);                
        if rec_len>z[1]:
            best_rec_dict[org] = [rec,rec_len,tr_len];

    perf = 0
    sum_ret = 0
    with open(Dest_File,'w') as write_file:
        for (org,z) in best_rec_dict.iteritems():
            write_file.write(org+'\t'+str(z[0])+'\t'+str(z[1]) + '\t' + str(z[2]) + '\n') # + str(rec[4])+'\n')



def reverse_analyzer(Repeats_File,Dest_File,fasta_file,N):
    tr_lengths = {}
    identity = ' '
    for lines in open(fasta_file):
        fields=lines.strip().split();
        if fields[0][0]=='>':
            identity = fields[0][1:]
            #pdb.set_trace()
            continue
        tr_lengths[identity] = tr_lengths.get(identity,0)+int(len(fields[0]))
        #tr_abundance[fields[0]] = float(fields[7])
        #sum_exp += float(fields[7])
        #sum_len += int(fields[1])
    #print(tr_lengths)
    #raw_input()
  
    f=open(Repeats_File,'r')
    lines = f.readlines()
    best_rec = []
    not_first_time = 0
    for (i,line) in enumerate(lines):
        #pdb.set_trace() 
        #print('lenth',len(lines))
        tokens = line.split()
        if tokens[0] == ">":
            #pdb.set_trace()
            if not_first_time:
                sorted_list = sorted(tempdict.iteritems(), key=operator.itemgetter(1),reverse=True)
                tr_len = tr_lengths.get(identity,-1)
                #tr_ab = tr_abundance.get(identity,0)
                #print(sorted_list)
                #raw_input()
                if len(sorted_list) > 0:
                    best_rec.append([identity, sorted_list[0][0],sorted_list[0][1],tr_len])
                else:
                    best_rec.append([identity, 'None', 0, tr_len])
            not_first_time = 1 
            identity = tokens[1]
            tempdict  = {}
        else:
            tempdict[tokens[0]] = tempdict.get(tokens[0],0)+int(tokens[3])
    perf = 0
    sum_ret = 0
    with open(Dest_File,'w') as write_file:
        for rec in best_rec:
            write_file.write(str(rec[0])+'\t'+str(rec[1])+'\t'+str(rec[2]) + '\t' + str(rec[3]) + '\n') # + str(rec[4])+'\n')




def pre_process_dual(in_file,out_file):
    f = open(out_file,'w')
    for line in open(in_file):

        fields=line.strip().split();
        #pdb.set_trace()
        if fields[0][0]=='>' and len(fields)>2 and fields[2]=='Reverse':
            pass
        else:
            f.write(line)

def main():
    import sys
    parallel_blat(sys.argv[1],sys.argv[2],sys.argv[3])




#analyzer("./OneHalf_L_100_N_1000000_algo_output/reconstr_per.txt","./OneHalf_L_100_N_1000000_algo_output/reconstr_log.txt","./OneHalf_L_100_N_1000000_algo_input/random_out.exp",1000000)

#false_positive('./WingWongTest_algo_output/TRANS/transabyss-final.fa','./WingWongTest_algo_output/TRANS/trans_pacbio_per.txt','abc')

