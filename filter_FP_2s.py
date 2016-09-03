import os, sys
def run_cmd(s1):
		print(s1); 
		os.system(s1) # + '>> temp_out.txt')


def write_filtered_tr(depth_file, in_tr_file, out_tr_file, log_file):
        import sys, pdb
        tr_hits = {}
        THRESH = 0.9
        for a in open(depth_file): # as depth_file:
                #a=depth_file.readlines()
                fields = a.strip().split()
                tr_hits[fields[0]] = tr_hits.get(fields[0],0) +1;

        tr_len = {}; tr_name = ''
        with open(out_tr_file,'w') as tr_file, open(log_file,'w') as lf: #output transcripts
           for line in open(in_tr_file): #in_tr_file
                fields = line.strip().split()
                if fields[0][0] == '>': tr_name = fields[0][1:]; continue
                tr_len[tr_name]=len(fields[0]); trh = tr_hits.get(tr_name,0)
                lf.write(tr_name + '\t' + str(trh) + '\t' + str(len(fields[0])) + '\n')
                if tr_hits.get(tr_name,0) >= len(fields[0]) * THRESH:
                        tr_file.write('>' + tr_name + '\n')
                        tr_file.write(fields[0] + '\n')


def convert_vector(l1):
	nl = [];
	for i in l1:
		if nl and nl[-1][1]==(i-1): 
			nl[-1][1]=i; continue
		else:
			nl.append([i,i])
	return nl



def break_seqs(tr_l,tr_r,l1):
    K = 25
    insert = 200
    L = 100
    MIN_SEQ = 200 #Only write sequences longer than MIN_SEQ
    #l1 = len(seq)
    seqs = []
    if len(tr_l)==0 or len(tr_r)==0: return seqs
    #Add end buffers
    if tr_l[0][0] < K: tr_l[0][0] =1
    if tr_l[-1][1] > l1-insert: tr_l[-1][1] = l1
    if tr_r[0][0] < insert: tr_r[0][0]=1
    if tr_r[-1][1] < K : tr_r[-1][1] = l1
    l_max =len(tr_l)-1; r_max = len(tr_r)-1
    l_ptr = 0; r_ptr =0
    while(1):
        print(str(l_ptr) + ',' + str(r_ptr))
        if l_ptr > l_max or r_ptr > r_max: break
        curr_l = tr_l[l_ptr]; curr_r = tr_r[r_ptr]
        print(curr_l)
        print(curr_r)
        l1 = max(curr_l[0],curr_r[0])
        u1 = min(curr_l[1],curr_r[1])
        if l1+ MIN_SEQ<=u1: seqs.append([l1,u1]) 
        if curr_l[1] <= u1: 
            l_ptr +=1
        else:
            tr_l[l_ptr][0] = l1
        if curr_r[1] <= u1:
            r_ptr +=1
        else: 
            tr_r[r_ptr][0] = l1
    return seqs



def break_seqs_new(tr_l,tr_r,l1):
    K = 25
    insert = 200
    L = 100
    MIN_SEQ = 200 #Only write sequences longer than MIN_SEQ
    #l1 = len(seq)
    seqs = []
    if len(tr_l)==0 or len(tr_r)==0: return seqs
    #Add end buffers
    if tr_l[0][0] < K: tr_l[0][0] =1
    if tr_l[-1][1] > l1-insert: tr_l[-1][1] = l1
    if tr_r[0][0] < insert: tr_r[0][0]=1
    if tr_r[-1][1] < K : tr_r[-1][1] = l1
    l_max =len(tr_l)-1; r_max = len(tr_r)-1
    l_ptr = 0; r_ptr =0
    while(1):
        print(str(l_ptr) + ',' + str(r_ptr))
        if l_ptr > l_max or r_ptr > r_max: break
        curr_l = tr_l[l_ptr]; curr_r = tr_r[r_ptr]
        print(curr_l)
        print(curr_r)
        l1 = max(curr_l[0],curr_r[0])
        u1 = min(curr_l[1],curr_r[1])
        if l1<=u1: 
        	#There is intersection
        	new_l1 = min(curr_l[0],curr_r[0])
        	new_u1 = max(curr_l[1],curr_r[1])
        	seqs.append([new_l1,new_u1])
        	l_ptr +=1; r_ptr+=1
        else:
            #No intersection
            if curr_l[0] < l1: 
                l_ptr +=1
            if curr_r[0] < u1:
                r_ptr +=1
    return seqs

		


def read_in(depth_file):
        tr_dict = {}
        curr_name = ''; curr_list = []
        for line in open(depth_file): # as left_depth_file:
                fields = line.strip().split()
                if fields[0] == curr_name:
                	curr_list.append(int(fields[1]))
                else:
                	nl = convert_vector(curr_list)
                	if curr_name: tr_dict[curr_name]=nl
                	curr_name = fields[0]
                	curr_list = [int(fields[1])]
        nl = convert_vector(curr_list)
        if curr_name: tr_dict[curr_name]=nl
        return tr_dict



def write_new_tr(ld_file, rd_file, in_tr_file, out_tr_file, log_file):
        import sys, pdb
        MIN_SEQ = 200 #Only write sequences longer than MIN_SEQ
        tr_left = read_in(ld_file)
        tr_right = read_in(rd_file)


        tr_len = {}; tr_name = ''
        with open(out_tr_file,'w') as tr_file, open(log_file,'w') as lf: #output transcripts
           for line in open(in_tr_file): #in_tr_file
                fields = line.strip().split()
                if fields[0][0] == '>': tr_name = fields[0][1:]; continue
                tr_len[tr_name]=len(fields[0]); 
                bs = break_seqs_new(tr_left.get(tr_name,[]),tr_right.get(tr_name,[]),tr_len.get(tr_name,[]))
                seq = fields[0]
                ctr = 0
                for (l,u) in bs:
                	if l + MIN_SEQ <= u: 
                		aug = '' if (l==1 and u==len(seq)) else ('_frag_' + str(ctr))
                		tr_file.write('>' + tr_name + aug + '\n')
                		tr_file.write(seq[l-1:u] + '\n')
                		ctr +=1
                lf.write(tr_name + '\t' + str(bs) + '\n')

def filter_FP(rec_fasta, read_1, read_2, out_dir, flags = '-f --ff'):
	#take rec_fasta, read_1, read_2 files of type. 
	#Output in out_fasta. Temp dir = out_dir. 
	#Output fasta is in out_dir/reconstructed.fasta 
	#Flags is '-f --ff' for fasta and forward-forward
	#Flags is '-q --fr' for fastq and forward-reverse
	hisat_dir = '' #include dir/ if specifying directory
	
	#Build hisat file
	run_cmd(hisat_dir + 'hisat-build ' + rec_fasta + ' ' + out_dir + '/rec.hisat')

	#Align
	run_cmd(hisat_dir + 'hisat --no-spliced-alignment --no-discordant '+ flags + ' -x ' + out_dir + '/rec.hisat -1 ' + read_1 +  ' -2 ' + read_2 + ' -S ' + out_dir + '/rec.sam' )
	
	#Process SAM / BAM file to get depth information
	run_cmd('samtools view -bS -f 0x2 ' + out_dir + '/rec.sam > ' + out_dir + '/rec.bam')
	#run_cmd('samtools sort '+  out_dir + '/rec.bam ' +  out_dir + '/rec_sort')
	#run_cmd('samtools depth ' +  out_dir + '/rec_sort.bam > ' +  out_dir + '/rec.depth')
	run_cmd('samtools view -b -f 64 ' + out_dir +  '/rec.bam > ' + out_dir + '/c1.bam')
	run_cmd('samtools view -b -f 128 ' + out_dir +  '/rec.bam > ' + out_dir + '/c2.bam')
	run_cmd('samtools view -b -F 16 ' + out_dir +  '/c1.bam > ' + out_dir + '/c1f.bam')
	run_cmd('samtools view -b -f 16 ' + out_dir +  '/c1.bam > ' + out_dir + '/c1r.bam')
	run_cmd('samtools view -b -f 32 ' + out_dir +  '/c2.bam > ' + out_dir + '/c2r.bam')
	run_cmd('samtools view -b -F 32 ' + out_dir +  '/c2.bam > ' + out_dir + '/c2f.bam')
	run_cmd('samtools merge ' + out_dir +  '/l.bam ' + out_dir + '/c1f.bam ' +out_dir + '/c2r.bam ' )
	run_cmd('samtools merge ' + out_dir +  '/r.bam ' + out_dir + '/c2f.bam ' +out_dir + '/c1r.bam ' )
	run_cmd('samtools sort '+  out_dir + '/l.bam ' +  out_dir + '/l_sort')
	run_cmd('samtools sort '+  out_dir + '/r.bam ' +  out_dir + '/r_sort')
	run_cmd('samtools depth ' +  out_dir + '/l_sort.bam > ' +  out_dir + '/l.depth')
	run_cmd('samtools depth ' +  out_dir + '/r_sort.bam > ' +  out_dir + '/r.depth')


	#Use the BAM file along with our tool to get new fasta file
	ld_file = out_dir + '/l.depth'; rd_file = out_dir + '/r.depth'
	in_tr_file = rec_fasta
	out_tr_file = out_dir + '/rec.fasta'; 
	log_file = out_dir + '/rec.log'
	write_new_tr(ld_file, rd_file, in_tr_file, out_tr_file, log_file)

	run_cmd('cp ' + rec_fasta + ' ' + out_dir +'/reconstructed_org.fasta')
	run_cmd('mv ' + out_tr_file + ' ' + out_dir + '/reconstructed.fasta')


if __name__ == '__main__':
	rec_fasta = sys.argv[1]
	read_1 = sys.argv[2]
	read_2 = sys.argv[3]
	out_dir = sys.argv[4]
	#flags = ' '.join(sys.argv[5:])
	filter_FP(rec_fasta, read_1, read_2, out_dir)


