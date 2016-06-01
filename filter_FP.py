import os, sys
def run_cmd(s1):
		print(s1); 
		os.system(s1) # + '>> temp_out.txt')


def write_filtered_tr(depth_file, in_tr_file, out_tr_file, log_file):
        import sys, pdb
        tr_hits = {}
        THRESH = 0.5
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
	run_cmd('samtools sort '+  out_dir + '/rec.bam ' +  out_dir + '/rec_sort')
	run_cmd('samtools depth ' +  out_dir + '/rec_sort.bam > ' +  out_dir + '/rec.depth')

	#Use the BAM file along with our tool to get new fasta file
	depth_file = out_dir + '/rec.depth'; 
	in_tr_file = rec_fasta
	out_tr_file = out_dir + '/rec.fasta'; 
	log_file = out_dir + '/rec.log'
	write_filtered_tr(depth_file, in_tr_file, out_tr_file, log_file)
	run_cmd('cp ' + rec_fasta + ' ' + out_dir +'/reconstructed_org.fasta')
	run_cmd('mv ' + out_tr_file + ' ' + out_dir + '/reconstructed.fasta')


if __name__ == '__main__':
	rec_fasta = sys.argv[1]
	read_1 = sys.argv[2]
	read_2 = sys.argv[3]
	out_dir = sys.argv[4]
	flags = ' '.join(sys.argv[5:])
	filter_FP(rec_fasta, read_1, read_2, out_dir, flags)


