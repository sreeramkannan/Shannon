import os, pdb
import sys, time
import os.path, tester
import multibridging
#from set_exp import set_exp 
#from filter_trans import filter_trans
# This script is used to call all the the steps of the algorithm.

# List of parameters
def run_MB_SF(arguments,inMem=False,contigs=[],weights=[],rps=[]):
	print(arguments)
	inDisk = not inMem
	arguments = arguments.strip().split()
	L = 100
	N = 8743351

	sn = ''
	sparsity = 0.5
	start_loc = 1
	stop_loc = 120000000


	paired_end = 0
	add_errors = 1
	double_stranded = 0

	sim = 0
	to_set_exp = 0
	generate_reads = 0
	trimmomatic = 0 #currently not enabled
	run_seecer =0
	run_jellyfish =0        # runs jellyfish to get kmers
	run_extension_corr =0   # runs contig based error correction to filter k1mers used to build kmer graph
	run_cpp = 0             # builds condensed kmer graph
	mb = 0                  # runs multibridging
	sparse_flow = 0
	parallelize_sf = 0
	generate_ref = 0
	compare_ans = 0         # compares reconstructed transcripts to reference transcripts to measure performance of assembler
	run_trinity = 0         # Runs a competing algorithm trinity     
	compare_trinity = 0
	plots = 0
	plots_express = 0
	run_rsem_eval = 0
	run_cuffinks = 0        # Runs a competing algorithm cufflinks
	compare_cufflinks = 0
	compare_soap =0
	compare_oasis=0
	compare_trans=0

	false_positive=0


	bowtie = 0 
	extract_bam = 0
	oracle_set = 0


	err_string = ' '
	pairedend_string = ' '
	mb_string = ' '
	trinity_string = ' '
	jellyfish_dir = ' jellyfish'	
	jellyfish_kmer_cutoff = 1
	seecer_dir = ' /data/sreeramk/packages/SEECER-0.1.3/SEECER/bin/run_seecer.sh'
	K_value = 24
	blat = 1

	# ----------------------------------------------------------------------------------


	# Sets parameters from terminal
	sample_name = None
	shannon_dir = ''
	only_k1 = ' --only_k1 '  #Default: write only k1mers
	only_reads = False  #default: only_reads = false
	nJobs = 1
	python_path = 'python'
	n_inp = arguments
	if len(n_inp)>1:
	    sample_name = arguments[0]
	    if '--run_alg' in n_inp:
	        mb = 1
	        sparse_flow = 1
	    if '--ds' in n_inp:
	        double_stranded = 1
	    if '--paired_end' in n_inp:
	        paired_end = 1
	    if '--compare' in n_inp:
	        compare_ans = 1
	    if '--kmer_size' in n_inp:
	        K_value = int(n_inp[n_inp.index('--kmer_size')+1])
	    if '--only_reads' in n_inp:
	    	only_reads = True
	    	run_jellyfish = 1
	    	run_extension_corr = 1
	    if '--nJobs' in n_inp:
	        nJobs = int(n_inp[n_inp.index('--nJobs')+1])
	    if '--both' in n_inp:
	    	only_k1 == ' ' #Dont write only k1mers
	        #print(K_value)
	    if '--dir_name'in n_inp:
	        directory_name = n_inp[n_inp.index('--dir_name')+1]
	    if '--shannon_dir' in n_inp:
		shannon_dir = n_inp[n_inp.index('--shannon_dir')+1]
	    if '--python_path' in n_inp:
	    	python_path = n_inp[n_inp.index('--python_path')+1]
	        
	if paired_end:
		F = 350 #Fragment size
		F_sd = 0 #Fragment size Standard Deviation
		sn = sn+'_F_'+str(F)
		pairedend_string = ' -p ' + str(F) + ',' + str(F_sd) + ' '
	else:
		F=L


	if add_errors:
		error_rate = 0.01
		err_string = ' -r ' + str(error_rate) + ' '
		sn = sn + '_ERR'
		mb_string += ' -e '

	if double_stranded:
		sn = sn + '_DS' 
		ds_string = ' '
		#mb_string += ' -d '
	else:
		ds_string = ' --stranded '
		trinity_string += ' --SS_lib_type FR '


	def run_cmd(s1):
		os.system(s1)

	sample_output_name = sample_name

	base_dir = '' 
	bed_file=' ./Genome/GSE51861_isoform.bed'
	ref_file=' ./Genome/hg19.fa' # reference chromosome
	exp_file= sample_name + 'algo_input/random_out.exp'
	reads_file = sample_name + 'algo_input/reads'
	if paired_end:
		reads_string = reads_file + '_1.fasta '+reads_file + '_2.fasta '
	else:	
		reads_string = reads_file + '.fasta '

	timer = {}
	timer['start'] = time.time()	

	if not os.path.exists(sample_name+'algo_input'):
	   	os.makedirs(sample_name+'algo_input')
	if not os.path.exists(sample_output_name+'algo_output'):
	   	os.makedirs(sample_output_name+'algo_output')
	if not os.path.exists(sample_output_name+'intermediate'):
	   	os.makedirs(sample_output_name+'intermediate')

	'''if mb or sparse_flow:
	        os.system('cp *.py ' + sample_output_name + 'algo_output/')'''


	if trimmomatic:
		run_cmd('java -jar /home/sreeramkannan/Packages/trinityrnaseq_r20131110/trinity-plugins/Trimmomatic-0.30/trimmomatic.jar SE -phred33 SRR453566_1.fastq reads_trimmomatic.fastq ILLUMINACLIP:/home/sreeramkannan/Packages/trinityrnaseq_r20131110/trinity-plugins/Trimmomatic-0.30/adapters/TruSeq2-SE.fa:2:30:10 LEADING:20 TRAILING:20 SLIDINGWINDOW:4:30 MINLEN:75 CROP:75 TOPHRED64')
		run_cmd('java -jar /home/sreeramkannan/Packages/trinityrnaseq_r20131110/trinity-plugins/Trimmomatic-0.30/trimmomatic.jar PE -phred33 reads_1.fastq reads_2.fastq reads_trim_1.fastq reads_trim_UP1.fastq reads_trim_2.fastq reads_trim_UP2.fastq LEADING:20 TRAILING:20 SLIDINGWINDOW:4:30 MINLEN:50 CROP:50 TOPHRED64')



	if run_jellyfish:
		run_jfs = ' '
		if double_stranded:
			run_jfs += ' -C '
		#run_cmd('rm '+sample_name+'algo_input/jelly*')  #Remove old jellyfish files
		run_cmd(jellyfish_dir+' count -m ' + str(K_value+1) + run_jfs+ ' -o ' + sample_name+'algo_input/jellyfish_p1_output -s 20000000 -c 4 -t ' + str(nJobs) + ' ' +reads_string)

		'''if os.path.isfile(sample_name+'algo_input/jellyfish_p1_output_1'):
			run_cmd(jellyfish_dir+' merge -o '+ sample_name+'algo_input/jellyfish_p1_output.jf ' + sample_name+'algo_input/jellyfish_p1_output\_*')
		else:
			run_cmd('mv ' + sample_name+'algo_input/jellyfish_p1_output_0 ' +sample_name+'algo_input/jellyfish_p1_output.jf')'''
		
		run_cmd(jellyfish_dir+' dump -c -t -L ' + str(jellyfish_kmer_cutoff) + ' ' + sample_name+'algo_input/jellyfish_p1_output > ' + sample_name+'algo_input/k1mer.dict_org')
		if (not run_extension_corr) and double_stranded:
	                tester.double_strandify(sample_name+'algo_input/k1mer.dict_org', sample_name+'algo_input/k1mer.dict')
	        if (not run_extension_corr) and (not double_stranded):
			run_cmd('mv ' + sample_name+'algo_input/k1mer.dict_org ' + sample_name+'algo_input/k1mer.dict')	

	if run_extension_corr:	
			if double_stranded:
				str_ec = ' -d '
			else: 
				str_ec = ' '
			run_cmd(python_path + ' ext_corr.py ' + str_ec + sample_name+'algo_input/k1mer.dict_org ' +sample_name+'algo_input/k1mer.dict 3 75')

	'''if run_jellyfish or run_extension_corr:
			run_cmd('python kp1mer_to_ kmer.py ' + sample_name+'algo_input/k1mer.dict ' + sample_name+'algo_input/kmer.dict')'''


	if run_cpp:
	        run_cmd('tr \'\\t\' \'\\n\' <' + sample_name+'algo_input/kmer.dict > ' + sample_name+'algo_input/kmer.dict_2l ') 
	        run_cmd('tr \'\\t\' \'\\n\' <' + sample_name+'algo_input/k1mer.dict > ' + sample_name+'algo_input/k1mer.dict_2l ')
		run_cmd('./condenser ' + sample_name+'algo_input/kmer.dict_2l '+sample_name+'algo_input/k1mer.dict_2l  ' + sample_name+'algo_input/'  + ' ' + str(K_value) + ' | tee ' + sample_name + '_cpp_terminal_output.txt')	

	with open(sample_name + '_terminal_output.txt', 'w') as f7:
	    f7.write(" \n")

	    
	if inDisk and mb:
		#run_cmd('rm '+sample_output_name+'intermediate/*')
		jf_s = ' '; cpp_s = ' ';
		use_jellyfish = 1; use_cpp = 0 #force set parameter for jellyfish
		if use_jellyfish:
			jf_s = '  '+ sample_name+'algo_input/kmer.dict '+sample_name+'algo_input/k1mer.dict '
		if use_cpp:
			cpp_s = '-c ' +sample_name+'algo_input/nodes.txt '+sample_name+'algo_input/edges.txt '
	        if not paired_end:                
			run_cmd(python_path + ' ' + shannon_dir + 'multibridging.py  -f --kmer=' +str(K_value) + mb_string + only_k1 +  jf_s + cpp_s + reads_file+'.fasta ' + sample_output_name + 'intermediate ' + ' | tee ' + sample_name + '_terminal_output.txt') # ' 2>&1 | tee ./' + sample_name + 'algo_input/log.txt')
		else:
			run_cmd(python_path + ' ' + shannon_dir + 'multibridging.py -f --kmer='+ str(K_value) + mb_string + only_k1 + jf_s + cpp_s + reads_file+'_1.fasta '+reads_file+'_2.fasta ' + sample_output_name+ 'intermediate ' + ' | tee ' + sample_name + '_terminal_output.txt') # 2>&1 | tee ./' + sample_name + 'algo_input/log.txt')
	elif inMem and mb: 
		jf_s = ''; cpp_s = ''
		#In Memory
		if paired_end:
			multibridging.main('-f --kmer='+ str(K_value) + mb_string + only_k1 + jf_s + cpp_s + reads_file+'_1.fasta '+reads_file+'_2.fasta ' + sample_output_name+ 'intermediate ', inMem, contigs, weights, rps)
		else:
			multibridging.main('-f --kmer='+ str(K_value) + mb_string + only_k1 + jf_s + cpp_s + reads_file+'.fasta ' + sample_output_name+ 'intermediate ', inMem, contigs, weights, rps)

	timer['after_mb'] = time.time()
	#timer['for_mb'] = timer['after_mb'] - timer['after_gen_reads']


	if sparse_flow:
	        reconstr_file = sample_output_name+'algo_output/reconstructed.fasta'
	        #run_cmd('rm '+sample_output_name+'algo_output/reconstructed_comp_*.fasta')
	        #run_cmd('rm '+sample_output_name+'algo_output/reconstructed.fasta')
	        
	        run_cmd(python_path + ' ' + shannon_dir + 'algorithm_SF.py -1 '+ sample_output_name)
	        ncomp = 0
	        iteration_string = " "
	        while os.path.isfile(sample_output_name + 'intermediate/nodes'+str(ncomp)+'.txt'):
	                if 0: #ncomp==0:  and ncomp!=2:  #testing purpose
	                	ncomp +=1
	                	continue
	                print('Component:',ncomp)
	                if not parallelize_sf:
	                        os.system(python_path + ' '  + shannon_dir + 'algorithm_SF.py ' + str(ncomp) + ' '+ sample_output_name) # + ' | tee ' + sample_name + '_' + str(ncomp) + '_terminal_output.txt')
	                iteration_string += str(ncomp) + " "
	                ncomp=ncomp+1
	        
		if parallelize_sf:
			os.system('parallel ' + python_path + ' ' + shannon_dir + 'algorithm_SF.py {} ' + sample_output_name+ " ::: " + iteration_string)
	        os.system("cat " + sample_output_name+'algo_output/reconstructed_comp_*.fasta' +  " >> " + reconstr_file)
	        #filter_trans(sample_output_name+'algo_output/reconstructed.fasta', sample_output_name+'algo_output/reconstructed_short.fasta', 200)
	    
	if sparse_flow:
	    if os.path.exists(directory_name+"/before_sp_log.txt"):
	        f_log = open(directory_name+"/before_sp_log.txt", 'a')
	    else:
	        f_log = open(directory_name+"/before_sp_log.txt", 'w')
	        
	    num_transcripts = 0
	    with open(sample_output_name+"algo_output/reconstructed.fasta", 'r') as reconstructed_transcripts:
	        num_transcripts = len(reconstructed_transcripts.readlines())
	    

	    f_log.write(str(time.asctime()) + ": " +sample_output_name + " has completed: " + str(num_transcripts) + " transcripts" + "\n")
	    print(str(time.asctime()) + ": " +sample_output_name + " has completed: " + str(num_transcripts) + " transcripts")
	    f_log.close()



	timer['after_sp'] = time.time()
	#timer['for_sp'] = timer['after_sp'] - timer['after_mb']

	if compare_ans:
		curr_ref = sample_output_name+'algo_output/reference.fasta'
		reconstr = sample_output_name+'algo_output/reconstructed.fasta'
		all_seq = sample_output_name+'algo_output/all_sequences.fasta'
		reconstr_per = sample_output_name+'algo_output/reconstr_per.txt'
		reconstr_rev_per = sample_output_name+'algo_output/reconstr_rev_per.txt'
		reconstr_log = sample_output_name+'algo_output/reconstr_log.txt'
		reconstr_rev_log = sample_output_name+'algo_output/reconstr_rev_log.txt'
		all_seq_per = sample_output_name+'algo_output/all_seq_per.txt'
		self_repeats = sample_output_name+'algo_output/self_repeats.txt'
		if not blat:
			run_cmd('mummer -maxmatch -l 80 ' + reconstr + ' ' + curr_ref + ' > ' + reconstr_per)
			tester.analyzer(reconstr_per,reconstr_log,exp_file,N)
			run_cmd('mummer -maxmatch -l 80 ' + curr_ref + ' ' + reconstr + ' > ' + reconstr_rev_per)
			tester.reverse_analyzer(reconstr_rev_per,reconstr_rev_log,reconstr,N)
		else:
			run_cmd(python_path + ' ' + shannon_dir + 'parallel_blat_python.py ' + reconstr + ' ' + curr_ref + ' ' + reconstr_per)
			tester.analyzer_blat_noExp(reconstr_per,reconstr_log,exp_file,N)
		if false_positive:
			tester.false_positive(reconstr,reconstr_per,reconstr_rev_log)	
	

if __name__ == '__main__':

    if len(sys.argv) == 1:
        arguments = ['kmers.dict', 'allowed_kmers.dict', '1', '1', '-d']
    else:
        arguments = sys.argv[1:]
    run_MB_SF('\t'.join(arguments))

