import os, pdb
import sys, time
import os.path, tester
from set_exp import set_exp 
from filter_trans import filter_trans

L = 100
N = 8743351

#sn='SRR059162'
sn='SRR453566'
sparsity = 0.5
start_loc = 1
stop_loc = 120000000

'''sn='ACAN'
sparsity=1
start_loc = 89346673 - 100
stop_loc = 89418584 + 100'''

paired_end = 0
add_errors = 1
double_stranded = 1

sim = 0
to_set_exp = 0
generate_reads = 0
trimmomatic = 0 #currently not enabled
run_seecer =0
run_jellyfish =0
run_extension_corr =0
run_cpp = 0
use_cpp = 0
mb = 0
sparse_flow = 0
parallelize_sf = 0
generate_ref = 0
compare_ans = 0
run_trinity = 0
compare_trinity = 0
plots = 0
plots_express = 1
run_rsem_eval = 0
run_cuffinks = 0
compare_cufflinks = 0
run_stringtie = 0
compare_stringtie = 0
compare_refShannon = 0


run_oasis = 0
run_express = 0

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

if paired_end:
	F = 350 #Fragment size
	F_sd = 0 #Fragment size Standard Deviation
	sn = sn+'_F_'+str(F)
	#sn = sn + '_PE'
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
	mb_string += ' -d '
else:
	ds_string = ' --stranded '
	trinity_string += ' --SS_lib_type FR '


def run_cmd(s1):
	print(s1); os.system(s1)

#sample_name = sn + '_L_'+str(L)+'_N_'+str(N)+'_'  #Always use underscore at end
#sample_name = 'SRR059162_'
#sample_name = 'SRR453566_bubble_'
#sample_name = 'WW_chr1_'
#sample_name = 'WW_NM_004930_'
#sample_name = 'WW_NR_109896_'
#sample_name = 'WW_NM_005427_'
#sample_name = 'WingWongTest_K24_'
#sample_name = 'WW_aggressive_'
#sample_name = 'WingWong604_'
#sample_name = 'WingWongTest_'
#sample_name = 'OneHalf_gen5_ERR_L_100_N_1000000_'
#sample_name = 'WingWongTest_K24_retreive_'
#sample_name = 'Tebola_'
#sample_name = 'Snyder/10M_1_'
#sample_name = 'Ebola/Full_'
#sample_name = 'D250_'
#sample_name = 'OR_'
#sample_name = 'Ebola_'
#sample_name = 'WW_100_chr15_'
#sample_name = 'WingWong_20M_'
#sample_name = 'Snyder_chr15_'
#sample_name = 'Snyder_chr15_'
#sample_name = 'S15_Quorum_'
#sample_name = 'S15_PE_'
#sample_name = 'S15_New2_'
#sample_name = 'WingWong_Uniform_R100_K25_'
#sample_name = 'S15_SE_try_'
#sample_name = 'Snyder/PE_'
#sample_name = 'WingWongTest_'
#sample_name = 'S15_p6/S15_p6_c1_1_'

#sample_name='NM_Test_r2_'
#sample_name = 'Snyder_chr15_merged_'
#sample_name = 'Snyder_group2/2Kmer_'
#sample_name = 'SPombe_'
sample_name = 'Snyder/Analysis_'
#sample_name = 'WW_Analysis_'
#sample_name = 'Snyder/Quorum_'
#sample_name = 'Snyder/QBlat_'
#sample_name = 'WW_NM_005427_Metis_remaining_'
#sample_name = 'WW_check_'
#sample_name = 'S15_Metis_p1r1_g2_remaining_'
#sample_name = 'S15_'
#sample_name = 'Quorum_20M_Trinity_'
#sample_name = 'WingWongTest_K24_retreive_rerun_'
#sample_name = 'Dummy_'
#sample_name = "S15_r_control/S15_r_control_c1_3"
#sample_name = "S15_r3_control_c1_3"

sample_output_name = sample_name#  + 'Try_' 

base_dir = './' #~/Full_Assembler/'
bed_file=' ./Genome/GSE51861_isoform.bed'
ref_file=' ./Genome/hg19.fa' # reference chromosome
exp_file='./' + sample_name + 'algo_input/random_out.exp'
reads_file = ' ./' + sample_name + 'algo_input/reads'
if paired_end:
	reads_string = reads_file + '_1.fasta '+reads_file + '_2.fasta '
else:	
	reads_string = reads_file + '.fasta '

timer = {}
timer['start'] = time.time()	

if not os.path.exists('./'+sample_name+'algo_input'):
   	os.makedirs('./'+sample_name+'algo_input')
if not os.path.exists('./'+sample_output_name+'algo_output'):
   	os.makedirs('./'+sample_output_name+'algo_output')
if not os.path.exists('./'+sample_output_name+'intermediate'):
   	os.makedirs('./'+sample_output_name+'intermediate')

if mb or sparse_flow:
        os.system('cp *.py ./' + sample_output_name + 'algo_output/')


if generate_reads:
   	if to_set_exp:
		set_exp(start_loc,stop_loc,sparsity,'loguniform',sample_name)
	#print('python ./Read_Simulator/gensimreads.py ' + ds_string + ' -e ' + exp_file + ' -n ' + str(N)+ ' -l ' + str(L) + ' ' + bed_file + ' > ' + reads_file + '.bed')
	#run_cmd('python ./Read_Simulator/gensimreads.py ' + ds_string + ' -e ' + exp_file + ' -n ' + str(N)+ ' -l ' + str(L) + pairedend_string + bed_file + ' > ' + reads_file + '.bed')
	#print('python ./Read_Simulator/getseqfrombed.py ' + err_string + reads_file + '.bed ' +  ref_file + ' > ' + reads_file + '.fa')
	run_cmd('python ./Read_Simulator/getseqfrombed.py ' + err_string + reads_file + '.bed ' + ref_file + ' > ' + reads_file + '.fa')

	if paired_end:	
		run_cmd('python ./Read_Simulator/splitfasta_file.py -o ' + reads_file)
		run_cmd( 'awk \'/^>/{print s? s"\\n"$0:$0;s="";next}{s=s sprintf("%s",$0)}END{if(s)print s}\' ' +reads_file+ '_1.fa > '+ reads_file +'_1.fasta')
		run_cmd( 'awk \'/^>/{print s? s"\\n"$0:$0;s="";next}{s=s sprintf("%s",$0)}END{if(s)print s}\' ' +reads_file+ '_2.fa > '+ reads_file +'_2.fasta')
	else:
		'''highly optional: to write a fasta file for reference'''
		run_cmd( 'awk \'/^>/{print s? s"\\n"$0:$0;s="";next}{s=s sprintf("%s",$0)}END{if(s)print s}\' ' +reads_file+ '.fa > '+ reads_file +'.fasta')

timer['after_gen_reads'] = time.time() 
timer['for_gen_reads'] = timer['after_gen_reads'] - timer['start']

if trimmomatic:
	run_cmd('java -jar /home/sreeramkannan/Packages/trinityrnaseq_r20131110/trinity-plugins/Trimmomatic-0.30/trimmomatic.jar SE -phred33 SRR453566_1.fastq reads_trimmomatic.fastq ILLUMINACLIP:/home/sreeramkannan/Packages/trinityrnaseq_r20131110/trinity-plugins/Trimmomatic-0.30/adapters/TruSeq2-SE.fa:2:30:10 LEADING:20 TRAILING:20 SLIDINGWINDOW:4:30 MINLEN:75 CROP:75 TOPHRED64')
	run_cmd('java -jar /home/sreeramkannan/Packages/trinityrnaseq_r20131110/trinity-plugins/Trimmomatic-0.30/trimmomatic.jar PE -phred33 reads_1.fastq reads_2.fastq reads_trim_1.fastq reads_trim_UP1.fastq reads_trim_2.fastq reads_trim_UP2.fastq LEADING:20 TRAILING:20 SLIDINGWINDOW:4:30 MINLEN:50 CROP:50 TOPHRED64')



if run_seecer:
	#run_cmd('bash ' + seecer_dir+ ' -t ' + base_dir+sample_name+'algo_input/seecer_tmp -k ' + str(K_value) + reads_string )
	run_cmd('bash ' + seecer_dir+ ' -t ' + base_dir+sample_name+'algo_input/seecer_tmp -k 16 ' + reads_string )
	if paired_end:
		run_cmd('mv ./'+sample_name+'algo_input/reads_1.fasta '+ './'+sample_name+'algo_input/reads_1_before_correction.fasta ')
		run_cmd('mv ./'+sample_name+'algo_input/reads_2.fasta '+ './'+sample_name+'algo_input/reads_2_before_correction.fasta ')
		run_cmd('mv ./'+sample_name+'algo_input/reads_1.fasta_corrected.fa '+'./'+sample_name+'algo_input/reads_1.fasta ')
		run_cmd('mv ./'+sample_name+'algo_input/reads_2.fasta_corrected.fa '+'./'+sample_name+'algo_input/reads_2.fasta ')

	else:
		run_cmd('mv ./'+sample_name+'algo_input/reads.fasta '+ './'+sample_name+'algo_input/reads_before_correction.fasta ')
		run_cmd('mv ./'+sample_name+'algo_input/reads.fasta_corrected.fa '+'./'+sample_name+'algo_input/reads.fasta ')


if run_jellyfish:
	run_jfs = ' '
	if double_stranded:
		run_jfs += ' --both-strands '
	run_cmd('rm ./'+sample_name+'algo_input/jelly*')  #Remove old jellyfish files
	#run_cmd(jellyfish_dir+' count -m ' + str(K_value)+ run_jfs+ ' -o ' + base_dir+sample_name+'algo_input/jellyfish_output -s 20000000 -c 4 -t 32 ' +reads_string)
	run_cmd(jellyfish_dir+' count -m ' + str(K_value+1) + run_jfs+ ' -o ' + base_dir+sample_name+'algo_input/jellyfish_p1_output -s 20000000 -c 4 -t 32 ' +reads_string)


	#pdb.set_trace()
	'''if os.path.isfile('./'+sample_name+'algo_input/jellyfish_output_1'):
		run_cmd(jellyfish_dir+' merge -o '+ base_dir+sample_name+'algo_input/jellyfish_output.jf ' + base_dir+sample_name+'algo_input/jellyfish_output\_*')
	else
		run_cmd('mv ' + base_dir+sample_name+'algo_input/jellyfish_output_0 ' +base_dir+sample_name+'algo_input/jellyfish_output.jf')'''

	if os.path.isfile('./'+sample_name+'algo_input/jellyfish_p1_output_1'):
		run_cmd(jellyfish_dir+' merge -o '+ base_dir+sample_name+'algo_input/jellyfish_p1_output.jf ' + base_dir+sample_name+'algo_input/jellyfish_p1_output\_*')
	else:
		run_cmd('mv ' + base_dir+sample_name+'algo_input/jellyfish_p1_output_0 ' +base_dir+sample_name+'algo_input/jellyfish_p1_output.jf')
	
	#run_cmd(jellyfish_dir+' dump -c -L ' + str(jellyfish_kmer_cutoff) + ' ' + base_dir+sample_name+'algo_input/jellyfish_output.jf > ' +base_dir+sample_name+'algo_input/kmer.dict')
	run_cmd(jellyfish_dir+' dump -c -t -L ' + str(jellyfish_kmer_cutoff) + ' ' + base_dir+sample_name+'algo_input/jellyfish_p1_output.jf > ' +base_dir+sample_name+'algo_input/k1mer.dict_org')
	if (not run_extension_corr) and double_stranded:
                tester.double_strandify(base_dir+sample_name+'algo_input/k1mer.dict_org', base_dir+sample_name+'algo_input/k1mer.dict')
        if (not run_extension_corr) and (not double_stranded):
		run_cmd('mv ' + base_dir+sample_name+'algo_input/k1mer.dict_org ' + base_dir+sample_name+'algo_input/k1mer.dict')	
if run_extension_corr:	
		if double_stranded:
			str_ec = ' -d '
		else: 
			str_ec = ' '
		run_cmd('python extension_correction_SIC.py ' + str_ec + base_dir+sample_name+'algo_input/k1mer.dict_org ' +base_dir+sample_name+'algo_input/k1mer.dict 3 25 12 0.5')

if run_jellyfish or run_extension_corr:
		run_cmd('pypy kp1mer_to_kmer.py ' + base_dir+sample_name+'algo_input/k1mer.dict ' + base_dir+sample_name+'algo_input/kmer.dict')


if run_cpp:
        run_cmd('tr \'\\t\' \'\\n\' <' + base_dir+sample_name+'algo_input/kmer.dict > ' + base_dir+sample_name+'algo_input/kmer.dict_2l ') 
        run_cmd('tr \'\\t\' \'\\n\' <' + base_dir+sample_name+'algo_input/k1mer.dict > ' + base_dir+sample_name+'algo_input/k1mer.dict_2l ')
        #run_cmd('tr \'  \' \'\\n\' <' + base_dir+sample_name+'algo_input/k1mer.dict > ' + base_dir+sample_name+'algo_input/k1mer.dict_2l ')
	run_cmd('./condenser ' + base_dir+sample_name+'algo_input/kmer.dict_2l '+base_dir+sample_name+'algo_input/k1mer.dict_2l  ' + base_dir+sample_name+'algo_input/ '  + str(K_value))	

if mb:
	#print('python multibridging.py -f '+reads_file+'.fasta ' + sample_name + 'intermediate')
	run_cmd('rm ./'+sample_output_name+'intermediate/*')
	jf_s = ' '; cpp_s = ' ';
	use_jellyfish = 0; #force set parameter for jellyfish
	if use_jellyfish:
		jf_s = ' -j '+ './'+sample_name+'algo_input/kmer.dict '+'./'+sample_name+'algo_input/k1mer.dict '
	if use_cpp:
		cpp_s = '-c ' +base_dir+sample_name+'algo_input/nodes.txt '+base_dir+sample_name+'algo_input/edges.txt '
        if not paired_end:                
		run_cmd(' pypy multibridging.py -f ' + mb_string + jf_s + cpp_s + '--kmer=' + str(K_value)+ reads_file+'.fasta ' + sample_output_name + 'intermediate ' ) # ' 2>&1 | tee ./' + sample_name + 'algo_input/log.txt')
	else:
		run_cmd('pypy multibridging.py -f '+ mb_string + jf_s + cpp_s + '--kmer=' + str(K_value) + reads_file+'_1.fasta '+reads_file+'_2.fasta ' + sample_output_name+ 'intermediate ') # 2>&1 | tee ./' + sample_name + 'algo_input/log.txt')

timer['after_mb'] = time.time()
timer['for_mb'] = timer['after_mb'] - timer['after_gen_reads']


if sparse_flow:
        reconstr_file = './'+sample_output_name+'algo_output/reconstructed.fasta'
        run_cmd('rm '+'./'+sample_output_name+'algo_output/reconstructed_comp_*.fasta')
        run_cmd('rm '+'./'+sample_output_name+'algo_output/reconstructed.fasta')
        
	run_cmd('python Algorithm_Fast_condense.py -1 '+ sample_output_name)
        ncomp = 0
        iteration_string = " "
        while os.path.isfile('./'+ sample_output_name + 'intermediate/nodes'+str(ncomp)+'.txt'):
                if 0: #ncomp==0:  and ncomp!=2:  #testing purpose
                	ncomp +=1
                	continue
                print('Component:',ncomp)
                #raw_input()
		if not parallelize_sf:
                	os.system('python Algorithm_Fast_condense.py ' + str(ncomp) + ' '+ sample_output_name)
                iteration_string += str(ncomp) + " "
                ncomp=ncomp+1
        
	if parallelize_sf:
		os.system('parallel python Algorithm_Fast_condense.py {} ' + sample_output_name+ " ::: " + iteration_string)
        os.system("cat " + reconstr_file + " " + sample_output_name+'algo_output/reconstructed_comp_*.fasta' +  " >> " + reconstr_file)
        filter_trans('./'+sample_output_name+'algo_output/reconstructed.fasta', './'+sample_output_name+'algo_output/reconstructed_short.fasta', 200)




timer['after_sp'] = time.time()
timer['for_sp'] = timer['after_sp'] - timer['after_mb']

if generate_ref:
	ref_source = './sim_input/chr15_gc_downloaded.fasta'
	curr_ref = './' +sample_name+'algo_input/reference.fasta'
	from ref_generate import ref_generate
	ref_generate(ref_source,curr_ref,start_loc,stop_loc,exp_file,L)

timer['after_generate_ref'] = time.time()
timer['for_ref'] = timer['after_generate_ref'] - timer['after_sp']

if compare_ans:
	curr_ref = './' +sample_output_name+'algo_output/reference.fasta'
	reconstr = './' +sample_output_name+'algo_output/reconstructed.fasta'
	all_seq = './' +sample_output_name+'algo_output/all_sequences.fasta'
	reconstr_per = './' +sample_output_name+'algo_output/reconstr_per.txt'
	reconstr_rev_per = './' +sample_output_name+'algo_output/reconstr_rev_per.txt'
	reconstr_log = './' +sample_output_name+'algo_output/reconstr_log.txt'
	reconstr_rev_log = './' +sample_output_name+'algo_output/reconstr_fp_log.txt'
	all_seq_per = './' +sample_output_name+'algo_output/all_seq_per.txt'
	self_repeats = './' +sample_output_name+'algo_output/self_repeats.txt'
	if not blat:
		run_cmd('mummer -maxmatch -l 80 ' + reconstr + ' ' + curr_ref + ' > ' + reconstr_per)
		tester.analyzer(reconstr_per,reconstr_log,exp_file,N)
		#run_cmd('mummer -maxmatch -l 80 ' + curr_ref + ' ' + curr_ref + ' > ' + self_repeats)
		run_cmd('mummer -maxmatch -l 80 ' + curr_ref + ' ' + reconstr + ' > ' + reconstr_rev_per)
		tester.reverse_analyzer(reconstr_rev_per,reconstr_rev_log,reconstr,N)
	else:
		run_cmd('python parallel_blat.py ' + reconstr + ' ' + curr_ref + ' ' + reconstr_per)
		pass
		#tester.analyzer_blat(reconstr_per,reconstr_log,exp_file,N)
		tester.analyzer_blat_noExp(reconstr_per,reconstr_log,exp_file,N)
		#tester.analyzer_blat_noExp( './' +sample_name+'algo_output/cufflinks_pacBio.psl','./' +sample_name+'algo_output/cufflinks_pacBio_log.txt',exp_file,N)
		#run_cmd('blat ' + curr_ref + ' ' + reconstr + '  ' + reconstr_rev_per)
		#tester.reverse_analyzer_blat(reconstr_rev_per,reconstr_rev_log,reconstr,N)
	if false_positive:
		long_ref = './' +sample_output_name+'algo_output/reference_long.fasta'
		print('I am here')
		reconstr_long_per = './' +sample_output_name+'algo_output/reconstr_long_per.txt'
		#run_cmd('python parallel_blat.py ' + reconstr + ' ' + long_ref + ' ' + reconstr_long_per)
		tester.false_positive(reconstr,reconstr_long_per,reconstr_rev_log)	
	


#analyzer("./OneHalf3_L_100_N_1000000_algo_output/Trinity/trinity_rev_per_dual.txt","./OneHalf3_L_100_N_1000000_algo_output/Trinity/trinity_rev_log.txt","./OneHalf3_L_100_N_1000000_algo_output/Trinity/trinity.fasta",1000000)



if run_express:
        curr_ref = './' +sample_output_name+'algo_output/reference.fasta'
	ref_index = './' +sample_output_name+'algo_output/ref_index'
	sam_file = './' +sample_output_name+'algo_output/hits.sam'
	#run_cmd('bowtie-build --offrate 1 ' + curr_ref +  ' ' + ref_index)
	if not paired_end:
		run_cmd('bowtie -p 30 -aS -X 800 --offrate 1 ' +ref_index + ' -f ' + reads_file + '.fasta > ' + sam_file)
	else:
		 run_cmd('bowtie -p 30 -aS -X 800 --offrate 1 ' +ref_index + ' -f -1 ' + reads_file + '_1.fasta -2 '+  reads_file + '_2.fasta > '+ sam_file)
	run_cmd('express ' + curr_ref + ' ' + sam_file) 


if run_trinity or compare_trinity:
	curr_ref = './' +sample_output_name+'algo_output/reference.fasta'
	run_cmd('mkdir '+ sample_output_name+'algo_output/Trinity/' )
	trinity_per = './' +sample_output_name+'algo_output/Trinity/trinity_per.txt'
	trinity_fasta = './' +sample_output_name+'algo_output/Trinity/Trinity.fasta'
	trinity_per_dual = './' +sample_output_name+'algo_output/Trinity/trinity_per_dual.txt'
	trinity_rev_per = './' +sample_output_name+'algo_output/Trinity/trinity_rev_per.txt'
	trinity_rev_per_dual = './' +sample_output_name+'algo_output/Trinity/trinity_rev_per_dual.txt'
	trinity_log = './' +sample_output_name+'algo_output/Trinity/trinity_log.txt'
	trinity_rev_log = './' +sample_output_name+'algo_output/Trinity/trinity_rev_log.txt'
	if not paired_end:
		if run_trinity:
			run_cmd('Trinity --seqType fa --JM 4G --CPU 25 --single ' + base_dir  + sample_name + 'algo_input/reads.fasta --output ' + base_dir + sample_output_name+ 'algo_output/Trinity/')
	else:
		if run_trinity:
			run_cmd('Trinity --seqType fa  --SS_lib_type FR --JM 4G --CPU 25 --left ' + base_dir  + sample_name + 'algo_input/reads_1.fasta --right ' + base_dir  + sample_name + 'algo_input/reads_2.fasta --output ' + base_dir + sample_output_name+ 'algo_output/Trinity/')

if compare_trinity:			
	if not blat:
		run_cmd('mummer -maxmatch -l 80 -b ' + trinity_fasta + ' ' + curr_ref + ' > ' + trinity_per)
		tester.pre_process_dual(trinity_per,trinity_per_dual)
		tester.analyzer(trinity_per_dual,trinity_log,exp_file,N)
		run_cmd('mummer -maxmatch -l 80 -b ' + curr_ref + ' ' + trinity_fasta + ' > ' + trinity_rev_per)
		tester.pre_process_dual(trinity_rev_per,trinity_rev_per_dual)
		tester.reverse_analyzer(trinity_rev_per_dual,trinity_rev_log,trinity_fasta,N)
	else:
		run_cmd('python parallel_blat.py ' + trinity_fasta + ' ' + curr_ref + ' ' + trinity_per)
		tester.analyzer_blat_noExp(trinity_per,trinity_log,exp_file,N)
		#run_cmd('blat ' + curr_ref + ' ' + trinity_fasta + ' ' + trinity_rev_per)
		#tester.reverse_analyzer_blat(trinity_rev_per,trinity_rev_log,trinity_fasta,N)
		pass
        if false_positive:
                long_ref = './' +sample_output_name+'algo_output/reference_long.fasta'
                trinity_long_per = './' +sample_output_name+'algo_output/Trinity/trinity_long_per.txt'
                run_cmd('python parallel_blat.py ' + trinity_fasta + ' ' + long_ref + ' ' + trinity_long_per)
                tester.false_positive(trinity_fasta,trinity_long_per,trinity_rev_log)
	



if run_cuffinks or compare_cufflinks:
	curr_ref = './' +sample_name+'algo_input/reference.fasta'
	genome_index = './Genome/Homo_sapiens/UCSC/hg19/Sequence/Bowtie2Index/genome '
	genome_index = './Genome/hg19 '
	cufflinks_dir = sample_output_name+'algo_output/Cufflinks/'
	run_cmd('mkdir '+ cufflinks_dir)
	cufflinks_per = './' +sample_output_name+'algo_output/Cufflinks/cufflinks_per.txt'
	cufflinks_fasta = './' +sample_output_name+'algo_output/Cufflinks/cufflinks.fasta'
	cufflinks_log = './' +sample_output_name+'algo_output/Cufflinks/cufflinks_log.txt'
        cufflinks_rev_log = './' +sample_output_name+'algo_output/Cufflinks/cufflinks_rev_log.txt'
	


if run_cuffinks:
	run_cmd('tophat -p 20 --no-coverage-search -o ' + cufflinks_dir + ' '+ genome_index + reads_string )
	run_cmd('cufflinks -p 20 -o ' + cufflinks_dir + ' ' + cufflinks_dir + 'accepted_hits.bam')
	run_cmd('gffread -w ' + cufflinks_fasta + ' -g ./Genome/hg19.fa ' + cufflinks_dir + 'transcripts.gtf')
	
if compare_cufflinks:
	#run_cmd('python parallel_blat.py ' + cufflinks_fasta + ' ' + curr_ref + ' ' + cufflinks_per)
	#tester.analyzer_blat_noExp(cufflinks_per,cufflinks_log,exp_file,N)
	if false_positive:
                long_ref = './' +sample_output_name+'algo_output/reference_long.fasta'
                cufflinks_long_per = './' +sample_output_name+'algo_output/Cufflinks/cufflinks_long_per.txt'
                run_cmd('python parallel_blat.py ' + cufflinks_fasta + ' ' + long_ref + ' ' + cufflinks_long_per)

                tester.false_positive(cufflinks_fasta,cufflinks_long_per,cufflinks_rev_log)
    


if run_stringtie or compare_stringtie:
        curr_ref = './' +sample_name+'algo_input/reference.fasta'
        genome_index = './Genome/Homo_sapiens/UCSC/hg19/Sequence/Bowtie2Index/genome '
        genome_index = './Genome/hg19 '
        stringtie_dir = sample_output_name+'algo_output/Cufflinks/'
        run_cmd('mkdir '+ stringtie_dir)
        stringtie_per = './' +sample_output_name+'algo_output/Cufflinks/stringtie_per.txt'
        stringtie_fasta = './' +sample_output_name+'algo_output/Cufflinks/stringtie.fasta'
        stringtie_log = './' +sample_output_name+'algo_output/Cufflinks/stringtie_log.txt'
        stringtie_rev_log = './' +sample_output_name+'algo_output/Cufflinks/stringtie_rev_log.txt'


if run_stringtie:
        #run_cmd('tophat -p 20 --no-coverage-search -o ' + stringtie_dir + ' '+ genome_index + reads_string )
        run_cmd('~/Packages/stringtie-1.0.4.Linux_x86_64/stringtie ' + stringtie_dir + 'accepted_hits.bam -f 0 -p 20 -o ' + stringtie_dir + 'stringtie.gtf ')
        run_cmd('gffread -w ' + stringtie_fasta + ' -g ./Genome/hg19.fa ' + stringtie_dir + 'stringtie.gtf')

if compare_stringtie:
        #run_cmd('python parallel_blat.py ' + stringtie_fasta + ' ' + curr_ref + ' ' + stringtie_per)
        #tester.analyzer_blat_noExp(stringtie_per,stringtie_log,exp_file,N)
        if false_positive:
                long_ref = './' +sample_output_name+'algo_output/reference_long.fasta'
                stringtie_long_per = './' +sample_output_name+'algo_output/cufflinks_long_per.txt'
                run_cmd('python parallel_blat.py ' + stringtie_fasta + ' ' + long_ref + ' ' + stringtie_long_per)

                tester.false_positive(stringtie_fasta,stringtie_long_per,stringtie_rev_log)




if compare_refShannon:
        curr_ref = './' +sample_output_name+'algo_output/reference.fasta'
	long_ref = './' +sample_output_name+'algo_output/reference_long.fasta'
        rs_reconstr = './' +sample_output_name+'algo_output/rs_reconstructed.fasta'
	rs_reconstr_long_per = './' +sample_output_name+'algo_output/rs_reconstr_long_per.txt'
        rs_reconstr_per = './' +sample_output_name+'algo_output/rs_reconstr_per.txt'
        rs_reconstr_rev_per = './' +sample_output_name+'algo_output/rs_reconstr_rev_per.txt'
        rs_reconstr_log = './' +sample_output_name+'algo_output/rs_reconstr_log.txt'
        rs_reconstr_rev_log = './' +sample_output_name+'algo_output/rs_reconstr_fp_log.txt'
	#run_cmd('python parallel_blat.py ' + rs_reconstr + ' ' + curr_ref + ' ' + rs_reconstr_per)
	#tester.analyzer_blat_noExp(rs_reconstr_per,rs_reconstr_log,exp_file,N)
	if false_positive:
		run_cmd('python parallel_blat.py ' + rs_reconstr + ' ' + long_ref + ' ' + rs_reconstr_long_per)
		tester.false_positive(rs_reconstr,rs_reconstr_long_per,rs_reconstr_rev_log)


if compare_soap:
        curr_ref = './' +sample_name+'algo_input/reference.fasta'
	curr_ref = './' +sample_name+'algo_input/reference.fasta'

	genome_index = './Genome/hg19 '
        soap_dir = sample_output_name+'algo_output/SOAP/'
        run_cmd('mkdir '+ soap_dir)
        soap_per = './' +sample_output_name+'algo_output/SOAP/soap_per.txt'
        soap_fasta = './' +sample_output_name+'algo_output/SOAP/SOAP_graph.contig'
        soap_log = './' +sample_output_name+'algo_output/SOAP/soap_log.txt'
        soap_rev_log = './' +sample_output_name+'algo_output/SOAP/soap_rev_log.txt'
	

if compare_soap:
        run_cmd('python parallel_blat.py ' + soap_fasta + ' ' + curr_ref + ' ' + soap_per)
        tester.analyzer_blat_noExp(soap_per,soap_log,exp_file,N)
        if false_positive:
                tester.false_positive(soap_fasta,soap_per,soap_rev_log)


if run_oasis or compare_oasis:
        curr_ref = './' +sample_name+'algo_input/reference.fasta'
        genome_index = './Genome/hg19 '
        oasis_dir = sample_name+'algo_output/OASIS/'
        run_cmd('mkdir '+ oasis_dir)
        oasis_per = './' +sample_name+'algo_output/OASIS/oasis_per.txt'
        oasis_fasta = './' +sample_name+'algo_output/OASIS/transcripts.fa'
        oasis_log = './' +sample_name+'algo_output/OASIS/oasis_log.txt'
        oasis_rev_log = './' +sample_name+'algo_output/OASIS/oasis_rev_log.txt'
	if not paired_end:
                if run_oasis:
                        run_cmd('velveth ' + base_dir + sample_output_name+ 'algo_output/OASIS/ ' + str(K_value+1) + ' ' + base_dir  + sample_name + 'algo_input/reads.fasta')
			run_cmd('velvetg '  + base_dir + sample_output_name+ 'algo_output/OASIS/ -read_trkg yes')
			run_cmd('/home/sreeramkannan/Packages/oases_0.2.8/oases ' + base_dir + sample_output_name+ 'algo_output/OASIS/')
        else:
                if run_oasis:
                        run_cmd('velveth ' + base_dir + sample_output_name+ 'algo_output/OASIS/ ' + str(K_value+1) + ' -shortPaired -fasta -separate  ' + base_dir  + sample_name + 'algo_input/reads_1.fasta ' + base_dir  + sample_name + 'algo_input/reads_2.fasta ')
			run_cmd('velvetg '  + base_dir + sample_output_name+ 'algo_output/OASIS/ -read_trkg yes') 
			run_cmd('/home/sreeramkannan/Packages/oases_0.2.8/oases ' + base_dir + sample_output_name+ 'algo_output/OASIS/')



if compare_oasis:
        #run_cmd('python parallel_blat.py ' + oasis_fasta + ' ' + curr_ref + ' ' + oasis_per)
        #tester.analyzer_blat_noExp(oasis_per,oasis_log,exp_file,N)
        if false_positive:
                tester.false_positive(oasis_fasta,oasis_per,oasis_rev_log)


if compare_trans:
        curr_ref = './' +sample_name+'algo_input/reference.fasta'
        genome_index = './Genome/hg19 '
        trans_dir = sample_name+'algo_output/TRANS/'
        run_cmd('mkdir '+ trans_dir)
        trans_per = './' +sample_name+'algo_output/TRANS/trans_per.txt'
        trans_fasta = './' +sample_name+'algo_output/TRANS/transabyss-final.fa'
        trans_log = './' +sample_name+'algo_output/TRANS/trans_log.txt'
        trans_rev_log = './' +sample_name+'algo_output/TRANS/trans_rev_log.txt'


if compare_trans:
        run_cmd('python parallel_blat.py ' + trans_fasta + ' ' + curr_ref + ' ' + trans_per)
        tester.analyzer_blat_noExp(trans_per,trans_log,exp_file,N)
        if false_positive:
                tester.false_positive(trans_fasta,trans_per,trans_rev_log)


if run_rsem_eval:
	r=1.57895017697314
	p=0.00112891449768464
	curr_ref = './' +sample_name+'algo_input/reference.fasta '
	reconstr = './' +sample_name+'algo_output/reconstructed.fasta ' 
	rec_short = './' +sample_name+'algo_output/rec_short.fasta ' 
	trinity_fasta = './' +sample_name+'algo_output/Trinity/Trinity.fasta '
	cufflinks_fasta = './' +sample_name+'algo_output/cufflinks.fasta '
	read_file = base_dir  + sample_name + 'algo_input/reads.fasta ' 
	rsem_string = 'rsem-calculate-expression --no-qualities -p 20 ' # + str(r) + ' ' +str(p) + ' ' + str(L) + ' 0 '
	#rsem_string = 'rsem-calculate-expression --no-qualities --calc-evaluation-score ' + str(r) + ' ' +str(p) + ' ' + str(L) + ' 0 '
	out_dir = './' +sample_name+'algo_output/'
	
	if 1: #os.path.exists(curr_ref):
		#run_cmd('rsem-prepare-reference --bowtie '+curr_ref+  out_dir + 'true_ref')
		run_cmd(rsem_string + read_file + out_dir + 'true_ref ' + out_dir + 'true_ref')
	
	if 0:
		run_cmd('rsem-prepare-reference --no-polyA '+reconstr+ out_dir + 'our_rec')
		run_cmd(rsem_string + read_file + out_dir + 'our_rec ' + out_dir + 'our_rec')
	
	if 0: #os.path.exists(trinity_fasta):
		run_cmd('rsem-prepare-reference --no-polyA '+trinity_fasta+ out_dir + 'trinity')
		run_cmd(rsem_string + read_file + out_dir + 'trinity ' + out_dir+ 'trinity')

	
	if 0: #os.path.exists(cufflinks_fasta):
		run_cmd('rsem-prepare-reference --no-polyA '+cufflinks_fasta+ out_dir + 'cufflinks')
		run_cmd(rsem_string + read_file + out_dir + 'cufflinks ' + out_dir+ 'cufflinks')
	'''run_cmd('rsem-prepare-reference --no-polyA '+rec_short+ out_dir + 'our_short')
	run_cmd(rsem_string + read_file + out_dir + 'our_short ' + out_dir + 'our_short')'''

if plots:
	reconstr_log = './' +sample_name+'algo_output/reconstr_log.txt'
	trinity_log = './' +sample_name+'algo_output/Trinity/trinity_log.txt'
	plot_file = './' +sample_name+'algo_output/plot_our_trinity.txt'
	reconstr_rev_log = './' +sample_name+'algo_output/reconstr_rev_log.txt'
	trinity_rev_log = './' +sample_name+'algo_output/Trinity/trinity_rev_log.txt'
	tester.performance_plot(reconstr_log,trinity_log,plot_file,L,200,N)
        tester.performance_plot('./' +sample_name+'algo_output/reconstr_log.txt','./' +sample_name+'algo_output/Cufflinks/cufflinks_log.txt','./' +sample_name+'algo_output/plot_our_cufflinks.txt',L,200,N)	

	#tester.reverse_performance_plot(reconstr_rev_log,trinity_rev_log,plot_file,L,200,N)

timer['after_comp'] = time.time()
timer['for_comp'] = timer['after_comp'] - timer['after_generate_ref']

print(timer)


if plots_express:
	reconstr_log = './' +sample_output_name+'algo_output/reconstr_log.txt'
	trinity_log = './' +sample_output_name+'algo_output/Trinity/trinity_log.txt'
	cufflinks_log = './' +sample_output_name+'algo_output/Cufflinks/cufflinks_log.txt'
        trans_log = './' +sample_output_name+'algo_output/TRANS/trans_log.txt'
        soap_log =  './' +sample_output_name+'algo_output/SOAP/soap_log.txt'
        oasis_log = './' +sample_output_name+'algo_output/OASIS/oasis_log.txt'
        curr_ref = './' +sample_output_name+'algo_output/reference.fasta'
	#oracle_fasta = './' +sample_output_name+'algo_output/oracle_10.exp'
	oracle_fasta = './' +sample_output_name+'algo_output/kmer_oracle.fasta'
        #oracle_fasta = './' +sample_output_name+'algo_output/GSE_oracle_set.fasta'

	express_file = './' +sample_output_name+'algo_output/results.xprs'
	ref_shannon_log  = './' +sample_output_name+'algo_output/ref_shannon_log.txt'
	kal_file = './' +sample_output_name+'algo_output/abundance.tsv'
	num_isoform_file = './' +sample_output_name+'algo_output/reference.isoform'
        plotOthers_file = './' +sample_output_name+'algo_output/plots_noOracle_oasisSoapTrans_kmeroracle_isogeq1.txt'
	plotMain_file = './' +sample_name+'algo_output/plots_NEW_isoAll.txt'
	iso_high=1e10;iso_low=1;
	use_oracle=True
	if use_oracle:
		pass
	else:
		oracle_fasta = curr_ref
	if sim:
		tester.performance_plot_express(oracle_fasta,exp_file,reconstr_log,trinity_log,cufflinks_log,plot_file,exp_file,'SIM',L,200,N,iso_low,iso_high,use_oracle)
		#tester.performance_plot_express(oracle_fasta,exp_file,oasis_log,soap_log,trans_log,plot_file,exp_file,'SIM',L,200,N,iso_low,iso_high,use_oracle)

        else:
                tester.performance_plot_express(oracle_fasta,num_isoform_file,reconstr_log,trinity_log,trinity_log,plotMain_file,kal_file,'KAL',L,200,N,iso_low,iso_high,use_oracle)

		#tester.performance_plot_express(oracle_fasta,num_isoform_file,reconstr_log,trinity_log,cufflinks_log,plotMain_file,express_file,'EXP',L,200,N,iso_low,iso_high,use_oracle)
		#tester.performance_plot_express(oracle_fasta,num_isoform_file,oasis_log,soap_log,trans_log,plotOthers_file,express_file,'EXP',L,200,N,iso_low,iso_high,use_oracle)

	#tester.performance_plot(reconstr_log,trinity_log,plot_file,L,200,N)
	#tester.reverse_performance_plot(reconstr_rev_log,trinity_rev_log,plot_file,L,200,N)


if extract_bam:
	'''samtools view accepted_hits.bam chr10:5,932,166-5,979,558 | awk '{OFS="\t"; print ">"$1"\n"$10}' - > reads_NM_005427.fasta
	samtools view accepted_hits.bam chr15 | awk '{OFS="\t"; print ">"$1"\n"$10}' - > reads_chr15.fasta
	samtools view accepted_hits.bam chr15'''
	pass

if bowtie:
	run_cmd('bowtie-build ref.fasta ref_index')
	run_cmd('bowtie -p 20 -f -S pacbio_index ../WingWong_algo_input/reads.fasta reads_pacbio.sam')
	run_cmd('bowtie ref_index -1 ../SPombe_PE_algo_input/reads_trim_1.fastq -2 ../SPombe_PE_algo_input/reads_trim_2.fastq reconstructed_hits_PE_reads')

if oracle_set:
	run_cmd('samtools view -bS reads_pacbio.sam > reads_pacbio.bam')
	run_cmd('samtools depth reads_pacbio.bam ')


# Running ALLPATHS-LG ECC
# /usr/bin/allpathsLG/ErrorCorrectReads.pl PHRED_ENCODING=64 UNPAIRED_READS_IN=reads_trimmomatic.fastq READS_OUT=reads_allpaths.fastq THREADS=50 PLOIDY=1 KEEP_KMER_SPECTRA=1
