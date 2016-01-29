import time
import sys
import re
import pdb,math
import os
import os.path
import numpy as np
import tester
from filter_trans import filter_trans

from kmers_for_component import kmers_for_component
from weight_updated_graph import weight_updated_graph
from process_concatenated_fasta import process_concatenated_fasta
    
# For jellyfish
double_stranded = True
run_jellyfish = True
paired_end = False # Automatically set if command line is used

# General, Can be read in from terminal
reads_files = ['~/Full_Assembler/SPombe_algo_input/reads.fasta'] # ['./S15_SE_algo_input/reads.fasta']
sample_name = "SPombe_test"

# For extension correction
run_extension_corr =True
hyp_min_weight = 3
hyp_min_length = 75
partition_size = 500
use_partitioning = True

# For kmers_for_component
use_second_iteration = True
get_partition_kmers = True
overload = 2
K = 24
penalty = 5

# General runtime options
run_parallel = False
compare_ans = False


# Everything beyond this point does not need to be set
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 


def run_cmd(s1):
        print(s1); os.system(s1)

exit_now = False
# Read input from terminal
n_inp = sys.argv
if '--help' in n_inp:
    with open('manual.md','r') as fin:
        print fin.read()
    exit_now = True
    sys.exit()

if '--compare' in n_inp:
    ind1 = n_inp.index('--compare')
    compare_ans = True
    ref_file = n_inp[ind1+1]
    ref_file = os.path.abspath(ref_file)
    n_inp = n_inp[:ind1]+n_inp[ind1+2:]    

if '-p' in n_inp:
    ind1 = n_inp.index('-p')
    nJobs = int(n_inp[ind1+1])
    n_inp = n_inp[:ind1]+n_inp[ind1+2:] 
    run_parallel = True

if '--partition' in n_inp:
    ind1 = n_inp.index('--partition')
    partition_size = int(n_inp[ind1+1])
    n_inp = n_inp[:ind1]+n_inp[ind1+2:] 
    
if '-K' in n_inp:
    ind1 = n_inp.index('-K')
    K = int(n_inp[ind1+1])
    n_inp = n_inp[:ind1]+n_inp[ind1+2:]



if '-o' in n_inp:
    ind1 = n_inp.index('-o')
    comp_directory_name = n_inp[ind1+1]
    comp_directory_name=os.path.abspath(comp_directory_name)
    n_inp = n_inp[:ind1]+n_inp[ind1+2:]
    if os.listdir(comp_directory_name):
	print('ERROR: Output directory specified with -o needs to be an empty or non-existent directory')
	exit_now = True
else:
    print('ERROR: Output directory needed. Use -o flag, which is mandatory.')
    exit_now = True

reads_files = []

if '--left' in n_inp and '--right' in n_inp:
    ind1 = n_inp.index('--left')
    reads_files.append(os.path.abspath(n_inp[ind1+1]))
    n_inp = n_inp[:ind1]+n_inp[ind1+2:]

    ind2 = n_inp.index('--right')
    reads_files.append(os.path.abspath(n_inp[ind2+1]))
    n_inp = n_inp[:ind2]+n_inp[ind2+2:]

elif '--single' in n_inp:
    ind1 = n_inp.index('--single')
    reads_files.append(os.path.abspath(n_inp[ind1+1]))
    n_inp = n_inp[:ind1]+n_inp[ind1+2:]

else:
    print("ERROR: Need to specify single-ended reads with --single or specify both --left and --right for paired-ended reads.")
# if len(n_inp)>1:
#     comp_directory_name = n_inp[1]
#     reads_files = [n_inp[2]]
#     if len(n_inp)>3:
#         reads_files.append(n_inp[3])
# else:


    ''''with open('manual.md','r') as fin:
        print fin.read()'''

if exit_now:
    print('Try running python shannon.py --help for a short manual')   
    sys.exit()


if len(reads_files) == 1:
    paired_end = False
elif len(reads_files) == 2:
    paired_end = True

run_quorum = False
if reads_files[0][-1] == 'q': #Fastq mode
    run_quorum = True
    
paired_end_flag = ""
if paired_end:
    paired_end_flag = " --paired_end "
    
    
jellyfish_dir = ' jellyfish'
jellyfish_kmer_cutoff = 1

# For extension correction
sample_name = comp_directory_name.split('/')[-1] + "_"
new_directory_name = comp_directory_name
sample_name_input = comp_directory_name + "/" + sample_name
comp_size_threshold = partition_size
run_og_components = False
get_og_comp_kmers = False

# For kmers_for_component
kmer_directory = sample_name_input+"algo_input" # "./S15_SE_algo_input" # K = 24
base_directory_name = comp_directory_name #"./S15_SE_contig_partition"
r1_contig_file_extension = "contigs.txt"
r1_new_kmer_tag = "r1"
r1_graph_file_extension = ".txt"
randomize = False 

if use_partitioning == False:
    partition_size = 100000000000000000000000000000000000000000000000000000000000000000
    comp_size_threshold = partition_size

run_cmd('mkdir ' + comp_directory_name)
run_cmd('mkdir ' + sample_name_input+ "algo_input")


#Run Quorum now
if run_quorum:
    run_cmd('python run_quorum.py ' + comp_directory_name + ' ' + '\t'.join(reads_files))
    if paired_end:
	reads_files = [comp_directory_name + '/corrected_reads_1.fa',comp_directory_name + '/corrected_reads_2.fa']
    else:
	reads_files = [comp_directory_name + '/corrected_reads.fa']

reads_string = ' '.join(reads_files)    

#og_algo_input = './S15_SE_algo_input' #'./WingWongTest_K24_algo_input' #'./S15_SE_algo_input'
#og_algo_output = './S15_SE_algo_output' #'./WingWongTest_K24_algo_output' #'./S15_SE_algo_output'
#og_algo_output = './SPombe_algo_input'
#ref_file = './SPombe_algo_input/reference.fasta'


#og_algo_input = '~/Full_Assembler/WW_NM_005427_algo_input'
#og_algo_output = '~/Full_Assembler/WW_NM_005427_algo_output'
#og_algo_input = './WingWongTest_K24_algo_input' #'./S15_SE_algo_input'
#og_algo_output = './WingWongTest_K24_algo_output' #'./S15_SE_algo_output'
#og_algo_input = './S15_SE_algo_input'
#og_algo_output = './S15_SE_algo_output'
#og_algo_input = "./Snyder/Full"
#og_algo_output = "./Snyder/Full"

# Creates output directory
#run_cmd('mkdir ' + comp_directory_name)
#run_cmd('mkdir ' + sample_name_input+ "algo_input")

# Runs Jellyfish
if run_jellyfish:
    K_value = K
    run_jfs = ' '
    if double_stranded:
        run_jfs += ' -C '
    run_cmd('rm '+sample_name_input+'algo_input/jelly*')  #Remove old jellyfish files
    run_cmd(jellyfish_dir+' count -m ' + str(K_value+1) + run_jfs+ ' -o ' + sample_name_input+'algo_input/jellyfish_p1_output.jf -s 20000000 -c 4 -t 32 ' +reads_string)

    '''if os.path.isfile(sample_name_input+'algo_input/jellyfish_p1_output_1'):
        run_cmd(jellyfish_dir+' merge -o ' + sample_name_input+'algo_input/jellyfish_p1_output.jf ' + sample_name_input+'algo_input/jellyfish_p1_output\_*')
    else:
        run_cmd('mv ' + sample_name_input+'algo_input/jellyfish_p1_output_0 ' +sample_name_input+'algo_input/jellyfish_p1_output.jf')'''

    run_cmd(jellyfish_dir+' dump -c -t -L ' + str(jellyfish_kmer_cutoff) + ' ' + sample_name_input+'algo_input/jellyfish_p1_output.jf > ' +sample_name_input+'algo_input/k1mer.dict_org')
    if (not run_extension_corr) and double_stranded:
        tester.double_strandify(sample_name_input+'algo_input/k1mer.dict_org', sample_name_input+'algo_input/k1mer.dict')
    if (not run_extension_corr) and (not double_stranded):
        run_cmd('mv ' + sample_name_input+'algo_input/k1mer.dict_org ' + sample_name_input+'algo_input/k1mer.dict')

# Runs error correction for k1mers (Deletes error k1mers) using contig approach
# and determines seperate groups of contigs that share no kmers (components)
if run_extension_corr:
    run_cmd('rm -rf ' + base_directory_name+"/component*contigs.txt")    

    if double_stranded:
        str_ec = ' -d '
    else: 
        str_ec = ' '
    #run_cmd('python extension_correction_SIC.py ' + str_ec + base_dir+sample_name+'algo_input/k1mer.dict_org ' +base_dir+sample_name+'algo_input/k1mer.dict 3 75 12 0.5')
    run_cmd('python extension_correction.py ' + str_ec + sample_name_input+'algo_input/k1mer.dict_org ' +sample_name_input+'algo_input/k1mer.dict ' + str(hyp_min_weight) + ' ' + str(hyp_min_length) + ' ' + comp_directory_name + " " + str(comp_size_threshold))

# Gets kmers from k1mers
if run_jellyfish or run_extension_corr:
    run_cmd('python kp1mer_to_kmer.py ' + sample_name_input+'algo_input/k1mer.dict ' + sample_name_input+'algo_input/kmer.dict')


# Runs gpmetis to partition components of size above "partition_size" into partitions of size "partition_size"
# Gets k1mers, kmers, and reads for each partition
[components_broken, new_components] = kmers_for_component(kmer_directory, reads_files, base_directory_name, r1_contig_file_extension, r1_new_kmer_tag, r1_graph_file_extension, get_og_comp_kmers, get_partition_kmers, double_stranded, paired_end, False, partition_size, overload, K)

# This counts remaining and non-remaining partitions for log.
num_remaining = 0
num_non_remaining = 0 
for part in new_components:
    if "remaining" in part:
        num_remaining += 1
    else:
        num_non_remaining += 1

# If "use_second_partition", rerun gpmetis with a penalization for contig edges broken in old partitioning
# This to give a new partitioning for each component of size above "partition_size"
# Gets k1mers, kmers, and reads for each partition
if use_second_iteration:
    r2_graph_file_extension = "r2.txt"
    r2_new_kmer_tag = "r2"
    r2_contig_file_extension = "r2.txt"

    for i in components_broken:
        partition_file = "/component" + str(i+1) + r1_graph_file_extension + ".part." + str(components_broken[i])
        og_graph_file = "/component" + str(i+1) + r1_graph_file_extension
        new_graph_file = "/component" + str(i+1) + r2_graph_file_extension
        contig_file = "/component" + str(i+1) + r1_contig_file_extension
        new_contig_file = "/component" + str(i+1) + r2_contig_file_extension 
        weight_updated_graph(base_directory_name, partition_file, og_graph_file, new_graph_file, contig_file, new_contig_file, penalty, randomize)

    get_og_comp_kmers = 0
    get_partition_kmers = 1
    [r2_components_broken, r2_new_components] = kmers_for_component(kmer_directory, reads_files, base_directory_name, r1_contig_file_extension, r2_new_kmer_tag, r2_graph_file_extension, get_og_comp_kmers, get_partition_kmers, double_stranded, paired_end, True, partition_size, overload, K)

    # This counts remaining and non-remaining partitions for log.
    for part in r2_new_components:
        num_non_remaining += 1

# This code updates the log
if os.path.exists(comp_directory_name+"/before_sp_log.txt"):
    f_log = open(comp_directory_name+"/before_sp_log.txt", 'a')
else:
    f_log = open(comp_directory_name+"/before_sp_log.txt", 'w')
f_log.write(str(time.asctime()) + ": " +"Number of simple Partitions: " + str(num_remaining)  + "\n")
f_log.write(str(time.asctime()) + ": " +"Number of complex Partitions: " + str(num_non_remaining)  + "\n")
f_log.close()
        
# parameters for main_server call
main_server_parameter_string = ""
main_server_og_parameter_string = ""

# Create directories for each partition where main_server.py will be run
for comp in new_components:
    dir_base = comp_directory_name + "/" + sample_name + "_c" + str(comp)
    run_cmd("mkdir " + dir_base + "algo_input")
    run_cmd("mkdir " + dir_base + "algo_output")
    if paired_end:
        run_cmd("mv " + base_directory_name + "/reads_c" + str(comp) + "_1.fasta " + dir_base + "algo_input/reads_1.fasta")
        run_cmd("mv " + base_directory_name + "/reads_c" + str(comp) + "_2.fasta " + dir_base + "algo_input/reads_2.fasta")
    else:
        run_cmd("mv " + base_directory_name + "/reads_c" + str(comp) + ".fasta " + dir_base + "algo_input/reads.fasta")

    run_cmd("mv " + base_directory_name + "/component" + comp +  r1_new_kmer_tag + "kmers_allowed.dict " + dir_base + "algo_input/kmer.dict")
    run_cmd("mv " + base_directory_name + "/component" + comp +  r1_new_kmer_tag + "k1mers_allowed.dict " + dir_base + "algo_input/k1mer.dict")
    main_server_parameter_string = main_server_parameter_string + dir_base + " " 
    
# Create directories for each partition where main_server.py will be run for second partitioning 
if use_second_iteration:
    for comp in r2_new_components:
        dir_base = comp_directory_name + "/" + sample_name + "_r2_c" + str(comp)
        run_cmd("mkdir " + dir_base + "algo_input")
        run_cmd("mkdir " + dir_base + "algo_output")
        if paired_end:
            run_cmd("mv " + base_directory_name + "/reads_r2_c" + str(comp)+"_1.fasta " + dir_base + "algo_input/reads_1.fasta")
            run_cmd("mv " + base_directory_name + "/reads_r2_c" + str(comp)+"_2.fasta " + dir_base + "algo_input/reads_2.fasta")
        else:
            run_cmd("mv " + base_directory_name + "/reads_r2_c" + str(comp)+".fasta " + dir_base + "algo_input/reads.fasta")        
        
        run_cmd("mv " + base_directory_name + "/component" + comp +  r2_new_kmer_tag + "kmers_allowed.dict " + dir_base + "algo_input/kmer.dict")
        run_cmd("mv " + base_directory_name + "/component" + comp +  r2_new_kmer_tag + "k1mers_allowed.dict " + dir_base + "algo_input/k1mer.dict")
        main_server_parameter_string = main_server_parameter_string + dir_base + " " 


child_names = [x[0][:-10] for x in os.walk(comp_directory_name) if x[0].endswith('algo_input') and not x[0].endswith('_algo_input') and not x[0].endswith('allalgo_input')]
main_server_parameter_string = ' '.join(child_names)

        
# Run main_server.py for each partition in parallel
if double_stranded:
    ds_string = "  --ds "
else:
    ds_string = "  "

if run_parallel:
    run_cmd("parallel -j " + str(nJobs) + " python run_MB_SF.py {} --run_alg " + ds_string + " --kmer_size " + str(K)  + " " + paired_end_flag + " --dir_name " + comp_directory_name + " ::: " + main_server_parameter_string)
else:
    for param_str in main_server_parameter_string.split():
	run_cmd("python run_MB_SF.py " + param_str + " --run_alg " + ds_string + " --kmer_size " + str(K)  + " " + paired_end_flag + " --dir_name " + comp_directory_name + " " + param_str)


# locates all reconstructed files          
reconstructed_files = {}
for comp in new_components:
    dir_base = comp_directory_name + "/" + sample_name + "_c" + str(comp)
    dir = dir_base + "algo_output"
    if comp[0] in reconstructed_files:
        reconstructed_files[comp[0]].append(dir + '/' + 'reconstructed.fasta')
    else:
        reconstructed_files[comp[0]] = [dir + '/' + 'reconstructed.fasta']

if use_second_iteration:
    for comp in new_components:
        if "remaining" in comp:
            continue
        dir_base = comp_directory_name + "/" + sample_name + "_r2_c" + str(comp)
        dir = dir_base + "algo_output"
        if comp[0] in reconstructed_files:
            reconstructed_files[comp[0]].append(dir + '/' + 'reconstructed.fasta')
        else:
            reconstructed_files[comp[0]] = [dir + '/' + 'reconstructed.fasta']

# Creates new directory with concatenation of all reconstructed files
dir_base = comp_directory_name + "/" + sample_name + "_all"
dir_out = dir_base + "algo_output"
run_cmd("mkdir " + dir_out)
temp_file = dir_out + "/" + "all_reconstructed.fasta"
temp_file_args = "" 
for comp in reconstructed_files:    
    for file_name in reconstructed_files[comp]:
        temp_file_args = temp_file_args + file_name + " "
        
run_cmd("cat " + temp_file_args + " > " + temp_file)
process_concatenated_fasta(temp_file, dir_out + "/reconstructed.fasta")

# Compares reconstructed file against reference
if compare_ans:
	run_cmd("cp " + ref_file + ' ' +  dir_out + "/reference.fasta")
	run_cmd("python run_MB_SF.py " + dir_base + " --compare ")

# updates log
if os.path.exists(comp_directory_name+"/before_sp_log.txt"):
    f_log = open(comp_directory_name+"/before_sp_log.txt", 'a')
else:
    f_log = open(comp_directory_name+"/before_sp_log.txt", 'w')
num_transcripts = 0
with open(dir_out + "/" + "all_reconstructed.fasta", 'r') as reconstructed_transcripts:
    num_transcripts = len(reconstructed_transcripts.readlines())
f_log.write(str(time.asctime()) + ": " +"Program complete: " + str(num_transcripts) + " transcripts reconstructed" + "\n")
f_log.close()

# Creates final output
run_cmd('mkdir '+ comp_directory_name + '/TEMP')
run_cmd('mv ' + comp_directory_name + '/* ' + comp_directory_name + '/TEMP')
run_cmd('mv ' + comp_directory_name + '/TEMP/before_sp_log.txt ' + comp_directory_name + '/log.txt')
run_cmd('mv ' +  comp_directory_name + "/TEMP/" + sample_name + "_allalgo_output/reconstructed.fasta " + comp_directory_name + '/shannon.fasta')
run_cmd('more ' +  comp_directory_name + "/TEMP/*output.txt > " +comp_directory_name + '/terminal_output.txt') 

if compare_ans:
   run_cmd('mv ' +   comp_directory_name + "/TEMP/" + sample_name + "_allalgo_output/reconstr_log.txt "  + comp_directory_name + '/compare_log.txt')  
