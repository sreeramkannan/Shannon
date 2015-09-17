import time
import sys
import re
import pdb,math
import os
import os.path
import numpy as np
import tester
from set_exp import set_exp 
from filter_trans import filter_trans

from kmers_for_component_p11_edit1 import kmers_for_component
from weight_updated_graph_f import weight_updated_graph
from process_concatenated_fasta_f import process_concatenated_fasta
    
# For jellyfish
double_stranded = True
run_jellyfish = True
paired_end = False # Automatically set if command line is used

# General, Can be read in from terminal
reads_files = ['~/Full_Assembler/SPombe_algo_input/reads.fasta'] # ['./S15_SE_algo_input/reads.fasta']
sample_name = "SPombe_test"

# For extension correction
run_extension_corr = True
hyp_min_weight = 3
hyp_min_length = 75
partition_size = 5000
use_partitioning = True

# For kmers_for_component
use_second_iteration = True
get_partition_kmers = True
overload = 2
K = 24
penalty = 5

run_parallel = False
compare_ans = False
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def run_cmd(s1):
        print(s1); os.system(s1)


# Read input from terminal
n_inp = sys.argv
if '--compare' in n_inp:
    ind1 = n_inp.index('--compare')
    compare_ans = True
    ref_file = n_inp[ind1+1]
    n_inp = n_inp[:ind1]+n_inp[ind1+2:]    

if '-p' in n_inp:
    ind1 = n_inp.index('-p')
    nJobs = int(n_inp[ind1+1])
    n_inp = n_inp[:ind1]+n_inp[ind1+2:] 
    run_parallel = True
     
if len(n_inp)>1:
    comp_directory_name = n_inp[1]
    reads_files = [n_inp[2]]
    if len(n_inp)>3:
        reads_files.append(n_inp[3])
else:
    with open('manual.md','r') as fin:
        print fin.read()   
    sys.exit()


if len(reads_files) == 1:
    paried_end = False
elif len(reads_files) == 2:
    paired_end = True
    
paired_end_flag = ""
if paired_end:
    paired_end_flag = " --paired_end "
    
if paired_end:
	reads_string = reads_files[0] + ' ' + reads_files[1]
else:	
	reads_string = reads_files[0]

    
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

    
run_cmd('mkdir ' + comp_directory_name)
run_cmd('mkdir ' + sample_name_input+ "algo_input")

if run_jellyfish:
    K_value = K
    run_jfs = ' '
    if double_stranded:
        run_jfs += ' -C '
    run_cmd('rm '+sample_name_input+'algo_input/jelly*')  #Remove old jellyfish files
    #run_cmd(jellyfish_dir+' count -m ' + str(K_value)+ run_jfs+ ' -o ' + base_dir+sample_name+'algo_input/jellyfish_output -s 20000000 -c 4 -t 32 ' +reads_string)
    run_cmd(jellyfish_dir+' count -m ' + str(K_value+1) + run_jfs+ ' -o ' + sample_name_input+'algo_input/jellyfish_p1_output.jf -s 20000000 -c 4 -t 32 ' +reads_string)

    '''if os.path.isfile(sample_name_input+'algo_input/jellyfish_p1_output_1'):
        run_cmd(jellyfish_dir+' merge -o ' + sample_name_input+'algo_input/jellyfish_p1_output.jf ' + sample_name_input+'algo_input/jellyfish_p1_output\_*')
    else:
        run_cmd('mv ' + sample_name_input+'algo_input/jellyfish_p1_output_0 ' +sample_name_input+'algo_input/jellyfish_p1_output.jf')'''

    #run_cmd(jellyfish_dir+' dump -c -L ' + str(jellyfish_kmer_cutoff) + ' ' + base_dir+sample_name+'algo_input/jellyfish_output.jf > ' +base_dir+sample_name+'algo_input/kmer.dict')
    run_cmd(jellyfish_dir+' dump -c -t -L ' + str(jellyfish_kmer_cutoff) + ' ' + sample_name_input+'algo_input/jellyfish_p1_output.jf > ' +sample_name_input+'algo_input/k1mer.dict_org')
    if (not run_extension_corr) and double_stranded:
        tester.double_strandify(sample_name_input+'algo_input/k1mer.dict_org', sample_name_input+'algo_input/k1mer.dict')
    if (not run_extension_corr) and (not double_stranded):
        run_cmd('mv ' + sample_name_input+'algo_input/k1mer.dict_org ' + sample_name_input+'algo_input/k1mer.dict')

        
if run_extension_corr:
    run_cmd('rm -rf ' + base_directory_name+"/component*contigs.txt")    

    if double_stranded:
        str_ec = ' -d '
    else: 
        str_ec = ' '
    #run_cmd('python extension_correction_SIC.py ' + str_ec + base_dir+sample_name+'algo_input/k1mer.dict_org ' +base_dir+sample_name+'algo_input/k1mer.dict 3 75 12 0.5')
    run_cmd('python extension_correction_SIC_p11.py ' + str_ec + sample_name_input+'algo_input/k1mer.dict_org ' +sample_name_input+'algo_input/k1mer.dict ' + str(hyp_min_weight) + ' ' + str(hyp_min_length) + ' ' + comp_directory_name + " " + str(comp_size_threshold))

if run_jellyfish or run_extension_corr:
    run_cmd('python kp1mer_to_kmer.py ' + sample_name_input+'algo_input/k1mer.dict ' + sample_name_input+'algo_input/kmer.dict')

'''if run_remaining_comps:
    run_cmd('mkdir ' + sample_name + 'remaining_algo_input')
    run_cmd('cp ' + comp_directory_name+"/remaining_k1mers.txt " + sample_name+'remaining_algo_input/k1mer.dict')
    run_cmd('python kp1mer_to_kmer.py ' + base_dir+sample_name+'remaining_algo_input/k1mer.dict ' + base_dir+sample_name+'remaining_algo_input/kmer.dict')
    run_cmd('cp ' + reads_string + ' ' + sample_name+'remaining_algo_input/')'''

[components_broken, new_components] = kmers_for_component(kmer_directory, reads_files, base_directory_name, r1_contig_file_extension, r1_new_kmer_tag, r1_graph_file_extension, get_og_comp_kmers, get_partition_kmers, double_stranded, paired_end, False, partition_size, overload, K)
#pdb.set_trace()

num_remaining = 0
num_non_remaining = 0 
for part in new_components:
    if "remaining" in part:
        num_remaining += 1
    else:
        num_non_remaining += 1

if use_second_iteration:
    r2_graph_file_extension = "r2.txt"
    r2_new_kmer_tag = "r2"
    r2_contig_file_extension = "r2.txt"

    for i in components_broken:
        #pdb.set_trace()
        partition_file = "/component" + str(i+1) + r1_graph_file_extension + ".part." + str(components_broken[i])
        og_graph_file = "/component" + str(i+1) + r1_graph_file_extension
        new_graph_file = "/component" + str(i+1) + r2_graph_file_extension
        contig_file = "/component" + str(i+1) + r1_contig_file_extension
        new_contig_file = "/component" + str(i+1) + r2_contig_file_extension 
        weight_updated_graph(base_directory_name, partition_file, og_graph_file, new_graph_file, contig_file, new_contig_file, penalty, randomize)

    #pdb.set_trace()
    get_og_comp_kmers = 0
    get_partition_kmers = 1
    [r2_components_broken, r2_new_components] = kmers_for_component(kmer_directory, reads_files, base_directory_name, r1_contig_file_extension, r2_new_kmer_tag, r2_graph_file_extension, get_og_comp_kmers, get_partition_kmers, double_stranded, paired_end, True, partition_size, overload, K)

    for part in r2_new_components:
        num_non_remaining += 1

if os.path.exists(comp_directory_name+"/before_sp_log.txt"):
    f_log = open(comp_directory_name+"/before_sp_log.txt", 'a')
else:
    f_log = open(comp_directory_name+"/before_sp_log.txt", 'w')

f_log.write(str(time.asctime()) + ": " +"Number of simple Partitions: " + str(num_remaining)  + "\n")
f_log.write(str(time.asctime()) + ": " +"Number of complex Partitions: " + str(num_non_remaining)  + "\n")

f_log.close()
        
       
#new_components = ["1_0", "1_1"]
#components_broken = ["1"]
main_server_parameter_string = ""
main_server_og_parameter_string = ""


for comp in new_components:
    dir_base = comp_directory_name + "/" + sample_name + "_c" + str(comp)
    run_cmd("mkdir " + dir_base + "algo_input")
    run_cmd("mkdir " + dir_base + "algo_output")
    #run_cmd("cp " + base_directory_name + "/reads_c" + str(comp) + ".fasta " + dir_base + "algo_input/reads.fasta")
    if paired_end:
        run_cmd("cp " + base_directory_name + "/reads_c" + str(comp) + "_1.fasta " + dir_base + "algo_input/reads_1.fasta")
        run_cmd("cp " + base_directory_name + "/reads_c" + str(comp) + "_2.fasta " + dir_base + "algo_input/reads_2.fasta")
    else:
        run_cmd("cp " + base_directory_name + "/reads_c" + str(comp) + ".fasta " + dir_base + "algo_input/reads.fasta")

    run_cmd("cp " + base_directory_name + "/component" + comp +  r1_new_kmer_tag + "kmers_allowed.dict " + dir_base + "algo_input/kmer.dict")
    run_cmd("cp " + base_directory_name + "/component" + comp +  r1_new_kmer_tag + "k1mers_allowed.dict " + dir_base + "algo_input/k1mer.dict")
    #run_cmd("cp " + og_algo_output + "/reference.fasta " + dir_base + "algo_output")
    main_server_parameter_string = main_server_parameter_string + dir_base + " " 
    #run_cmd("python main_server_WWY_Kayvon.py " + dir_base + " run_alg")
    
if use_second_iteration:
    for comp in r2_new_components:
        dir_base = comp_directory_name + "/" + sample_name + "_r2_c" + str(comp)
        run_cmd("mkdir " + dir_base + "algo_input")
        run_cmd("mkdir " + dir_base + "algo_output")
        #run_cmd("cp " + base_directory_name + "/reads_r2_c" + str(comp)+".fasta " + dir_base + "algo_input/reads.fasta")
        if paired_end:
            run_cmd("cp " + base_directory_name + "/reads_r2_c" + str(comp)+"_1.fasta " + dir_base + "algo_input/reads_1.fasta")
            run_cmd("cp " + base_directory_name + "/reads_r2_c" + str(comp)+"_2.fasta " + dir_base + "algo_input/reads_2.fasta")
        else:
            run_cmd("cp " + base_directory_name + "/reads_r2_c" + str(comp)+".fasta " + dir_base + "algo_input/reads.fasta")        
        
        run_cmd("cp " + base_directory_name + "/component" + comp +  r2_new_kmer_tag + "kmers_allowed.dict " + dir_base + "algo_input/kmer.dict")
        run_cmd("cp " + base_directory_name + "/component" + comp +  r2_new_kmer_tag + "k1mers_allowed.dict " + dir_base + "algo_input/k1mer.dict")
        #run_cmd("cp " + og_algo_output + "/reference.fasta " + dir_base + "algo_output")
        main_server_parameter_string = main_server_parameter_string + dir_base + " " 
        #run_cmd("python main_server_WWY_Kayvon.py " + dir_base + " run_alg")

'''if run_remaining_comps:
    main_server_parameter_string = main_server_parameter_string + sample_name+'remaining_'     '''
        
#pdb.set_trace()
        
# stdbuf -oL  python filename | tee output_log.txt       
if double_stranded:
    ds_string = "  --ds "
else:
    ds_string = "  "

if run_parallel:
    run_cmd("parallel -j " + str(nJobs) + " python main_server_WWY_p11_edit1.py {} --run_alg " + ds_string + " --kmer_size " + str(K)  + " " + paired_end_flag + " --dir_name " + comp_directory_name + " ::: " + main_server_parameter_string)
else:
    for param_str in main_server_parameter_string.split():
	run_cmd("python main_server_WWY_p11_edit1.py " + param_str + " --run_alg " + ds_string + " --kmer_size " + str(K)  + " " + paired_end_flag + " --dir_name " + comp_directory_name + " " + param_str)


#run_cmd("parallel stdbuf -oL python main_server_WWY_Kayvon.py {} run_alg | tee " + { + "_terminal_output.txt ::: " + main_server_parameter_string)
        
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
        dir_base = comp_directory_name + "/" + sample_name + "_r2_c" + str(comp)
        dir = dir_base + "algo_output"
        if comp[0] in reconstructed_files:
            reconstructed_files[comp[0]].append(dir + '/' + 'reconstructed.fasta')
        else:
            reconstructed_files[comp[0]] = [dir + '/' + 'reconstructed.fasta']

'''
for comp in reconstructed_files:
    dir_base = new_directory_name + "_concat_c" + str(comp)
    dir = dir_base + "algo_output"
    run_cmd("mkdir " + dir)
    temp_file = dir + "/" + comp + "_reconstructed.fasta"
    temp_file_args = ""
    for file in reconstructed_files[comp]:
        temp_file_args = temp_file_args + file + " " 
    temp_file_args = temp_file_args + sample_name+'remaining_algo_output/reconstructed.fasta'
    run_cmd("cat " + temp_file_args + " > " + temp_file)
    process_concatenated_fasta(temp_file, dir + "/reconstructed.fasta")
    run_cmd("cp " + og_algo_output + "/reference.fasta " + dir + "/")

for comp in reconstructed_files:
    dir_base = new_directory_name + "_concat_c" + str(comp)
    run_cmd("python main_server_WWY_Kayvon.py " + dir_base + " compare ")'''

dir_base = comp_directory_name + "/" + sample_name + "_all"
dir_out = dir_base + "algo_output"
run_cmd("mkdir " + dir_out)
temp_file = dir_out + "/" + "all_reconstructed.fasta"
temp_file_args = "" 
for comp in reconstructed_files:    
    for file_name in reconstructed_files[comp]:
        temp_file_args = temp_file_args + file_name + " "
        
#temp_file_args = temp_file_args #+ sample_name #+'remaining_algo_output/reconstructed.fasta'
run_cmd("cat " + temp_file_args + " > " + temp_file)
process_concatenated_fasta(temp_file, dir_out + "/reconstructed.fasta")

if compare_ans:
	run_cmd("cp " + ref_file + ' ' +  dir_out + "/")
	run_cmd("python main_server_WWY_p11_edit1.py " + dir_base + " --compare ")


if os.path.exists(comp_directory_name+"/before_sp_log.txt"):
    f_log = open(comp_directory_name+"/before_sp_log.txt", 'a')
else:
    f_log = open(comp_directory_name+"/before_sp_log.txt", 'w')
    
num_transcripts = 0
with open(dir_out + "/" + "all_reconstructed.fasta", 'r') as reconstructed_transcripts:
    num_transcripts = len(reconstructed_transcripts.readlines())

f_log.write(str(time.asctime()) + ": " +"Program complete: " + str(num_transcripts) + " transcripts reconstructed" + "\n")

f_log.close()

run_cmd('mkdir '+ comp_directory_name + '/TEMP')
run_cmd('mv ' + comp_directory_name + '/* ' + comp_directory_name + '/TEMP')
run_cmd('mv ' + comp_directory_name + '/TEMP/before_sp_log.txt ' + comp_directory_name + '/log.txt')
run_cmd('mv ' +  comp_directory_name + "/TEMP/" + sample_name + "_allalgo_output/reconstructed.fasta " + comp_directory_name + '/shannon.fasta')
run_cmd('more ' +  comp_directory_name + "/TEMP/*output.txt > " +comp_directory_name + '/terminal_output.txt') 

if compare_ans:
   run_cmd('mv ' +   comp_directory_name + "/TEMP/" + sample_name + "_allalgo_output/reconstr_log.txt "  + comp_directory_name + '/compare_log.txt')  
 
#run_cmd("mv *"+sample_name+"* "+comp_directory_name+"/")
#run_cmd("mv "+comp_directory_name+"/"+sample_name+"algo_input ./"+sample_name+"algo_input")
# mv all directories to one directory at the end.
# Need to concatenate all reconstructed files (including remaining) to compare them    
# Make contig graph generation seperate from error correction.  
# What is the difference in SIC_correction that Sreeram was talking about?
# When component is size 0 do something so it doesn't crash 
# Seperate reads for each component
 
    
