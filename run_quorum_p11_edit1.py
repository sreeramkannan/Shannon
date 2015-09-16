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

reads_files = []
paired_end = False
n_inp = sys.argv
if len(n_inp)>2:
    base_dir = sys.argv[1]
    reads_files = [sys.argv[2]]
    if len(n_inp)>3:
        reads_files.append(sys.argv[3])

#pdb.set_trace()        
        
if len(reads_files) == 2:
    paired_end = True
    
def run_cmd(s1):
	print(s1); os.system(s1)
    
def find_common_string(s1, s2):
    s3 = ""
    i = 0
    for each in s1: 
        if each == s2[i]:
            s3 = s3 + each
        i += 1
    return s3

if paired_end == False:
    # Run Quorum Normally
    base_file = base_dir + "/corrected_reads.fasta"
    run_cmd("/home/sreeramkannan/Packages/quorum-1.0.0/quorum --prefix " + base_file + " " + reads_files[0])
    
elif paired_end == True:  
    reads1_file = reads_files[0]
    reads2_file = reads_files[1]
    with open(reads1_file, 'r') as r1_file:
        with open(reads2_file, 'r') as r2_file:
        
            new_reads1_file = base_dir + "/new_reads1.fastq"
            with open(new_reads1_file, 'w') as new_r1_file:
                r1_reads = []
                j = 0
                i = 0
                old_name = ""
                new_name = ""
                for line in r1_file:
                    if i % 4 == 0:
                        #pdb.set_trace()
                        #new_name = find_common_string(r1_reads[j], line.split()[0])
                        new_name = line.split()[0] + "_1"
                        r1_reads.append(line.split()[0])
                        old_name = line.split()[0]
                        new_r1_file.write(line.replace(old_name, new_name))
                        #new_r2_file.write(line[i+1].replace(old_name, new_name))
                        #new_r2_file.write(line[i+2].replace(old_name, new_name))
                        #new_r2_file.write(line[i+3].replace(old_name, new_name))
                        j += 1
                    else:
                        new_r1_file.write(line.replace(old_name, new_name))
                    i += 1


            new_reads2_file = base_dir + "/new_reads2.fastq"
            with open(new_reads2_file, 'w') as new_r2_file:
                j = 0
                i = 0
                old_name = ""
                new_name = ""
                for line in r2_file:
                    if i % 4 == 0:
                        #pdb.set_trace()
                        #new_name = find_common_string(r1_reads[j], line.split()[0])
                        new_name = r1_reads[j] + "_2"
                        old_name = line.split()[0]
                        new_r2_file.write(line.replace(old_name, new_name))
                        #new_r2_file.write(line[i+1].replace(old_name, new_name))
                        #new_r2_file.write(line[i+2].replace(old_name, new_name))
                        #new_r2_file.write(line[i+3].replace(old_name, new_name))
                        j += 1
                    else:
                        new_r2_file.write(line.replace(old_name, new_name))
                    i += 1
    pdb.set_trace()
    base_file = base_dir + "/quorum_output"    
    run_cmd("/home/sreeramkannan/Packages/quorum-1.0.0/quorum  --prefix " + base_file + " " +  new_reads1_file + " " + new_reads2_file)
    
    with open(base_file + ".fa", 'r') as quorum_output: 
        with open(base_dir + "/corrected_reads_1.fasta", 'w') as cr1:
            with open(base_dir + "/corrected_reads_2.fasta", 'w') as cr2:
                left_ends = {}
                right_ends = {}
                i = 0
                write_next_cr1 = False
                write_next_cr2 = False
                prev_name = None
                add_next_to_left_ends = False
                add_next_to_right_ends = False
                #pdb.set_trace()
                for line in quorum_output:
                    if i % 2 == 0:
                        #pdb.set_trace()
                        if line.split()[0][-2:] == "_1":
                            if line.split()[0][:-2] + "_2" in right_ends:
                                cr1.write(line)
                                write_next_cr1 = True
                                cr2.write(right_ends[line.split()[0][:-2] + "_2"][0])
                                cr2.write(right_ends[line.split()[0][:-2] + "_2"][1])
                                right_ends.pop(line.split()[0][:-2] + "_2")
                            elif line.split()[0][:-2] + "_2" not in right_ends:
                                left_ends[line.split()[0]] = [line]
                                prev_name = line.split()[0]
                                add_next_to_left_ends = True
                        elif line.split()[0][-2:] == "_2":
                            #pdb.set_trace()
                            if line.split()[0][:-2] + "_1" in left_ends:
                                cr1.write(left_ends[line.split()[0][:-2] + "_1"][0])
                                cr1.write(left_ends[line.split()[0][:-2] + "_1"][1])
                                cr2.write(line)
                                #cr2.write(line[i+1])
                                write_next_cr2 = True
                                left_ends.pop(line.split()[0][:-2] + "_1")
                            elif line.split()[0][:-2] + "_1" not in left_ends:
                                right_ends[line.split()[0]] = [line]
                                prev_name = line.split()[0]
                                add_next_to_right_ends = True
                    else:
                        if write_next_cr1:
                            cr1.write(line)
                            write_next_cr1 = False

                        if write_next_cr2:
                            cr2.write(line)
                            write_next_cr2 = False   

                        if add_next_to_right_ends:
                            right_ends[prev_name].append(line)
                            add_next_to_right_ends = False
                        if add_next_to_left_ends:
                            left_ends[prev_name].append(line)   
                            add_next_to_left_ends = False
                    i += 1
                #pdb.set_trace()
                with open(base_dir + "/unpaired_reads.fasta", 'w') as unpaired_reads:
                    for each in left_ends:
                        unpaired_reads.write(left_ends[each][0])
                        unpaired_reads.write(left_ends[each][1])
                    
                    for each in right_ends:
                        unpaired_reads.write(right_ends[each][0])
                        unpaired_reads.write(right_ends[each][1])

            
            
            
    