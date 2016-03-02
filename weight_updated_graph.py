import time
import sys
import re
import pdb,math
import os
import os.path
import numpy as np

def weight_updated_graph(directory, partition_file, og_graph_file, new_graph_file, contig_file, new_contig_file, penalty = 5, randomize = False):
    ''' This function creates a gpmetis input graph where the edge weights
    between contigs is increased if the edge was broken in the previous partitioning.
    This makes it so paths that were broken in the previous partitioning
    are more likely to be present in the next partitioning since gpmetis avoids
    breaking edges with high weight according to it's optimization.
    '''
    f1 = open(directory + og_graph_file, 'r')
    f2 = open(directory + partition_file, 'r')
    f3 = open(directory + new_graph_file, 'w')

    lines1 = f1.readlines()
    lines2 = f2.readlines()


    if not randomize:
        f3.write(lines1[0])
        i = 0
        for line in lines1:
            if i != 0:
                contig = i
                new_line = ""
                tokens = line.split()
                j = 0
                for contig2 in tokens:
                    if j%2 == 0:
                        new_line = new_line + contig2 + "\t" 
                        if int(lines2[i-1]) != int(lines2[int(contig2)-1]):
                            new_line = new_line + str(penalty*int(tokens[j+1])) + "\t"
                        else:
                            new_line = new_line + tokens[j+1] + "\t" 
                    j += 1
                f3.write(new_line + '\n')
            i += 1
        
    if randomize:
        old_contig_f = open(directory+contig_file, 'r')
        new_contig_f = open(directory+new_contig_file, 'w')
        old_contig_lines = old_contig_f.readlines()
            
        f3.write(lines1[0])
        new_order = np.random.permutation(len(lines1)-1)
        for i in new_order:
            line = lines1[i+1]
            new_line = ''
            tokens = line.split()
            j = 0
            for contig2 in tokens:
                if j%2 == 0:
                    new_line = new_line + str(new_order[int(contig2)-1]+1) + "\t" 
                    if int(lines2[i]) != int(lines2[int(contig2)-1]):
                        new_line = new_line + tokens[j+1] + "\t"  # No penalty
                    else:
                        new_line = new_line + tokens[j+1] + "\t" 
                j += 1   
            f3.write(new_line + '\n')
        
        for i in new_order:
            new_contig_f.write(old_contig_lines[i])

            
