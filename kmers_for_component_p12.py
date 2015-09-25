import time
import sys
import re
import pdb,math
import os
import os.path
from collections import Counter

def run_cmd(s1):
    print(s1); os.system(s1) 

class LocalCounter():
    '''Used for printing the number of kmers processed'''
    def __init__(self, name, report_length):
        self.name = name
        self.count = 0
        self.report_length = report_length

    def increment(self):
        self.count += 1
        if self.count % self.report_length == 0:
            print "{:s}: {:s}, processed {:d}".format(time.asctime(), self.name, self.count)

def reverse_complement(bases):
    """Return the reverse complement of BASES. Assumes BASES is
    all uppercase.

    >>> Read.reverse_complement("ATCGGGG")
    'CCCCGAT'
    """
    replacements = [('A', 't'), ('T', 'a'), ('C', 'g'), ('G', 'c')]
    for ch1, ch2 in replacements:
        bases = re.sub(ch1, ch2, bases)
    return bases[::-1].upper()

def kmers_for_component(kmer_directory, reads_files, directory_name, contig_file_extension, new_kmer_tag, graph_file_extension, get_og_comp_kmers, get_partition_kmers, double_stranded, paired_end = False, second_iteration = False,  partition_size = 500, overload = 1.5, K = 24):
    """This fuction runs gpmetis on the components above a threshold size.  It then creates a dictionary called 
    k1mers2component {k1mer : component}.  It then sends any reads that share a kmer with a component to that component.
    It then cycles through the k1mer file and outputs the k1mers along with their weights to a file for each component.
    It then creates a kmers2component dictionary.  It then outputs a kmers file for each component. 
    """

    
    if os.path.exists(directory_name+"/before_sp_log.txt"):
        f_log = open(directory_name+"/before_sp_log.txt", 'a')
    else:
        f_log = open(directory_name+"/before_sp_log.txt", 'w')
        
    # Counts number of components above size threshold
    Num_Components = 0
    i = 1
    while os.path.exists(directory_name+"/component" + str(i) + "contigs.txt"):
        i += 1
        Num_Components += 1  

    # Counts number of components below size threshold
    Num_Remaining_Components = 0
    i = 1
    while os.path.exists(directory_name+"/remaining_contigs" + str(i) + ".txt"):
        i += 1
        Num_Remaining_Components += 1   
       
    if get_partition_kmers:
        ufactor = int(1000.0*overload - 1000.0)
        components_broken = {}
        temp_string = ""
        
        # Runs gpmetis
        for i in range(Num_Components):
            with open(directory_name+"/component" + str(i+1) + contig_file_extension , 'r') as f:
                lines = f.readlines()
                num_contigs = len(lines)
                Partitions = int(math.ceil(float(num_contigs)/float(partition_size)))
                temp_string += "Component " + str(i) + ": " + str(Partitions) + " partitions, "  
                if len(lines) >= 2:
                    if ufactor == 0:
                        print("gpmetis " +directory_name+"/component" + str(i+1) + graph_file_extension + " " + str(Partitions))
                        run_cmd("gpmetis " +directory_name+"/component" + str(i+1) +  graph_file_extension + " " + str(Partitions) )
                        components_broken[i] = Partitions
                    elif ufactor > 0:
                        print("gpmetis " + " -ufactor="+str(ufactor) + " "+directory_name+"/component" + str(i+1) + graph_file_extension + " " + str(Partitions))
                        run_cmd("gpmetis " + " -ufactor="+str(ufactor) + " "+directory_name+"/component" + str(i+1) + graph_file_extension + " " +str(Partitions))
                        components_broken[i] = Partitions          
        
        f_log.write(str(time.asctime()) + ": " + "gpmetis for partitioning " + str(int(second_iteration)+1) + " is complete -- " + temp_string + "\n")
        
        new_components = {}
        k1mers2component = {}
        
        # Builds k1mer2component dictionary
        for i in components_broken:
            with open(directory_name+"/component" + str(i+1) + contig_file_extension, 'r') as f_contigs:
                with open(directory_name+"/component" + str(i+1) + graph_file_extension + ".part." + str(components_broken[i]) , 'r') as f_component:
                    contig_lines = f_contigs.readlines()
                    j = 0
                    for line in f_component:
                        tokens = line.split()
                        comp = str(i+1) + "_" + tokens[0]
                        contig = contig_lines[j].split()[0]
                        if comp not in new_components:
                            new_components[comp] = [contig]
                        else:
                            new_components[comp].append(contig)
                        for each in range(len(contig)-(K+1) + 1):
                            k1mers2component[contig[each:each+(K+1)]] = [comp, 0.0]
                        j += 1
                        
                         
        # Adds remaining components to k1mers2component
        if not second_iteration:
            for i in range(Num_Remaining_Components):
                with open(directory_name+"/remaining_contigs"+str(i+1)+".txt", 'r') as remaining_contigs_file:
                    if 1:
                        lines = remaining_contigs_file.readlines()
                        comp = "remaining"+str(i+1)
                        j = 0
                        for line in lines:
                            tokens = line.split()
                            contig = tokens[0]
                            if comp not in new_components:
                                new_components[comp] = [contig]
                            else:
                                new_components[comp].append(contig)
                            for each in range(len(contig)-(K+1) + 1):
                                k1mers2component[contig[each:each+(K+1)]] = [comp, 0.0]
                            j += 1                             

        f_log.write(str(time.asctime()) + ": " + "k1mers2component dictionary created " + "\n")

                            

        iter_tag = "_c"
        if second_iteration:
            iter_tag = "_r2_c"
        
        if paired_end == False:    
            comp2reads = {}
            comp2reads_reversed = {}
            readfile = open(reads_files[0], 'r')
            rlines = readfile.readlines()
            j = 0
            for line in rlines:
                if line.split()[0][0] != ">":
                    current_comps = Counter()
                    read = line.split()[0]
                    i = 0
                    while i < len(read) - (K+1):
                        k1mer = read[i:i+(K+1)]
                        if k1mer in k1mers2component:
                            #pdb.set_trace()
                            comp = k1mers2component[k1mer][0]
                            for each in comp:
                                current_comps[each] +=1

                        i += K
                    k1mer = read[-K:]
                    if k1mer in k1mers2component:
                        #pdb.set_trace()
                        comp = k1mers2component[k1mer][0]
                        for each in comp:
                            current_comps[each] +=1
                    if current_comps:
                        comp = max(current_comps.iterkeys(), key=(lambda key: current_comps[key]))
                        if comp in comp2reads:
                            comp2reads[comp].add(j)
                        else:
                            comp2reads[comp]=set(); comp2reads[comp].add(j)



                    if double_stranded:
                        current_comps_reversed = Counter()
                        rc_read = reverse_complement(read)
                        i = 0
                        while i < len(rc_read) - (K+1):
                            k1mer = rc_read[i:i+(K+1)]
                            if k1mer in k1mers2component:
                                #pdb.set_trace()
                                comp = k1mers2component[k1mer][0]
                                for each in comp:
                                    current_comps_reversed[each] +=1
                            i += K
                            
                        k1mer = rc_read[-K:]
                        if k1mer in k1mers2component:
                            #pdb.set_trace()
                            comp = k1mers2component[k1mer][0]
                            for each in comp:
                                current_comps_reversed[each] +=1

                            
                    if current_comps_reversed:
                        comp = max(current_comps_reversed.iterkeys(), key=(lambda key: current_comps_reversed[key]))    
                        if comp in comp2reads_reversed:
                                comp2reads_reversed[comp].add(j)
                        else:
                                comp2reads_reversed[comp] = set(); comp2reads_reversed[comp].add(j)
                    
                j += 1

            for comp in new_components:
                with open(directory_name+"/reads"+iter_tag+str(comp)+".fasta", 'w') as reads_part_file:
                    if comp in comp2reads:
                        for read_ind in comp2reads[comp]:
                            reads_part_file.write(rlines[read_ind-1])
                            reads_part_file.write(rlines[read_ind])
                    if comp in comp2reads_reversed: 
                        for read_ind in comp2reads_reversed[comp]:
                            read_id = rlines[read_ind-1].split()[0]
                            read_id_rem = "\t".join(rlines[read_ind-1].split()[1:])
                            read = rlines[read_ind]
                            reads_part_file.write(read_id + "_reversed \t" + read_id_rem + "\n")  
                            reads_part_file.write(reverse_complement(read) + " \n")
        

        
        elif paired_end == True:
            comp2reads = {}
            comp2reads_reversed = {}
            readfile1 = open(reads_files[0], 'r')
            readfile2 = open(reads_files[1], 'r')
            r1lines = readfile1.readlines()
            r2lines = readfile2.readlines()
            j = 0
            for line1 in r1lines:
                read1 = line1.split()[0] ## shouldn't this be line.split()[1]
                line2 = r2lines[j]
                read2 = line2.split()[0]
                read1_reversed = reverse_complement(read1)
                read2_reversed = reverse_complement(read2)
                current_comps = Counter()               
                i = 0
                while i < len(read1) - (K+1):
                    k1mer = read1[i:i+(K+1)]
                    if k1mer in k1mers2component:
                        #pdb.set_trace()
                        comp = k1mers2component[k1mer][0];
                        for each in comp:
                            current_comps[each] +=1
                        
                    i += K
                    
                k1mer = read1[-K:]
                if k1mer in k1mers2component:
                    #pdb.set_trace()
                    comp = k1mers2component[k1mer][0]
                    for each in comp:
                        current_comps[each] +=1
                i = 0
                while i < len(read2_reversed) - (K+1):
                    k1mer = read2_reversed[i:i+(K+1)]
                    if k1mer in k1mers2component:
                        #pdb.set_trace()
                        comp = k1mers2component[k1mer][0]
                        for each in comp:
                            current_comps[each] +=1
                        
                    i += K
                    
                k1mer = read2_reversed[-K:]
                if k1mer in k1mers2component:
                    #pdb.set_trace()
                    comp =  
                    for each in comp:
                        current_comps[each] +=1
                    
                if current_comps:
                    comp = max(current_comps.iterkeys(), key=(lambda key: current_comps[key]))
                    if comp in comp2reads:
                        comp2reads[comp].add(j)
                    else:
                        comp2reads[comp]=set(); comp2reads[comp].add(j)
 

                    ## You will have to add both reads in the pair in the right direction if either read is added to the comp.
                    
                if double_stranded:
                    i = 0; current_comps_reversed = Counter()
                    while i < len(read1_reversed) - (K+1):
                        k1mer = read1_reversed[i:i+(K+1)]
                        if k1mer in k1mers2component:
                            #pdb.set_trace()
                            comp = k1mers2component[k1mer][0]
                            for each in comp:
                                current_comps_reversed[each] += 1
                            
                        i += K
                        
                    k1mer = read1_reversed[-K:]
                    if k1mer in k1mers2component:
                        #pdb.set_trace()
                        comp = k1mers2component[k1mer][0]; 
                        for each in comp:
                            current_comps_reversed[each] += 1
                    i = 0            
                    while i < len(read2) - (K+1):
                        k1mer = read2[i:i+(K+1)]
                        if k1mer in k1mers2component:
                            #pdb.set_trace()
                            comp = k1mers2component[k1mer][0]
                            for each in comp:
                                current_comps_reversed[each] += 1
                            #current_comps_reversed[comp] = current_comps_reversed.get(comp,0) + 1
                        i += K
                        
                    k1mer = read2[-K:]
                    if k1mer in k1mers2component:
                        #pdb.set_trace()
                        comp = k1mers2component[k1mer][0]
                        for each in comp:
                            current_comps_reversed[each] += 1

                
                    if current_comps_reversed:
                        comp = max(current_comps_reversed.iterkeys(), key=(lambda key: current_comps_reversed[key]))
                        if comp in comp2reads_reversed:
                                comp2reads_reversed[comp].add(j)
                        else:
                                comp2reads_reversed[comp] = set(); comp2reads_reversed[comp].add(j)

                
                    ## You will have to add both reads in the pair in the right direction if either read is added to the comp.
              
                
                j += 1
    
                
            for comp in new_components:
                with open(directory_name+"/reads"+iter_tag+str(comp)+"_1.fasta", 'w') as reads_part_file1:
                    with open(directory_name+"/reads"+iter_tag+str(comp)+"_2.fasta", 'w') as reads_part_file2:
                        if comp in comp2reads:
                            for read_ind in comp2reads[comp]:
                                reads_part_file1.write(r1lines[read_ind-1])
                                reads_part_file1.write(r1lines[read_ind])
                                reads_part_file2.write(r2lines[read_ind-1])
                                reads_part_file2.write(r2lines[read_ind])

                        if double_stranded:
                            if comp in comp2reads_reversed: 
                                for read_ind in comp2reads_reversed[comp]:
                                    # I flipped the order of the reads in pair for reverse complement
                                    id1 = r1lines[read_ind-1].split()[0]
                                    id1_rem = "\t".join(r1lines[read_ind-1].split()[1:])
                                    read1 = r1lines[read_ind]
                                    reads_part_file2.write(id1 + "_reversed\t" + id1_rem + "\n")  
                                    reads_part_file2.write(reverse_complement(read1) + " \n")

                                    id2 = r2lines[read_ind-1].split()[0]
                                    id2_rem = "\t".join(r2lines[read_ind-1].split()[1:])
                                    read2 = r2lines[read_ind]
                                    reads_part_file1.write(id2 + "_reversed \t" + id2_rem + "\n")  
                                    reads_part_file1.write(reverse_complement(read2) + " \n")
       
        f_log.write(str(time.asctime()) + ": " + "reads partititoned " + "\n")        
        
        c2 = LocalCounter("k1mer Streaming", 10**6)

        # Cycles through k1mer file and adds k1mer weights to k1mers2component dictionary
        with open(kmer_directory + "/k1mer.dict" , 'r') as f:
            for line in f:
                if len(line) == 0: continue
                c2.increment()
                k1mer, weight = line.split()
                k1mer = k1mer.upper()
                weight = (float(weight))
                if k1mer in k1mers2component:
                    k1mers2component[k1mer][1] += weight
                if double_stranded:
                    r_k1mer = reverse_complement(k1mer)
                    if r_k1mer in k1mers2component:
                        k1mers2component[r_k1mer][1] += weight

        f_log.write(str(time.asctime()) + ": " + "k1mers weights stored in dictionary" + "\n")

        
        # Writes out k1mers with weights for each partition
        for comp in new_components:
            with open(directory_name+"/component" + comp +  new_kmer_tag + "k1mers_allowed.dict" , 'w') as new_file:
                for contig in new_components[comp]:
                    for i in range(len(contig)-(K+1) + 1):
                        k1mer = contig[i:i+(K+1)]
                        new_file.write(k1mer + "\t" + str(k1mers2component[k1mer][1]) + "\n")

        f_log.write(str(time.asctime()) + ": " + "k1mers written to file " + "\n")
                       
                        
            
        del(k1mers2component)
        
        kmers2component = {}
        
        # builds kmers2components dictionary that will store kmer weights
        for i in components_broken:
            with open(directory_name+"/component" + str(i+1) + contig_file_extension, 'r') as f_contigs:
                with open(directory_name+"/component" + str(i+1) + graph_file_extension + ".part." + str(components_broken[i]) , 'r') as f_component:
                    contig_lines = f_contigs.readlines()
                    j = 0
                    for line in f_component:
                        tokens = line.split()
                        comp = str(i+1) + "_" + tokens[0]
                        contig = contig_lines[j].split()[0]
                        for each in range(len(contig)-(K) + 1):
                            kmers2component[contig[each:each+(K)]] = 0.0
                        j += 1
                        
        # adds kmers from remaining components if not second iterarion
        if not second_iteration:
            for i in range(Num_Remaining_Components):
                with open(directory_name+"/remaining_contigs"+str(i+1)+".txt", 'r') as remaining_contigs_file:
                    if 1:
                        lines = remaining_contigs_file.readlines()
                        comp = "remaining"+str(i+1)

                        j = 0
                        for line in lines:
                            tokens = line.split()
                            contig = tokens[0]
                            for each in range(len(contig)-(K) + 1):
                                kmers2component[contig[each:each+(K)]] = 0.0
                            j += 1                            

        f_log.write(str(time.asctime()) + ": " + "kmers2component dictionary created " + "\n")
                            
        c2 = LocalCounter("kmer Streaming", 10**6)

        # gets kmer weights from kmer file and adds them to kmers2component
        with open(kmer_directory + "/kmer.dict" , 'r') as f:
            for line in f:
                if len(line) == 0: continue
                c2.increment()
                kmer, weight = line.split()
                kmer = kmer.upper()
                weight = (float(weight))
                if kmer in kmers2component:
                    kmers2component[kmer] += weight
                if double_stranded:
                    r_kmer = reverse_complement(kmer)
                    if r_kmer in kmers2component:
                        kmers2component[r_kmer] += weight

        f_log.write(str(time.asctime()) + ": " + "kmers weights stored in dictionary" + "\n")
                        
        # Writes out kmers and weights for each partition          
        for comp in new_components:
            with open(directory_name+"/component" + comp +  new_kmer_tag + "kmers_allowed.dict" , 'w') as new_file:
                for contig in new_components[comp]:
                    for i in range(len(contig)-(K) + 1):
                        kmer = contig[i:i+(K)]
                        new_file.write(kmer + "\t" + str(kmers2component[kmer]) + "\n")

        f_log.write(str(time.asctime()) + ": " + "kmers written to file " + "\n")
                        
        return [components_broken, new_components]
