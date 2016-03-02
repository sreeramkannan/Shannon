import time
import sys
import re
import pdb,math
import numpy
from sets import Set

BASES = ['A', 'G', 'C', 'T']
r=K-1
rmer_to_contig = {}
contig_to_rmer = {}

cmer_to_contig = {}

def reverse_complement(bases):
    """Return the reverse complement of BASES. Assumes BASES is
    all uppercase.
    """
    replacements = [('A', 't'), ('T', 'a'), ('C', 'g'), ('G', 'c')]
    for ch1, ch2 in replacements:
        bases = re.sub(ch1, ch2, bases)
    return bases[::-1].upper()

class Counter():
    def __init__(self, name, report_length):
        self.name = name
        self.count = 0
        self.report_length = report_length

    def increment(self):
        self.count += 1
        if self.count % self.report_length == 0:
            print "{:s}: {:s}, processed {:d}".format(time.asctime(), self.name, self.count)


def find_kmers(contig,k,ds):
    '''find the kmers of a contig of size k with ds denoting double stranded)'''
    rmer_list = []
    for i in range(0, len(contig)-r+1):
        rmer_list.append(contig[i:i+r])
        if ds:
            rmer_list.append(reverse_complement(contig[i:i+r]))
    return rmer_list



def argmax(lst, key):
    """Returns the element x in LST that maximizes KEY(x).
    """
    best = lst[0]
    for x in lst[1:]:
        if key(x) > key(best):
            best = x
    return best


def duplicate_check_ends(contigs,contig_name, ds, f = 0.9):
        contig = contigs[contig_name];
        qName = contig_name; qLen = len(contig);
        first_rmer = contig[:r]; last_rmer = contig[-r:]
        if first_rmer in rmer_to_contig and last_rmer in rmer_to_contig:
                


def duplicate_check(contigs, contig_name, ds, f = 0.99):
    # To add: if rmer in the contig multiple times, only increment the dup-contig once for each time its in dup-contig
    dup_count = {}
    contig = contigs[contig_name]
    qName = contig_name; qLen = len(contig); 
    max_till_now = 0; max_contig_index = -1
    for rmer in find_kmers(contig,r,ds):
        if rmer in rmer_to_contig:
            for dup in rmer_to_contig[rmer]:
                if dup == contig_name:
                    'Only take duplicates other than curr contig'
                    continue
                if dup not in dup_count:
                    dup_count[dup] = 1
                else:
                    dup_count[dup] += 1
                if dup_count[dup] >= max_till_now:
                    max_till_now = dup_count[dup]; max_contig_name = dup

    imp_dups = {} #Duplicates that are considered
    for dup in dup_count:
        if dup_count[dup] >= f*(qLen-r):
            imp_dups[dup] = numpy.zeros(len(contig))

    imp_dups_set = Set(list(imp_dups.keys()))

    a = numpy.zeros(len(contig))
    for i in range(len(contig)-r+1):
        rmer = contig[i:i+r]
        if rmer in rmer_to_contig:
            for dup in imp_dups_set.intersection(rmer_to_contig[rmer]):
                imp_dups[dup][i:i+r] = 1
        if ds:
            rmer_rc=reverse_complement(rmer)
            if rmer_rc in rmer_to_contig:
                for dup in imp_dups_set.intersection(rmer_to_contig[rmer]):
                    imp_dups[dup][i:i+r] = 1
    for dup in imp_dups: #for dup in dup_count:
        tName = dup; tLen = len(contigs[tName])
        score = sum(imp_dups[dup])
        if score >= f*qLen:
            if (tLen > qLen) or ((tLen == qLen) and (tName > qName)):
                return True   
    return False
           
        
        
def find_reps(infile, outfile,ds):
    contigs = {}
    contig_no = 0 
    for line in open(infile):
        if line[0]=='>':
            temp = line.strip().split(); curr_name = temp[0][1:]
            contig_no += 1
            if contig_no % 1000 == 0:
                print ("processed " + str(contig_no) + " sequences")
            continue
        else:
            curr_contig = line.strip()
            contigs[curr_name] = curr_contig
            for (i,rmer) in enumerate(find_kmers(curr_contig,r,ds)):
                if rmer in rmer_to_contig:
                    rmer_to_contig[rmer].add([curr_name,i])
                else:
                    rmer_to_contig[rmer]= Set([curr_name,i])
    with open(outfile,'w') as out_file:
	contig_no = 0
        for contig_name in contigs:
	    contig_no +=1
            duplicate_suspect = duplicate_check(contigs,contig_name, ds)
	    if contig_no % 1000 == 0:
		print("written " + str(contig_no) + " sequences")
            if not duplicate_suspect:
                out_file.write('>' + contig_name+'\n')
                out_file.write(contigs[contig_name]+'\n')
             
                

def main():
    if len(sys.argv) == 1:
        arguments = ['asd', 'in_fasta', 'out_fasta', '-d']
    else:
        arguments = sys.argv
    #pdb.set_trace()
    ds = '-d' in arguments
    arguments = [a for a in arguments if len(a) > 0 and a[0] != '-'][1:]

    infile, outfile = arguments[:2]
    #pdb.set_trace()
    find_reps(infile, outfile, ds)

if __name__ == '__main__':
    #c1 = Counter("Loading", 10**6)
    #c2 = Counter("Correction", 10**6)
    main()
