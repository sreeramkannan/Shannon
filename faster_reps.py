import time
import sys
import re
import pdb,math
import numpy
from sets import Set

BASES = ['A', 'G', 'C', 'T']
r=24
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


def duplicate_check_ends(contigs,contig_name, rc, f = 0.99):
        contig = contigs[contig_name];
        if rc:
                contig = reverse_complement(contig) 
        qName = contig_name; qLen = len(contig);
        first_rmer = contig[:r]; last_rmer = contig[-r:]
        if first_rmer in rmer_to_contig and last_rmer in rmer_to_contig:
                contig_dict = {}
                for c,p in rmer_to_contig[first_rmer]:
                        if c == contig_name:
                                continue

                        if c in contig_dict:
                                contig_dict[c][0] = p
                        else:
                                contig_dict[c] = [p,-1]  
                for c,p in rmer_to_contig[last_rmer]:
                        if c == contig_name:
                                continue
                        if c in contig_dict:
                                contig_dict[c][1] = p
                        else:
                                contig_dict[c] = [-1,p]
                for c in contig_dict:
                        if contig_dict[c][0] >= 0 and contig_dict[c][1]>=0:
                                diff = contig_dict[c][1] - contig_dict[c][0]
                                if abs(diff - (len(contig)-r)) < 3: #we are assuming that if two kmers match and if the length is the same, then the sequences must be the same
                                        if len(contig) < len(contigs[c]) or ( (len(contig) == len(contigs[c])) and contig_name > c):
                                                print(contig_name + ' is already in '+ c)
                                                #print(contig)
                                                #print(contigs[c][contig_dict[c][0]: contig_dict[c][1]+r])
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
            for i in range(len(curr_contig)-r+1): #enumerate(find_kmers(curr_contig,r,False)):
                rmer = curr_contig[i:i+r]
                if rmer in rmer_to_contig:
                    rmer_to_contig[rmer].append([curr_name,i])
                else:
                    rmer_to_contig[rmer]= [[curr_name,i]]
    with open(outfile,'w') as out_file:
	contig_no = 0
        for contig_name in contigs:
	    contig_no +=1
            duplicate_suspect = duplicate_check_ends(contigs,contig_name, False)
            if ds:
                duplicate_suspect = duplicate_suspect or duplicate_check_ends(contigs,contig_name, True)
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
