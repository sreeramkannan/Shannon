import time
import sys
import re
import pdb,math
import numpy

BASES = ['A', 'G', 'C', 'T']
correct_errors = True
rmer_to_contig = {}
contig_to_rmer = {}
cmer_to_contig = {}

class Counter():
    def __init__(self, name, report_length):
        self.name = name
        self.count = 0
        self.report_length = report_length

    def increment(self):
        self.count += 1
        if self.count % self.report_length == 0:
            print "{:s}: {:s}, processed {:d}".format(time.asctime(), self.name, self.count)

c1 = Counter("Loading", 10**6)
c2 = Counter("Correction", 10**6)


def reverse_complement(bases):
    """Return the reverse complement of BASES. Assumes BASES is
    all uppercase.
    """
    replacements = [('A', 't'), ('T', 'a'), ('C', 'g'), ('G', 'c')]
    for ch1, ch2 in replacements:
        bases = re.sub(ch1, ch2, bases)
    return bases[::-1].upper()

def argmax(lst, key):
    """Returns the element x in LST that maximizes KEY(x).
    """
    best = lst[0]
    for x in lst[1:]:
        if key(x) > key(best):
            best = x
    return best

def load_kmers(infile, double_stranded):
    """Loads the list of K-mers and copycounts and determines K.
    Returns (kmers, K).
    """
    #c1 = Counter("Loading", 10**6)
    

    kmers = {}
    with open(infile) as f:
        for line in f:
            if len(line) == 0: continue
            c1.increment()
            kmer, weight = line.split()
            kmer = kmer.upper()
            weight = (float(weight))
            kmers[kmer] = kmers.get(kmer,0)+weight
            if double_stranded:
                kmers[reverse_complement(kmer)] = kmers.get(reverse_complement(kmer),0)+weight

    K = len(kmers.keys()[0])
    return kmers, K

def extend(start, extended, unextended, traversed, kmers, K):
    last = unextended(start)
    tot_weight = 0; tot_kmer=0
    extension = []
    
    while True:
        next_bases = [b for b in BASES if extended(last, b) in kmers and extended(last, b) not in traversed]
        if len(next_bases) == 0: return [extension,tot_weight,tot_kmer]
        
        next_base = argmax(next_bases, lambda b : kmers[extended(last, b)])
        extension.append(next_base); tot_weight += kmers[extended(last,next_base)]; tot_kmer +=1
        last = extended(last, next_base)
        traversed.add(last)
        last = unextended(last)
        c2.increment()

def extend_right(start, traversed, kmers, K):
    return extend(start, lambda last, b: last + b, lambda kmer: kmer[-(K - 1):],
        traversed, kmers, K)

def extend_left(start, traversed, kmers, K):
    return extend(start, lambda last, b: b + last, lambda kmer: kmer[:K - 1],
        traversed, kmers, K)

def duplicate_check(contig, r = 12, f = 0.5):
    # To add: if rmer in the contig multiple times, only increment the dup-contig once for each time its in dup-contig
    dup_count = {}
    max_till_now = 0; max_contig_index = -1
    for i in range(0, len(contig)-r+1):
        if contig[i:i+r] in rmer_to_contig:
            for dup in rmer_to_contig[contig[i:i+r]]:
                if dup not in dup_count:
                    dup_count[dup] = 1
                else:
                    dup_count[dup] += 1
                if dup_count[dup] >= max_till_now:
                    max_till_now = dup_count[dup]; max_contig_index = dup
    a = numpy.zeros(len(contig))
    for i in range(0, len(contig)-r+1):
        if contig[i:i+r] in rmer_to_contig:
            if max_contig_index in rmer_to_contig[contig[i:i+r]]:
                a[i:i+r] = 1

    if 1: #for dup in dup_count:
        if sum(a)> f* float(len(contig)):
            return True
        else:
            return False
           


def trim_polyA(contig):
    '''Trim polyA tails if atleast last minLen bases are polyA or first minLen bases are polyT'''
    minLen = 12; startLen = 0; endLen = 0; startPt = 0
    maxMiss = 1; startMiss = 0; endMiss = 0; endPt = 0
    for i in range(len(contig)):
        if contig[i] != 'T':
            startMiss+=1;
        else:
            startPt = startLen  + 1
        if startMiss > maxMiss:
            break
        startLen +=1
    
    totLen = len(contig)
    

    for j in range(len(contig)):
        if contig[totLen-1-j] != 'A':
            endMiss +=1;
        else:
            endPt = endLen +1
        if endMiss > maxMiss:
            break;
        endLen +=1

    if startLen < minLen:
        startLen = 0
 
    if endLen < minLen:
        endLen = 0
 
    return contig[startPt:totLen-endPt]
        
        
def run_correction(infile, outfile, min_weight, min_length,double_stranded):
    f_log = open(outfile+"__ec_log.txt", 'w')
    #pdb.set_trace()
    print "{:s}: Starting Kmer error correction..".format(time.asctime())
    f_log.write("{:s}: Starting..".format(time.asctime()) + "\n")
    kmers, K = load_kmers(infile, double_stranded)
    print "{:s}: {:d} K-mers loaded.".format(time.asctime(), len(kmers))
    f_log.write("{:s}: {:d} K-mers loaded.".format(time.asctime(), len(kmers)) + "\n")
    heaviest = sorted(kmers.items(), key=lambda kv: kv[1])
    heaviest = [(k, w) for k, w in heaviest if w >= min_weight]
    traversed, allowed = set(), set()
    f1 = open(outfile+'_contig','w')
    contig_index = 0;
    contig_connections = {}
    components_numbers = {}
    contigs = ["buffer"]
    #pdb.set_trace()
    while len(heaviest) > 0:
        start_kmer, _ = heaviest.pop()
        if start_kmer in traversed: continue
        traversed.add(start_kmer)

        [right_extension,right_kmer_wt,right_kmer_no] = extend_right(start_kmer, traversed, kmers, K)
        [left_extension,left_kmer_wt,left_kmer_no] = extend_left(start_kmer, traversed, kmers, K)
        tot_wt = right_kmer_wt + left_kmer_wt + kmers[start_kmer];
        tot_kmer = right_kmer_no + left_kmer_no + 1;
        avg_wt = tot_wt/max(1,tot_kmer);
        contig = ''.join(reversed(left_extension)) + start_kmer + ''.join(right_extension)
        #contig = trim_polyA(contig)

        r = 12
        duplicate_suspect = duplicate_check(contig, r) #Check whether this contig is significantly represented in a earlier contig. 
        #pdb.set_trace()
        # The line below is the hyperbola error correction.
        if (not correct_errors) or (len(contig) >= min_length and len(contig)*math.pow(avg_wt,1/4.0) >= 2*min_length*math.pow(min_weight,1/4.0) and not duplicate_suspect):
            f1.write("{:s}\n".format(contig));  contig_index+=1
            contigs.append(contig)
            if contig_index not in contig_connections:
                contig_connections[contig_index] = {}
            # For adding new kmers
            for i in range(len(contig) - K +1):
                allowed.add(contig[i:i+K])
            # For error correction
            for i in range(len(contig) -r + 1):
                if contig[i:i+r] in rmer_to_contig:
                    rmer_to_contig[contig[i:i+r]].append(contig_index)
                else:
                    rmer_to_contig[contig[i:i+r]] = [contig_index]

    with open(outfile, 'w') as f:
        for kmer in allowed:
            f.write("{:s}\t{:d}\n".format(kmer, int(kmers[kmer])))
    del kmers
    f_log.write("{:s}: {:d} K-mers written to file.".format(time.asctime(), len(allowed)) + " \n")


    f1.close()
    print "{:s}: {:d} K-mers remaining after error correction.".format(time.asctime(), len(allowed))
    f_log.write("{:s}: {:d} K-mers remaining after error correction.".format(time.asctime(), len(allowed)) + " \n")
    f_log.close()



   
def extension_correction(arguments):
    #pdb.set_trace()
    double_stranded = '-d' in arguments
    arguments = [a for a in arguments if len(a) > 0 and a[0] != '-']

    infile, outfile = arguments[:2]
    min_weight, min_length = int(arguments[2]), int(arguments[3])
    #pdb.set_trace()
    allowed_kmer_dict = run_correction(infile, outfile, min_weight, min_length, double_stranded)
    return allowed_kmer_dict

if __name__ == '__main__':
    #c1 = Counter("Loading", 10**6)
    #c2 = Counter("Correction", 10**6)
    if len(sys.argv) == 1:
        arguments = ['kmers.dict', 'allowed_kmers.dict', '1', '1', '-d']
    else:
        arguments = sys.argv[1:]

    extension_correction(arguments)
 
