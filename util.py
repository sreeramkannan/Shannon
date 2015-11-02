import time
import random
import re
import os

BASES = 'AGCT'
L = 100

class Alignment:
    def __init__(self, g_name, blocks, rev, seq, weight):
        self.g_name = g_name
        self.blocks = blocks
        self.weight = weight
        self.seq = seq
        self.reversed = rev
    def start(self):
        return self.blocks[0][0]
    def end(self):
        start, length = self.blocks[-1]
        return start + length - 1
    def bases(self):
        for start, length in self.blocks:
            for i in range(start, start+length):
                yield i
class Read:
    def __init__(self, name, seq):
        self.name = name
        self.alignments = []
        self.seq = seq
class Counter(dict):
    def __missing__(self, key):
        return float()
    def max(self):
        max_key, max_value = None, -1
        for k, v in self.items():
            if v >= max_value:
                max_key, max_value = k, v
        return max_key
    def total(self):
        return sum(v for _, v in self.items())
class Clock:
    def __init__(self):
        self.last = time.time()
    def time(self):
        cur = time.time()
        elapsed = cur - self.last
        self.last = cur
        return elapsed

def random_seq(length):
    return ''.join(random.choice(BASES) for _ in xrange(length))
def corrupt(seq, p):
    def base(b):
        if random.random() < p:
            return random.choice(BASES)
        else:
            return b
    return ''.join(base(b) for b in seq)
def random_read(seq, length):
    start = random.randrange(len(seq) - length + 1)
    return seq[start:start + length]
def distance(seq1, seq2):
    return sum(seq1[i] != seq2[i] for i in xrange(len(seq1)))
def reverse_complement(s):
    comp = {'A':'T', 'T':'A', 'G':'C', 'C':'G'}
    return ''.join(comp[c] for c in s[::-1])

def cigar_to_blocks(start, cigar):
    """Converts cigarstring into list of blocks as (start, length).
    """
    blocks = []
    for block in re.findall(r'\d+[MNS]', cigar):
        length, match = int(block[:-1]), block[-1] == 'M'
        if match:
            blocks.append((start, length))
        if block[-1] != 'S':
            start += length
    return blocks
def read_sam(filename, target=None):
    """Loads reads from a SAM file.
    Reads are normalized so that their alignments go left-to-right, whether
    they were reverse-complemented or not.

    '___A' and '___B' are the first and second reads of a mate pair; A may
    be to the left or to the right of B.
    """
    reads = {}
    with open(filename) as f:
        for line in f:
            if line[0] == '@': continue
            fields = line.split()
            name, flags, g_name, start, cigar, seq = (
                fields[0], int(fields[1]), fields[2],
                int(fields[3]) - 1, fields[5], fields[9])
            if target and g_name != target: continue

            rev = bool((flags >> 4) & 1)
            if (flags >> 6) & 1:
                name += 'A'
            elif (flags >> 7) & 1:
                name += 'B'
            else: #Hack to deal with misnamed input reads
                name += str(hash(seq))

            if name not in reads:
                reads[name] = Read(name, seq)
            elif reads[name].seq != seq and reads[name].seq != reverse_complement(seq):
                print("Error: Read '{}' has seqs {} and {}".format(
                    name, reads[name].seq, seq))
                exit()

            if ((flags >> 2) & 1) or not re.match(r'(\d+[MNS])+\Z', cigar):
                continue

            a = Alignment(g_name, cigar_to_blocks(start, cigar), rev, seq, 1.0)
            reads[name].alignments.append(a)

    return list(read for _, read in reads.items())
def split_reads(reads):
    """Splits reads into unmapped, unique, and multiply-mapped reads.
    """
    unmapped, uniques, multiples = [], [], []
    for read in reads:
        if len(read.alignments) > 1:
            multiples.append(read)
        elif len(read.alignments) == 1:
            uniques.append(read)
        else:
            unmapped.append(read)
    return unmapped, uniques, multiples
def filter_mates(read_list):
    """Given a list of reads named as "read000A" and
    "read000B", returns two lists of mate pairs. The first list has
    uniquely-mapped consistent pairs, and the second list has all other
    pairs.

    Discards incongruous mappings. A consistent mapping is one where the
    right read is reverse-complemented, and the left read is not.
    
    Each pair in UNIQUES is a pair of alignments, normalized so the first
    alignment is to the left of the second alignment. Each pair in OTHERS
    is a pair of reads.
    """
    reads = {}
    for read in read_list:
        reads[read.name] = read

    uniques, others = [], []

    for name, read in reads.items():
        if name[-1] != 'A': continue
        mate_name = name[:-1] + 'B'
        if mate_name not in reads:
            print('Error: {} not in reads'.format(mate_name))
            exit()
        mate = reads[mate_name]

        if (len(read.alignments) != 1) or (len(mate.alignments) != 1):
            others.append((read, mate))
            left, right = len(read.alignments), len(mate.alignments)
            continue

        a, a2 = read.alignments[0], mate.alignments[0]

        if a.start() > a2.start():
            a, a2 = a2, a
        if a.reversed or not a2.reversed:
            others.append((read, mate))
            continue

        uniques.append((a, a2))

    return uniques, others

def ensure_path(filename):
    pardir = os.path.dirname(filename)
    if not os.path.exists(pardir):
        os.makedirs(pardir)
def to_fasta(filename, seq, seq_name):
    seqs_to_fasta(filename, [(seq_name, seq)])
def seqs_to_fasta(filename, seqs):
    LINE_LENGTH = 50
    with open(filename, 'w') as f:
        for seq_name, seq in seqs:
            f.write(">{}\n".format(seq_name))
            for i in xrange(0, len(seq), LINE_LENGTH):
                f.write("{}\n".format(seq[i:i+LINE_LENGTH]))
def reads_to_weighted(filename, reads):
    """READS should be a list of (weight, seq) tuples. """
    name = ''.join(random.choice('abcdefghij') for _ in range(10))
    ensure_path(filename)
    with open(filename, 'w') as f:
        for i, read in enumerate(reads):
            weight, read = read
            f.write(">read_{}_{} {}\n".format(name, i, weight))
            f.write("{}\n".format(read))
def reads_to_fasta(filename, reads):
    name = ''.join(random.choice('abcdefghij') for _ in range(10))
    ensure_path(filename)
    with open(filename, 'w') as f:
        for i, read in enumerate(reads):
            f.write(">read_{}_{}\n".format(name, i))
            f.write("{}\n".format(read))
def to_fastq(filename, reads):
    ensure_path(filename)
    with open(filename, 'w') as f:
        for i, read in enumerate(reads):
            n, read = read
            f.write("@{}-{}\n".format(n, i))
            f.write("{}\n".format(read))
            f.write("+\n")
            f.write("{}\n".format('~' * len(read)))
def from_fasta(filename):
    sequences = []
    seq = []
    next_name = '_'
    with open(filename) as f:
        for line in f:
            if line[0] == '>':
                sequences.append((next_name, ''.join(seq)))
                seq = []
                next_name = line.split()[0][1:]
            else:
                seq.append(line.strip().upper())
    sequences.append((next_name, ''.join(seq)))
    return sequences[1:]

def read_weights(reads):
    hits, weights = Counter(), Counter()
    for read in reads:
        for alignment in read.alignments:
            hits[alignment.start()] += 1
            hits[alignment.end()] += 1
            for start, length in alignment.blocks:
                for i in xrange(start, start + length):
                    weights[i] += 1
    return hits, weights

def kmers(read, K):
    for i in xrange(len(read) - K + 1):
        yield read[i:i+K]

def compare(f1, f2):
    s1 = from_fasta(f1)[0]
    s2 = from_fasta(f2)[0]
    print "{} differences.".format(distance(s1, s2))

