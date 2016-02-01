import collections
import sys
import time

def load_kp1mers(infile):
    """Loads the list of K+1-mers and copycounts and determines K.
    Returns (kmers, K).
    """
    kp1mers = []
    with open(infile) as f:
        for line in f.readlines():
            if len(line) == 0: continue
            kp1mer, weight = line.split()
            kp1mer = kp1mer.upper()
            weight = int(float(weight))
            #weight = float(weight) ## usually int(weight)
            kp1mers.append((kp1mer, weight))
    if not kp1mers:
        return None, None
    K = len(kp1mers[0][0]) - 1
    return kp1mers, K

def convert(infile, outfile):
    print "{:s}: Creating Kmers from (K+1)-mers..".format(time.asctime())
    kp1mers, K = load_kp1mers(infile)
    if kp1mers is None:
        print('NULL kp1mers file')
        with open(outfile,'w') as f:
            f.write()
        return
    print "{:s}: {:d} K+1-mers loaded.".format(time.asctime(), len(kp1mers))

    kmers = collections.Counter()
    for kp1mer, weight in kp1mers:
    	k1, k2 = kp1mer[:K], kp1mer[-K:]
    	kmers[k1] += weight
    	kmers[k2] += weight

    print "{:s}: {:d} K-mers produced.".format(time.asctime(), len(kmers))
    with open(outfile, 'w') as f:
        for kmer, weight in kmers.items():
            f.write("{:s}\t{:d}\n".format(kmer, weight))
            #f.write("{:s}\t{:f}\n".format(kmer, weight)) ## usually {:d}

def main():
    if len(sys.argv) == 1:
        arguments = ['', 'kp1mers.dict', 'kmers.dict']
    else:
        arguments = sys.argv
    infile, outfile = arguments[1], arguments[2]
    convert(infile, outfile)

if __name__ == '__main__':
    main()
