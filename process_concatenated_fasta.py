#file = 'S15_SE_c1cat_rerun_rand_algo_output/reconstructed_pre.fasta'
#new_file = 'S15_SE_c1cat_rerun_rand_algo_output/reconstructed.fasta'



def process_concatenated_fasta(file, new_file, ds):
    # Makes sure all transcripts names are unique so that blat can be run to compare output to reference
    D = {'A':'T','C':'G','G':'C','T':'A','N':'N'}
    reverse_complement = lambda x: ''.join([D[B] for B in x][::-1])

    f = open(file, 'r')
    new_f = open(new_file, 'w')
    lines = f.readlines()
    contigs = {}
    names_seen = {}
    last_name = ''
    for line in lines:
        tokens = line.split()
        if tokens[0][0] == ">":
            if tokens[0] in names_seen:
                last_name = (tokens[0] + "_" + str(names_seen[tokens[0]]) + '\t'.join(tokens[1:]) + '\n')
                names_seen[tokens[0]] += 1
            else:
                last_name=(line)
                names_seen[tokens[0]] = 1
        elif len(line) >  200:
            curr_contig=line.strip()
            if curr_contig in contigs: continue
            if ds and reverse_complement(curr_contig) in contigs: continue
            contigs[curr_contig]=1
	    new_f.write(last_name)
            new_f.write(line)


def main():
    if len(sys.argv) == 1:
        arguments = ['', 'in_file.fasta', 'out_file.fasta']
    else:
        arguments = sys.argv

    ds = '-d' in arguments
    arguments = [a for a in arguments if len(a) > 0 and a[0] != '-'][1:]

    infile, outfile = arguments[1], arguments[2]
    process_concatenated_fasta(infile, outfile, ds)

if __name__ == '__main__':
    main()
