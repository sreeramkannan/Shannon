#file = 'S15_SE_c1cat_rerun_rand_algo_output/reconstructed_pre.fasta'
#new_file = 'S15_SE_c1cat_rerun_rand_algo_output/reconstructed.fasta'

def process_concatenated_fasta(file, new_file):
    # Makes sure all transcripts names are unique so that blat can be run to compare output to reference
    f = open(file, 'r')
    new_f = open(new_file, 'w')
    lines = f.readlines()
    names_seen = {}
    for line in lines:
        tokens = line.split()
        if tokens[0][0] == ">":
            if tokens[0] in names_seen:
                new_f.write(tokens[0] + "_" + str(names_seen[tokens[0]]) + '\n')
                names_seen[tokens[0]] += 1
            else:
                new_f.write(line)
                names_seen[tokens[0]] = 1
        else:
            new_f.write(line)