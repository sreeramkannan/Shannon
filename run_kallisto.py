import rc_gnu, os

def run_cmd(s1):
        print(s1); 
        os.system(s1) # + '>> temp_out.txt')


def filter_using_Kallisto(rec_file, ab_file, out_file,COV_CUTOFF,L):
    accepted_trans = set()
    with open(ab_file) as f:
        f.readline()
        for line in f:
            name, _, el, ec, weight = line.split()
            cov = float(ec)/float(el)*L
            if cov>=COV_CUTOFF: accepted_trans.add(name)
    write_now = True
    with open(out_file,'w') as out_file_ptr:
        for line in open(rec_file):
            fields = line.strip().split()
            if fields and fields[0][0] == '>': write_now = (fields[0][1:] in accepted_trans)
            if write_now: out_file_ptr.write(line)

def run_kallisto(rec_file,out_dir,reads_files,ds,kallisto_dir='kallisto',nJobs=1):
    kal_params = ' '
    if len(reads_files)==1: kal_params += ' -l 200 -s 20  --single '
    if not ds: kal_params += ' --fr-stranded '
    kal_index = rec_file + '_kal.index '
    run_cmd('kallisto index -i ' + kal_index + ' ' + rec_file)
    run_cmd('kallisto quant -i ' + kal_index + ' -o ' + out_dir + ' -t ' + str(nJobs) + kal_params + ' '.join(reads_files))
    kal_ab_file=out_dir+'/abundance.tsv'
    return kal_ab_file







if __name__ == '__main__':
    import sys
    arguments = sys.argv[1:]
    ss = '-s' in arguments; ds = not ss

    L = 100
    if '-L' in arguments:
        ind = arguments.index('-L')
        L = int(arguments[ind+1])
        arguments = arguments[:ind] + arguments[ind+2:]

    COV_CUTOFF = 2.5
    if '--cov' in arguments:
        ind = arguments.index('--cov')
        COV_CUTOFF = float(arguments[ind+1])
        arguments = arguments[:ind] + arguments[ind+2:]

    arguments = [a for a in arguments if a[0]!='-']


    reconstr_file = arguments[0]
    out_dir = arguments[1]
    reads_files = arguments[2:]
    kal_file = run_kallisto(reconstr_file,out_dir,reads_files,ds)
    filter_using_Kallisto(reconstr_file,kal_file,reconstr_file+'_kal',COV_CUTOFF,L)



