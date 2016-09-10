import os, subprocess, sys,pdb,time, math, multiprocessing, copy
D={'A':'T','C':'G','G':'C','T':'A','N':'N'}
reverse_complement = lambda x: ''.join([D[B] for B in x][::-1])

def run_cmd(cmd):
    print(cmd)
    os.system(cmd)


def cut_file(in_name,out_name,line_start,line_end):
    os.system('awk \'NR > ' + str(line_end) + ' { exit } NR >= ' + str(line_start) +  '\' '+ in_name + ' > ' + out_name )


def rc_gnu(infile,tempfile,outfile,nCPU):
    chunks = nCPU
    file_length = float(subprocess.check_output('grep -c \'>\' ' + infile,shell=True))
    split_size = int(math.ceil(float(file_length)/chunks))
    infile_piece = open(tempfile+'_1','w'); piece_no = 1;  curr_seqs = []
    for line in open(infile):
	curr_seqs.append(line); 
	if len(curr_seqs)==split_size*2: 
		infile_piece.write(''.join(curr_seqs))
		infile_piece.close()
		piece_no +=1
		infile_piece = open(tempfile+'_'+str(piece_no),'w')
	        curr_seqs = []; 
    infile_piece.write(''.join(curr_seqs))
    infile_piece.close()


    '''for n in range(chunks):
        if n==chunks-1:
            cut_file(infile,infile+'_'+str(n+1),2*(n)*split_size+1,2*file_length)
        else:
            cut_file(infile,infile+'_'+str(n+1),2*(n)*split_size+1,2*(n+1)*split_size)'''
    chunks = piece_no; c_range = range(chunks)
    x = [int(i)+1 for i in c_range]
    c_str = " ".join(map(str,x))
    os.system('parallel python rc_s.py ' + tempfile + '_{} ' + outfile + '_{} ' + ' ::: ' + c_str  )

    file_list = ' '.join([outfile+'_'+str(i+1) for i in range(chunks)])

    os.system('cat ' + file_list+' > ' + outfile)
    os.system('rm ' + outfile+'_* ')
    os.system('rm ' + tempfile+'_*  ')




    #pdb.set_trace()

def main():
    if len(sys.argv) == 1:
        arguments = ['asd', 'in_fasta', 'out_fasta', '-d']
    else:
        arguments = sys.argv[1:]

    infile, tempfile, outfile = arguments[:3]
    nCPU = int(arguments[3])
    rc_gnu(infile,tempfile,outfile,nCPU)

if __name__ == '__main__':
    #c1 = Counter("Loading", 10**6)
    #c2 = Counter("Correction", 10**6)
    main()


