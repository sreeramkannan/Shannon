import os, subprocess, sys,pdb,time, math, multiprocessing, copy
import run_parallel_cmds
D={'A':'T','C':'G','G':'C','T':'A','N':'N'}
reverse_complement = lambda x: ''.join([D[B] for B in x][::-1])

def run_cmd(cmd):
    print(cmd)
    #os.system(cmd)


def cut_file(in_name,out_name,line_start,line_end):
    os.system('awk \'NR > ' + str(line_end) + ' { exit } NR >= ' + str(line_start) +  '\' '+ in_name + ' > ' + out_name )


def rc_gnu(infile,tempfile,outfile,nCPU,python_path='python ',shannon_dir=''):
    chunks = nCPU
    if chunks ==1: 
        run_cmd(python_path + ' ' + shannon_dir + 'rc_s.py ' + infile + ' ' + outfile); 
        return


    file_length = float(subprocess.check_output('grep -c \'>\' ' + infile,shell=True))
    split_size = int(math.ceil(float(file_length)/chunks))
    infile_piece = open(tempfile+'_1','w'); piece_no = 1;  curr_seqs = []
    for line in open(infile):
	curr_seqs.append(line); 
	if len(curr_seqs)==split_size*2:
                infile_piece = open(tempfile+'_'+str(piece_no),'w') 
		infile_piece.write(''.join(curr_seqs))
		infile_piece.close()
		piece_no +=1
	        curr_seqs = []; 
    if curr_seqs:
        infile_piece = open(tempfile+'_'+str(piece_no),'w') 
        infile_piece.write(''.join(curr_seqs))
        infile_piece.close()
    else:
        piece_no-=1



    '''for n in range(chunks):
        if n==chunks-1:
            cut_file(infile,infile+'_'+str(n+1),2*(n)*split_size+1,2*file_length)
        else:
            cut_file(infile,infile+'_'+str(n+1),2*(n)*split_size+1,2*(n+1)*split_size)'''
    chunks = piece_no; c_range = range(chunks)
    x = [int(i)+1 for i in c_range]
    c_str = " ".join(map(str,x))
    cmds = []
    for i in range(chunks):
        cmds.append(python_path + ' ' + shannon_dir + 'rc_s.py ' + tempfile + '_' + str(i+1) + ' ' + outfile + '_' + str(i+1))
    run_parallel_cmds.run_cmds(cmds,chunks)

    #os.system('parallel --bibtex ' + python_path + ' ' + shannon_dir + 'rc_s.py ' + tempfile + '_{} ' + outfile + '_{} ' + ' ::: ' + c_str  )

    file_list = ' '.join([outfile+'_'+str(i+1) for i in range(chunks)])

    run_cmd('cat ' + file_list+' > ' + outfile)
    run_cmd('rm ' + outfile+'_* ')
    run_cmd('rm ' + tempfile+'_*  ')




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


