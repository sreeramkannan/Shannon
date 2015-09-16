def filter_trans(fin,fout,length):
	ntrans = 0
	tr_name = ''
	with open(fout,'w') as file_out:
		for line in open(fin,'r'):
			if line[0]=='>':
				fields=line.strip().split()
				tr_name=fields[0][1:]
				to_write = line[:]
			elif len(line) >=length:
					file_out.write(to_write)
					file_out.write(line)




def filter_duplicates(fin,fout):
	'''Remove duplicates of a FASTA file with the same name'''
	tr_dict = {}
	write = 1
	with open(fout,'w') as file_out:
		for line in open(fin,'r'):
			if line[0]=='>':
				fields=line.strip().split()
                                tr_name=fields[0][1:]
				if tr_dict.get(tr_name):
					write = 0
					continue
				else:
					write=1
					file_out.write(line[:])
			else:
				if write:
					file_out.write(line[:])
						


#filter_trans('./OneHalf_NewSE_L_100_N_1000000_algo_output/reconstructed.fasta','./OneHalf_NewSE_L_100_N_1000000_algo_output/rec_short.fasta',200)
