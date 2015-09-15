import subprocess
import os
import pdb
import math

def cut_file(in_name,out_name,line_start,line_end):
	os.system('awk \'NR > ' + str(line_end) + ' { exit } NR >= ' + str(line_start) +  '\' '+ in_name + ' > ' + out_name )

def parallel_blat(target_fasta,query_fasta,out_file,QUERY_SPLIT):
	'''Function takes in target,query and output file. parallelizes blat by running GNU parallel
	- Currently only parallelizes on query space
	- Also assumes that query fasta file takes two lines per sequence (not wrapped)'''
	target_length = float(subprocess.check_output('grep -c \'>\' ' + target_fasta,shell=True))
	query_length = float(subprocess.check_output('grep -c \'>\' ' + query_fasta,shell=True))
	os.system( 'awk \'/^>/{print s? s"\\n"$0:$0;s="";next}{s=s sprintf("%s",$0)}END{if(s)print s}\' ' +query_fasta + ' > '+ query_fasta +'_nospace')
	#os.system( 'awk \'/^>/{print s? s"\\n"$0:$0;s="";next}{s=s sprintf("%s",$0)}END{if(s)print s}\' ' +target_fasta + ' > '+ target_fasta)
	query_fasta = query_fasta + '_nospace'
	#TARGET_SPLIT = 1
	#QUERY_SPLIT = 4
	#Alernately
	#QUERY_SPLIT = min(int(math.ceil(float(query_length)/float(target_length))),50)
	#QUERY_SPLIT = max(int(math.ceil(float(query_length)/float(target_length))),500) 
	#QUERY_SPLIT = int(min(QUERY_SPLIT,query_length))
	#QUERY_SPLIT= 100
	#pdb.set_trace()
	print('Query length: ' +str(query_length) + ' Target length: ' + str(target_length) + ' Query Split: ' + str(QUERY_SPLIT))
	split_size = int(math.floor(float(query_length)/QUERY_SPLIT))
	'''if split_size % 2 !=0:
		split_size +=1'''
	'''if query_length <= float(target_length):
		print('Cannot parallelize on query. Running Vanilla Blat')
		os.system('blat -noHead ' + target_fasta + ' ' +  query_fasta + ' ' + out_file)	
		return'''

	for n in range(QUERY_SPLIT):
		if n==QUERY_SPLIT-1:
                	cut_file(query_fasta,query_fasta+'_'+str(n+1),2*(n)*split_size+1,2*query_length)
		else:
			cut_file(query_fasta,query_fasta+'_'+str(n+1),2*(n)*split_size+1,2*(n+1)*split_size)
	#pdb.set_trace()
	q_range = range(QUERY_SPLIT)
	x = [int(i)+1 for i in q_range]
	q_str = " ".join(map(str,x))
	os.system('rm ' + out_file + '_*' )
	print('parallel blat -noHead ' + target_fasta + ' ' + query_fasta + '_{} ' +out_file + '_{} ::: ' + q_str )
	os.system('time parallel blat -noHead ' + target_fasta + ' ' + query_fasta + '_{} ' +out_file + '_{} ::: ' + q_str  )
	#os.system('parallel blat ' + target_fasta + ' ' + query_fasta + '_{} ' +out_file + '_{} ::: {1..' + str(QUERY_SPLIT) + '}' )
	#os.system('sort -k 10 ' + out_file + '_* > ' + out_file)
	os.system('cat ' + out_file + '_* > ' + out_file)
	os.system('rm ' + out_file + '_*' )
	os.system('rm ' + query_fasta + '_*' )

def main():
    import sys
    if len(sys.argv) == 5:
	Query_split = int(sys.argv[4])
    else: 
	Query_split = 100
    parallel_blat(sys.argv[1],sys.argv[2],sys.argv[3],Query_split)


if __name__ == '__main__':
    main()

