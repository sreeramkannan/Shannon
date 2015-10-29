import os, pdb
import sys, time

def run_cmd(s1):
	print(s1); os.system(s1)

def find_representatives(file_in,file_out):
	'''Find representative transcripts from file_in and return them in file_out.'''
	blat_file = file_out+ '_temp_blat.psl'
	thresh_factor = 0.9 
	os.system('python parallel_blat.py '+ file_in + ' '+ file_in + ' ' + blat_file + ' 100')
	allTr ={}
	#blat_file = '../Shannon/Alpha/Cancer/shannon_temp.psl'
	with open(blat_file,'r') as fp:
		i = 0
		for line in fp:
			if i <0:
				i+=1; continue;
			tokens = line.split();qName = tokens[9]; tName = tokens[13]; qLen = int(tokens[10]); tLen = int(tokens[14]); score = int(tokens[0]); 
			#allTr[qName] = allTr.get(qName,0) +1; 
			#allTr[tName] = allTr.get(tName,0) +1; 
			if qName == tName:
				continue
			if qLen < 200:
				allTr[qName] = -1
			if score > thresh_factor*qLen:
				if (tLen > qLen) or ((tLen == qLen) and (tName > qName)):
					allTr[qName] = -1
					#pdb.set_trace()
			i+=1
	# curr_head = ''
	# with open(file_in,'r') as fp_in:
	# 	for line in fp_in:
	# 		tokens=line.strip().split();
	# 		if tokens[0][0]=='>':
	# 			tr_name=tokens[0][1:]
	# 			val = allTr[tr_Name]


	with open(file_out,'w') as fp_out:
		write_now = False
		with open(file_in,'r') as fp_in:
			for line in fp_in:
				tokens=line.strip().split();
				if tokens[0][0]=='>':
					tr_name=tokens[0][1:]
					if allTr.get(tr_name,0) < 0:
						write_now = False
					else:
						write_now = True
				if write_now:
					fp_out.write(line)


	# with open(file_out,'w') as fp_out:
	# 	for tr,val in allTr.iteritems():
	# 		if val < 0:
	# 			fp_out.write('>'+tr+'\n')
	# 			fp_out.write()

#find_representatives('./SRR/Shannon/shannon.fasta','./SRR/Shannon/shannon_short.fasta')
#find_representatives('../Shannon/Alpha/Cancer/shannon.fasta','../Shannon/Alpha/Cancer/shannon_short.fasta')
find_representatives('./Yeast/Shannon.fasta','./Yeast/Shannon_reps.fasta')

#find_representatives('../josephhui/shannon/ww/analysis/reconstructed.fa','../josephhui/shannon/ww/analysis/reconstructed_trimed.fa')
#find_representatives('./Quorum_20M_p20r2u1000/Quorum_20M_p20r2u1000_allalgo_output/reconstructed.fasta','./Quorum_20M_p20r2u1000/Quorum_20M_p20r2u1000_allalgo_output/reps.fasta')
#find_representatives('./Snyder/Quorum_algo_input/c2.fasta','./Snyder/Quorum_algo_input/c2_reps.fasta')
#find_representatives('./Ebola_algo_output/Run2/Trinity/inchworm.K25.L25.DS.fa','./short_inch.fa')
