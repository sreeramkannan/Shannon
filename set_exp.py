
def set_exp(range_start,range_end,sparsity,option,sample_name):
  import random
  #range_start = 20557339;
  #range_end = 20563434;
  #option ='loguniform'

  nline =0
  nfield=7
  totalweight=0 
  allids={};
  allidl=[];
  allexp=[];


  output_file = open('./'+ sample_name +'algo_input/random_out.exp','w')
  for lines in open('./sim_input/random.exp'):
    nline=nline+1;
    if lines[0]=='#':
      output_file.write(lines)
      continue;
    fields=lines.strip().split();
    #print(fields)
    if len(fields)<nfield+1:
      print('Error: the annotation file should include at least '+str(nfield+1)+' fields.');
      sys.exit();
    [l,r] = fields[4].split(':')
    [ch_start,ch_end] = r.split('-')
    ch_start = int(ch_start)
    ch_end = int(ch_end)
    tr_length = float(fields[1])
    if ch_start >= range_start and ch_end <= range_end:
      if option == 'equal':
        exp_val = 1
      elif option == 'loguniform':
        exp_val = 10**(random.uniform(-3,0))
      elif option == 'lognormal':
        exp_val = 10**(random.normalvariate(0,1))
      else:
        exp_val = float(fields[nfield])
    else:
      exp_val = 0
    
    exp_val = exp_val * tr_length
    if random.uniform(0,1)<(1-sparsity): #sparsity is the fraction of non-zeros
      exp_val=0
    vec = fields[0]+'\t'+fields[1]+'\t'+fields[2]+'\t'+fields[3]+'\t'+fields[4]+'\t'+fields[5]+'\t'+fields[6]+'\t'+str(exp_val)+'\n'
    output_file.write(vec)
    #print(ch_start,ch_end)

  output_file.close() 
  print('Read %d lines of the annoatation' % nline);

