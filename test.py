D = {'A':'T','C':'G','G':'C','T':'A'}

rc = lambda x: ''.join([D[B] for B in x[::-1]])
comp = lambda x:''.join([D[B] for B in x])

reverse_complement = lambda x: ''.join([{'A':'T','C':'G','G':'C','T':'A', 'N':'N'}[B] for B in x][::-1])

rc2 = lambda x: ''.join([{'A':'T','C':'G','G':'C','T':'A', 'N':'N'}[B] for B in x][::-1])

import time
a=time.time()
for i in range(100000):
   r=comp('ACAAAAAAAAAAAAAAAAAAAAAAAAAGACCGGATCAACCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC')
b = time.time()
print(b-a)



a=time.time()
for i in range(100000): 
   r=rc('ACAAAAAAAAAAAAAAAAAAAAAAAAAAGACCGGATCAACCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC')
b = time.time()
print(b-a)

