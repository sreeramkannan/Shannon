import pdb, heapq
def local_dot(a,b):
	n = len(a)
	sum=0
	for i in range(n):
		sum += a[i]*b[i]
	return sum


def nth_largest(n, iter):
    return heapq.nlargest(n, iter)[-1]


def path_decompose(a,b,a_true,b_true,overwrite_norm,P,use_GLPK, sparsity= False):
	#mb_check = 1 #if this parameter is set to 1, if the data can be set using MB, then it
	from cvxopt import matrix, solvers, spmatrix
	import numpy, copy
	from numpy import linalg as LA
	from numpy import array
	from numpy import zeros
	import operator
	solvers.options['msg_lev'] = 'GLP_MSG_OFF'
	'''print('path decompose:')
	print('a:',a)
	print('b:',b)
	print('P')
	print(P)'''
	m = len(a)
	n = len(b)
	if m==0 or n==0:
		return [[],0]
	'''if m==1 or n==1:
		answer = array(matrix(1.,(m,n)))
		print('Answer')
		print(matrix(1.,(m,n)))
		return [answer,0]'''
	#pdb.set_trace()
	#print('m,n:'+str(m)+','+str(n))
	'''SET UP THE problem correctly'''
	#Can use least squres, but currenlty using a hack to fix
	
	#a = max(a,0)
	#b = max(b,0)
	for i in range(m):
		a[i]=max(a[i],1e-10)
	for j in range(n):
		b[j]=max(b[j],1e-10)
		
	
	'''if sparsity:
		an = nth_largest(sparsity,a)
		bn = nth_largest(sparsity,b)
	else: 	
		an =0; bn = 0

	for i in range(m):
                a[i]=max(a[i],1e-10) if a[i]>=an else 0
        for j in range(n):
                b[j]=max(b[j],1e-10) if b[j]>=bn else 0'''



	if sum(a)>sum(b):
		const = (sum(a)-sum(b))/n
		b = [k+const for k in b ]
		#b[-1] += sum(a)-sum(b)
	else:
		const = (sum(b)-sum(a))/m
		a = [k+const for k in a ]


	if m==1:
		answer = array(matrix(b,(m,n)))
		return [answer,0]
	elif n==1:
		answer = array(matrix(a,(m,n)))
		return [answer,0] 


	
	#raw_input()
		
		#a[-1] += sum(b)-sum(a)





	A = matrix(0.,(m+n-1,m*n))
	p=matrix(0.,(m*n,1))
	
	for i in range(m):
		for j in range(n):
			A[i,j*m+i] = 1.;  #check indexing
			p[j*m+i] = P[i,j];	
	
	for j in range(n-1):  #Must be range(n) to get full mattrix
		for i in range(m):
			A[m+j,j*m+i] = 1.;
	z = a +b
	#z.append(b)
	


	#print(A)
	rhs = matrix(z)
	rhs = rhs[0:n+m-1,0]
	#print(rhs)

	weight = LA.norm(a,1)
	eps = 0.001
	tol = eps*weight  #test for significance
	sparsity_factor = 0.4  #very aggressive curently - revert to 0.1 later
	removal_factor = 0.4;
	scale = max(max(rhs),1e-100) * 0.01;
	#print(rhs)
	trials = int(round(min(2*m*n*max(m,n),100)))
	curr_min = m*n +1
	curr_ans = [];
	curr_err = 0;
	curr_mult = 0; #multiplicity of current solution
	curr_on_unknown = 0;
	for ctr in range(trials):
		#print("Iam here")
		c=matrix(abs(numpy.random.normal(0,1,(m*n,1))))
		for i in range(m*n):
			c[i] = c[i]*p[i]

		#print('Cost function')
		#print(copy)

		#print(c)
		G = spmatrix(-1.,range(m*n),range(m*n))
		h = matrix(0.,(m*n,1))
		#print(G)
		#pdb.set_trace()
		if use_GLPK: 
			sol = solvers.lp(c =c ,G=G,h=h,A=A,b=rhs/scale,solver='glpk')
		else:		
                        sol = solvers.lp(c =c ,G=G,h=h,A=A,b=rhs/scale)
		temp_sol = array(sol['x'])*scale
		#print('m,n',m,n)
		another_sol = copy.deepcopy(temp_sol)
		if overwrite_norm:
			#the true values of copy count are used to decide thresholding behavior
			a = a_true[:]
			b = b_true[:]
		for i in range(m):
			for j in range(n):
				if another_sol[j*m+i]<sparsity_factor*min(a[i],b[j]):
					another_sol[j*m+i]=0
					
				if temp_sol[j*m+i]<removal_factor*min(a[i],b[j]) or temp_sol[j*m+i]<tol: # temporarily disabled
					temp_sol[j*m+i]=0
					another_sol[j*m+i]=0

				if another_sol[j*m+i]<0:
					another_sol[j*m+i]=0
					temp_sol[j*m+i]=0
		
		
		s=0
		for i in range(m*n):
			if p[i] > 0:  #Only couont the paths that are not supported
				s=s+numpy.count_nonzero(another_sol[i])
		#indices = (another_sol<noise_thresh)
		#print('an_sol:',another_sol)
		#print('temp_sol:',temp_sol)
		#another_sol[indices] = 0; 
		
		#temp_sol = another_sol  #whether you want to exterminate small weight paths or not

		# Correct the solution if it has any incoming edge unassigned
		# for i in range(m):
		# 	if sum(another_sol[i][:]) == 0:

		# 	for j in range(n):



		#row_totals = [ sum(x1) for x1 in my_list ]
		#col_totals = [ sum(x2) for x2 in zip(*my_list) ]




		#print('where')
		#print(temp_sol)
		

		#s = numpy.count_nonzero(another_sol)
		#print(s)

                if s<curr_min:
                        curr_min = s
                        curr_ans = temp_sol[:]
                        curr_mult = 0
                        curr_on_unknown = local_dot(array(p),temp_sol)

                else:
                        if s==curr_min:
                                #curr_ans = temp_sol
                                if LA.norm(curr_ans-temp_sol) > tol:
                                        curr_mult = curr_mult +1
                                if curr_ans==[] or ( abs(sum(sum(temp_sol))-sum(sum(curr_ans)))<tol and (local_dot(array(p),temp_sol) < curr_on_unknown) ) or sum(sum(temp_sol)) > sum(sum(curr_ans)):
                                        curr_ans = temp_sol[:]
                                        curr_on_unknown = local_dot(array(p),temp_sol)



		'''if s<curr_min:
		 	curr_min = s
		 	curr_ans = temp_sol
		 	curr_mult = 0
		 	curr_err = local_dot(array(p),temp_sol)

		else:
		 	if s==curr_min:
		 		#curr_ans = temp_sol
		 		if LA.norm(curr_ans-temp_sol) > tol:
		 			curr_mult = curr_mult +1
		 		if local_dot(array(p),temp_sol) < curr_err:
	 				curr_ans = temp_sol
					curr_err = local_dot(array(p),temp_sol)'''	 				
	#print(curr_ans)

	answer = matrix(0.,(m,n))
	if len(curr_ans) < m*n:
		pdb.set_trace()
	for i in range(m):
		for j in range(n):
			answer[i,j]=float(curr_ans[j*m+i])
	#print('Answer:')
	#print(answer)
    
	answer = array(answer)
	non_unique = 0
	if curr_mult > 1:
		non_unique =1
    
	if sparsity != False:
		if m*n > sparsity:    
			tmp_dict = {}   
			for i in range(m):
				for j in range(n):
					tmp_dict[(i,j)]=answer[i, j]
			#pdb.set_trace()
			sorted_tmp = sorted(tmp_dict.items(), key=operator.itemgetter(1))[::-1]
			sorted_tmp = sorted_tmp[:sparsity]
        
			new_ans = zeros((m, n))
			for ind in sorted_tmp:
				new_ans[ind[0][0], ind[0][1]] = ind[1]
			answer = new_ans 
	return [answer,non_unique]




				

'''a1=5.; a2=7.; a3=9.; a4 = 11.; a5= 15; a6 = 13;
a = [a1,a2,a3,a4+a5+a6]
b = [a1+a2+a3,a4,a5,a6]

a = [a1,a2,a3]
b = [a1,a2+a3]
import numpy;
P = numpy.array([[1, 1],[1, 1],[1,1]])

[a,z] = path_decompose(a,b,a,b,False,P,False,3)
print(z)
print(a)'''
