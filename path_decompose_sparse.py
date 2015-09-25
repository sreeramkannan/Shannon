import pdb, heapq
def local_dot(a,b):
	## This function takes the dot product between two vectors.  
	n = len(a)
	sum=0
	for i in range(n):
		sum += a[i]*b[i]
	return sum


def nth_largest(n, iter):
    return heapq.nlargest(n, iter)[-1]


def path_decompose(a,b,a_true,b_true,overwrite_norm,P,use_GLPK, sparsity= False):
	'''This function takes in a node's information in the attempt to decompose it into the lowest number of paths that
	accounts for the flow constraints.  
	Thhe algorithm uses many trials of a randomizaed optimization model and takes the best result seen.  
	**Do not use over_write_norm becuase a_true and b_true surrently have normalization instead of copycount. 
	a: a vector of the copytcounts of the in-edges for the node.
	b: a vector of the copycounts of the out-edges for the node.  
	a_true: a vector that should have the copycount values for the in-edges but does not.  (Don't Use)
	b_true: a vector that should have the copycount values for the out-edges but does not.  (Don't Use)
	decides whether or not to  use a_true and b_true.
	This is a matrix of in-edges versus out-edges that has a 0 if there is a known path and a 1 otherwise.  
	'''

    #mb_check = 1 #if this parameter is set to 1, if the data can be set using MB, then it
	from cvxopt import matrix, solvers, spmatrix
	import numpy, copy
	from numpy import linalg as LA
	from numpy import array
	from numpy import zeros
	import operator
	solvers.options['msg_lev'] = 'GLP_MSG_OFF'
	m = len(a)  ## a is a vector of the current in edge copycount values.
	n = len(b)  ## b is a vector of the current out edge copycount values.
	## Trivial case 
	if m==0 or n==0:
		return [[],0]

	## Make all in flow values non-zero.
	for i in range(m):
		a[i]=max(a[i],1e-10)
	for j in range(n):
		b[j]=max(b[j],1e-10)
		
	

	## If the flow in does not equal the flow out, make them equal.
	if sum(a)>sum(b):
		const = (sum(a)-sum(b))/n
		b = [k+const for k in b ]
	else:
		const = (sum(b)-sum(a))/m
		a = [k+const for k in a ]

	## Trivial cases.
	if m==1:
		answer = array(matrix(b,(m,n)))
		return [answer,0]
	elif n==1:
		answer = array(matrix(a,(m,n)))
		return [answer,0] 





	A = matrix(0.,(m+n-1,m*n))  ## A is used to enforce flow constraint on decomposition.
	p=matrix(0.,(m*n,1))  ## This vector tells whether there is a known path for the pair of nodes or not.  
	
	for i in range(m):
		for j in range(n):
			A[i,j*m+i] = 1.;  #check indexing
			p[j*m+i] = P[i,j];	
	
	for j in range(n-1):  #Must be range(n) to get full mattrix
		for i in range(m):
			A[m+j,j*m+i] = 1.;
	z = a +b
	



	rhs = matrix(z)
	rhs = rhs[0:n+m-1,0]  ## this vector is used to enforce flow constraints 

	weight = LA.norm(a,1) ## (Not used)
	eps = 0.001
	tol = eps*weight  #test for significance.  Used for various purposes
	sparsity_factor = 0.4  #very aggressive curently - revert to 0.1 later
	removal_factor = 0.4;
	scale = max(max(rhs),1e-100) * 0.01;
	#print(rhs)
	trials = int(round(min(2*m*n*max(m,n),100))) ## Number of randomized trials.
	curr_min = m*n +1  ## The lowest number of non known paths seen so far.
	curr_ans = [];     ## The best solution seen so far.
	curr_err = 0;
	curr_mult = 0;         ## multiplicity of current solution
	curr_on_unknown = 0;   ## amount of flow on non "known paths".
	for ctr in range(trials):  ## randomize the the coefficients for the known non=known path flow values.
		c=matrix(abs(numpy.random.normal(0,1,(m*n,1))))
		for i in range(m*n):
			c[i] = c[i]*p[i]

		## G and h are used to make sure that the flow values are non-negative.  
		G = spmatrix(-1.,range(m*n),range(m*n))
		h = matrix(0.,(m*n,1))
		if use_GLPK: 
			sol = solvers.lp(c =c ,G=G,h=h,A=A,b=rhs/scale,solver='glpk')
		else:		
                        sol = solvers.lp(c =c ,G=G,h=h,A=A,b=rhs/scale)
		temp_sol = array(sol['x'])*scale
		another_sol = copy.deepcopy(temp_sol)
		if overwrite_norm: ## Do not use right now because a_true and b_true are wrong.
			#the true values of copy count are used to decide thresholding behavior
			a = a_true[:]
			b = b_true[:]
            
		## This loop basically sets the flow values equal to 0 under certain conditions.              
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
		
		
		s=0 ## s equals how many non-zero flows we are sending down non "known paths" in the temp solution
		for i in range(m*n):
			if p[i] > 0:  #Only couont the paths that are not supported
				s=s+numpy.count_nonzero(another_sol[i])

				## if the temporary solution is less than the current minimum solution, replace current minimum solution with
				## the temporary solution.
                if s<curr_min:
                        curr_min = s  ## current minimum value of non "known paths" in solution 
                        curr_ans = temp_sol[:]  ## answer that attains it.  
                        curr_mult = 0  ## This is how many times a solution with this many non-zero non-known paths is seen.
                        curr_on_unknown = local_dot(array(p),temp_sol)  ## This says how much flow we are sending down non "known paths"

                else:
                        if s==curr_min:
                                if LA.norm(curr_ans-temp_sol) > tol:  ## Determines whethee to classify the solutions as different.
                                        curr_mult = curr_mult +1
                                if curr_ans==[] or ( abs(sum(sum(temp_sol))-sum(sum(curr_ans)))<tol and (local_dot(array(p),temp_sol) < curr_on_unknown) ) or sum(sum(temp_sol)) > sum(sum(curr_ans)):
										## These are a few conditions that make it the temporary solution:
										## 1: curr_sol is empty.  2: total flow difference is below a threshold AND less flow is going down unknown paths.  3: total flow is greater than curr_ans total flow.
                                        curr_ans = temp_sol[:]
                                        curr_on_unknown = local_dot(array(p),temp_sol)




	answer = matrix(0.,(m,n))
	if len(curr_ans) < m*n:
		pdb.set_trace()
	for i in range(m):
		for j in range(n):
			answer[i,j]=float(curr_ans[j*m+i])
    
    # Uniqueness consideration
	answer = array(answer)
	non_unique = 0
	if curr_mult > 1:
		non_unique =1
    
    # This makes the flow solution more sparse
	if sparsity != False:
		if m*n > sparsity:    
			tmp_dict = {}   
			for i in range(m):
				for j in range(n):
					tmp_dict[(i,j)]=answer[i, j]
			sorted_tmp = sorted(tmp_dict.items(), key=operator.itemgetter(1))[::-1]
			sorted_tmp = sorted_tmp[:sparsity]
        
			new_ans = zeros((m, n))
			for ind in sorted_tmp:
				new_ans[ind[0][0], ind[0][1]] = ind[1]
			answer = new_ans 
	return [answer,non_unique]




				

'''
#Test Case

a1=5.; a2=7.; a3=9.; a4 = 11.; a5= 15; a6 = 13;
a = [a1,a2,a3,a4+a5+a6]
b = [a1+a2+a3,a4,a5,a6]

a = [a1,a2,a3]
b = [a1,a2+a3]
import numpy;
P = numpy.array([[1, 1],[1, 1],[1,1]])

[a,z] = path_decompose(a,b,a,b,False,P,False,3)
print(z)
print(a)'''
