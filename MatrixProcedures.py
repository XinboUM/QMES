'''
Development notes:
1. This class has a property "expand". Its value should be decided decided in the PrepSystem method.
    expand = False only when we know that matrix A0 is Hermitian in prior, and M is a power of 2.
    The user should pass in isHermitian as a parameter in the constructor instead of expand. 
    Changing this can make the code easier to use: the user only need to know if the matrix A0 is Hermitian, QMES would be a black box.
2. Check the notes for PrepSystem
'''

import math
import numpy as np
import random
import bisect

EPSILON = 1e-10
i = complex(0.0,1.0)

class MatrixSystem:
    '''
    The MatrixSystem class contains the information of the matrix and the right-hand-side (RHS) vector in the matrix equation.
    '''

    def __init__(self, M=1, expand=True):
        self.expand=expand      # flag that controls if the matrix equation will be expanded according to (63) and (64)
        self.M=M                # the given arbitrary matrix A0 is MxM
        self.n=0                # the number of qubits needed to encode a 
        self.N=0                # N = 2^n, the system matrix A after possible expansions of A0
        self.bnorm=0.0          # the length of the RHS vector |b>
        self.d=0.0              # shifting constant to ensure no negative diagonal elements in the system matrix A
        self.ap=0.0             # 
        self.X=0.0              # N|Ajk|max
        self.C=0.0              # 
        self.A0_indices=[]      # a 2-dimentional list that contains the [row][col] index for non-zero elements in A
        self.A0_elements=[]     # a 
        self.b0=[]
        self.A_indices=[]
        self.A_elements=[]
        self.b=[]
    
    def FileInit(self):
        f=open('./MS.txt','r')
        lines=f.readlines()
        f.close()
        jp=0
        for c in range(len(lines)):
            l=lines[c].split(',')
            j=int(l[0])-1
            k=int(l[1])-1
            while(j>jp):
                jp=jp+1
                self.A0_indices.append([])
                self.A0_elements.append([])
            self.A0_indices[jp].append(k)
            self.A0_elements[jp].append(float(l[2]))
        f=open('./b0.txt','r')
        lines=f.readlines()
        f.close()
        for c in range(len(lines)):
            self.b0.append(float(l[0]))
    
    # parallel strip system
    def MoMInit(self):
        l=2.0/float(self.M)
        p=[]
        
        for j in range(self.M):
            if(j<self.M/2):
                self.b0.append(1.0)
                p.append([-0.5,l*j-0.5])
            else:
                p.append([0.5,l*(j-self.M/2)-0.5])
                self.b0.append(-1.0)
        
        for j in range(self.M):
            self.A0_indices.append([])
            self.A0_elements.append([])
            for k in range(self.M):
                self.A0_indices[j].append(k)
                if(j==k):
                    self.A0_elements[j].append(-l/(1.0)*(math.log(l)-1.5))
                else:
                    self.A0_elements[j].append(-l/(1.0)*math.log(math.sqrt((p[j][0]-p[k][0])**2.0+(p[j][1]-p[k][1])**2.0)))
    

    def RandInit(self,D=1):
        '''
        Initialize the matrix A0 and RHS vector b0 in the matrix equation randomly, where A0 is a Hermitian matrix that has roughly D elements in each row.
        
        Parameter:
            D (int): controls the sparsity of the matrix A0.
        '''

        # fill_prob = float(D-1)/float(self.M)  # remove one since we'll definitely put an element on the diagonal
        
        # reserve a spot for each diagonal element of A0
        for j in range(self.M):
            self.A0_indices.append([j])
            self.A0_elements.append([])
        
        # for each row, reserve D-1 random column entries while ignoring recurrent indices such that each row has approximately D non-zero elements
        # (statistics on "D-1 nonzeros" could be better)
        for j in range(self.M):
            for c in range(D-1):
                k = random.randint(0, self.M-1)
                if(not (k in self.A0_indices[j])):
                    bisect.insort(self.A0_indices[j],k)
                    bisect.insort(self.A0_indices[k],j)
          
        # populate the matrix
        for j in range(self.M):
            for c in range(len(self.A0_indices[j])):
                k = self.A0_indices[j][c]
                if(k==j):
                    x = 2.0*(random.random()-0.5)
                    self.A0_elements[j].append(complex(x,0.0))
                elif(k>j):
                    x=2.0*(random.random()-0.5)
                    y=2.0*(random.random()-0.5)
                    self.A0_elements[j].append(complex(x,y))
                    self.A0_elements[k].append(complex(x,-y))  # order works automatically
        
        # initialize b0 to a fully random vector
        for j in range(self.M):
            x = 2.0*(random.random()-0.5)
            y = 2.0*(random.random()-0.5)
            self.b0.append(np.complex(x,y))
    
    def PrepSystem(self):
        '''
        Prepares the given arbitrary matrix A0 and RHS vector for the quantum algorithm.
        In this code, the expansion expands A0 to A according to eq(66). The case that the A0 is Hermitian is not considered, i.e., eq(67) is not implemented.
        If A0 is not expanded, i.e., expand=False, we assume that M is a power of 2.
        There are two possible operations:
        1. expanding A0 and |b> according to the fact if A0 is Hermitian and the value of M
        2. shifting A0 according to the diagonal elements of A0

        Optimization and Implementation: 
        1. This function can be optimized by computing bnorm at first. That way we can populate the system matrix A and RHS vector b with properly scaled values.  
        2. eq(67) should be implemented.
        3. The process of finding X and d should be written together.
        '''

        # determine the appropriate size of the system matrix A for quantum operation
        if(self.expand):       
            self.n = math.ceil(math.log(self.M,2)-EPSILON)+1
            self.N = int(math.pow(2.0,self.n)+EPSILON)
        else:                  
            self.n = math.ceil(math.log(self.M,2)-EPSILON)
            self.N = self.M      # This assumes that M is a power of 2
        
        self.A_indices=[]
        self.A_elements=[]
        self.b=[]
        
        # set the elements of A and b from A0 and b0 according to eq(66)
        if(self.expand):
            #----------
            # expand b0 into b
            #----------
            for j in range(self.M):
                self.b.append(self.b0[j])
            for j in range(self.M,self.N):
                self.b.append(0.0)
            
            #----------
            # expand A0 into A
            #----------
            for j in range(self.N):
                self.A_indices.append([])
                self.A_elements.append([])

            # put A0 and A0dagger into A 
            for j in range(self.M):                     
                for c in range(len(self.A0_indices[j])):# c is the non-zero element column index in the jth row of A0
                    k = self.A0_indices[j][c]
                    self.A_indices[j].append(k+self.M)
                    self.A_indices[k+self.M].append(j)
                    self.A_elements[j].append(self.A0_elements[j][c])
                    self.A_elements[k+self.M].append(np.conjugate(self.A0_elements[j][c]))

            # put the identity matrix I_N-2M into the bottom right corner of A
            for j in range(2*self.M, self.N):            
                self.A_indices[j].append(j)
                self.A_elements[j].append(1.0)
                
        # copy A0 and b0 to A and b, respectively      
        else:
            for j in range(self.M):
                self.b.append(self.b0[j])
                self.A_indices.append([])
                self.A_elements.append([])
                for c in range(len(self.A0_indices[j])):
                    self.A_indices[j].append(self.A0_indices[j][c]);
                    self.A_elements[j].append(self.A0_elements[j][c]);
        
        # normalize b and rescale A according to bnorm
        self.bnorm = math.sqrt(sum([np.absolute(x)*np.absolute(x) for x in self.b]))
        for j in range(self.N):
            self.b[j] = self.b[j]/self.bnorm
            for c in range(len(self.A_elements[j])):
                self.A_elements[j][c] = self.A_elements[j][c]/self.bnorm
        
        # find and apply the diagonal offset 
        if(self.expand): # no offset needed if system is expanded
            self.d=0.0
        else:            # if system is not expanded, use maximal in magnitude diagonal element of A0
            self.d=0.0
            # Find the maximal diagonal magnitude. O(N)
            for j in range(self.N):
                for c in range(len(self.A_indices[j])):
                    if(self.A_indices[j][c]==j):
                        if(self.d<np.absolute(self.A_elements[j][c])):
                            self.d = np.absolute(self.A_elements[j][c])
                    elif(self.A_indices[j][c]>j):
                        break
            # Apply offset to diagonal elements of A. O(N)
            for j in range(self.N):
                for c in range(len(self.A_indices[j])):
                    if(self.A_indices[j][c]==j):
                        self.A_elements[j][c] = self.A_elements[j][c]+self.d
                    elif(self.A_indices[j][c]>j):
                        break
        
        # calculate X (get upper bound on maximal in magnitude element of A first). 
        # Complexity: O(nnz of A), where O(N) <= O(nnz of A) <= O(N^2)
        self.ap = 0.0
        for j in range(self.N):
            for elem in self.A_elements[j]:
                if(self.ap<np.absolute(elem)):
                    self.ap = np.absolute(elem)
        self.X = np.absolute(float(self.N)*self.ap + EPSILON)
    
    def CompareClassical(self,sol):
        A0=[]
        for j in range(self.M):
            A0.append([])
            for k in range(self.M):
                A0[j].append(complex(0.0,0.0))
        for j in range(self.M):
            for c in range(len(self.A0_indices[j])):
                k=self.A0_indices[j][c]
                A0[j][k]=self.A0_elements[j][c]
        A0inv=np.linalg.inv(A0)
        sol_class=A0inv.dot(self.b0)
        
        # get the relative global phase
        t1=np.angle(sol[0])
        t2=np.angle(sol_class[0])
        trel=t2-t1
        
        # compare
        print('----------------------------')
        print('Solution Comparison: ')
        print()
        print('%24s|%24s'%('Quantum','Classical'))
        print('-------------------------------------------------')
        for j in range(self.M):
            print('%10f+i(%10f)|%10f+i(%10f)'%((sol[j]*np.exp(i*trel)).real,(sol[j]*np.exp(i*trel)).imag,sol_class[j].real,sol_class[j].imag))
        print()
        print('Relative Errors:')
        summ=0.0
        for j in range(self.M):
            summ=summ+np.absolute(1.0-sol[j]/sol_class[j]*np.exp(i*trel))
            print(np.absolute(1.0-sol[j]/sol_class[j]*np.exp(i*trel)))
        print()
        print('Average Relative Error: ', summ/float(self.M))
    
    def PrintMatrix(self):
        ''''
        Print the matrix A0 and A
        '''
        # Print A0
        for j in range(self.M):
            c=0
            for k in range(self.M):
                if(c<len(self.A0_indices[j])):
                    if(self.A0_indices[j][c]==k):
                        print(" %5.2f+(%5.2f) "%(self.A0_elements[j][c].real,self.A0_elements[j][c].imag),end='')
                        c=c+1
                    else:
                        print("  0.00+( 0.00) ",end='')
                else:
                    print("  0.00+( 0.00) ",end='')
            print("\n")
        print("\n")
        print("\n")
        # Print A
        for j in range(self.N):
            c=0
            for k in range(self.N):
                if(c<len(self.A_indices[j])):
                    if(self.A_indices[j][c]==k):
                        print(" %5.2f+(%5.2f) "%(self.A_elements[j][c].real,self.A_elements[j][c].imag),end='')
                        c=c+1
                    else:
                        print("  0.00+( 0.00) ",end='')
                else:
                    print("  0.00+( 0.00) ",end='')
            print("\n")

