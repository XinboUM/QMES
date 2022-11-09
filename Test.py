import sys
import time
import qiskit as qk
import MatrixProcedures as mp
import QAlgs as qa
import QWOps as qw
import numpy as np

#==========
# Encodign an arbitrary quantum state
#==========


# tstart=time.time()
# sim = qk.BasicAer.get_backend('statevector_simulator')

# # the choice of nq_phase affects the accuracy of QPE
# nq_phase=6
# MM=2

# # initialize the matrix equation, and prepare it for the quantum procedure
# print("Building system... ", end='') # end with space rather than an endline

# msystem=mp.MatrixSystem(M=MM,expand=False)
# msystem.RandInit(D=2)
# msystem.PrepSystem()
# print("Done.")
# msystem.PrintMatrix()


# print('numpy complex square root of 3-4i = ', np.sqrt(3-4J))

# # date
# from datetime import datetime

# now = datetime.now()
# current_time = now.strftime("%Y_%m_%d_%H_%M_%S")

# print(current_time)


# # set up the output directory
# import os 
# outfile_dir = "./Output"
# if not os.path.isdir(outfile_dir):
#     os.makedirs(outfile_dir)

# # output to file 
# outfile_name = current_time + '.txt'
# output_file = open(outfile_dir+'/'+outfile_name,'w')

# output_file.write("Test for the output system")
# output_file.close()

#%%

# a = [1, 2, 3, 4]

# print(a, ', ', len(a))

# a.append([])

# print(a, ', ', len(a))

# # %%
# import bisect
# print(bisect.bisect_left(a, 0))


# %%
M = 10
D = 3
matrix_sys = mp.MatrixSystem(M)
matrix_sys.RandInit(D)

for i in range(len(matrix_sys.A0_indices)):
    print(matrix_sys.A0_indices[i])
    print(matrix_sys.A0_elements[i])

matrix_sys.PrintMatrix()
# %%

import math

EPSILON = 1e-10
print(math.ceil(math.log(matrix_sys.M,2)-EPSILON)+1)


# %%
