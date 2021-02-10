import numpy as np
import matplotlib.pyplot as plt

####### OpenMP parallel force with 3rd Newton law ########

##natoms = 108

M = 4
numthreads = np.arange(M)

width=0.35

t_omp_newt = (392.314, 402.297, 397.708, 403.149)

plt.bar(numthreads,t_omp_newt,color='r')
plt.xlabel('numthreads')
plt.ylabel('Time (s)')
plt.title('Scaling natoms=108')
plt.xticks(numthreads, ('1', '2', '4', '4'))

plt.show()


