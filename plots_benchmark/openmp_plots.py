import numpy as np
import matplotlib.pyplot as plt

####### OpenMP parallel force with 3rd Newton law ########

##natoms = 108

labels = ['1', '2', '4', '8']

t_omp_newt_108 = [7.496, 7.694, 8.141, 8.569]
#t_omp_newt_2916 = [392.314, 402.297, 397.708, 403.149]

x1 = np.arange(len(labels))
width = 0.35

fig1, ax1 = plt.subplots()

rects1 = ax1.bar(x1, t_omp_newt_108, width, color='b')
#rects2 = ax1.bar(x1 + width/2, t_omp_newt_2916, width, label='2916 atoms', color='r')

ax1.set_title('OPENMP 3rd Newton law natoms=108 (Time vs number of threads)')
ax1.set_ylabel('Time (s)')
ax1.set_xticks(x1)
ax1.set_xticklabels(labels)
#ax1.legend()
fig1.tight_layout()
plt.savefig('openmp_implementation_108')
plt.show()

################ natoms = 2916 ####################

labels = ['1', '2', '4', '8']

#t_omp_newt_108 = [7.496, 7.694, 8.141, 8.569]
t_omp_newt_2916 = [392.314, 402.297, 397.708, 403.149]

x2 = np.arange(len(labels))
width = 0.35

fig2, ax2 = plt.subplots()

#rects1 = ax1.bar(x1 - width/2, t_omp_newt_108, width, label='108 atoms', color='b')
rects2 = ax2.bar(x2, t_omp_newt_2916, width, color='b')

ax2.set_title('OPENMP 3rd Newton law natoms=2916 (Time vs number of threads)')
ax2.set_ylabel('Time (s)')
ax2.set_xticks(x2)
ax2.set_xticklabels(labels)
#ax2.legend()
fig2.tight_layout()
plt.savefig('openmp_implementation_2916')
plt.show()
