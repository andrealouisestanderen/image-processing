import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import colors

t_small_img_mpi = [0.247, 0.336, 0.446, 0.618, 0.710]
t_large_img_mpi = [4.047, 3.152, 3.283, 3.607, 4.017]
t_small_img_omp = [0.087, 0.085, 0.085, 0.084, 0.093]
t_large_img_omp = [1.911, 1.991, 1.963, 1.913, 1.907]


df = pd.DataFrame({
    'x_values': range(2, 11, 2),
    'y1_values': t_small_img_mpi,
    'y2_values': t_small_img_omp,
    'y3_values': t_large_img_mpi,
    'y4_values': t_large_img_omp
})

plt.plot('x_values', 'y1_values', data=df,
         marker='', color=colors.GREEN, linewidth=2, label='MPI')
plt.plot('x_values', 'y2_values', data=df,
         marker='', color=colors.RED, linewidth=2, label='OpenMP')
plt.legend()
plt.xlabel('Number of processors')
plt.ylabel('Execurion time')
plt.title('Speed-up over processors of parallel processing of a small image')
plt.show()
plt.close()


plt.plot('x_values', 'y3_values', data=df,
         marker='', color=colors.GREEN, linewidth=2, label='MPI')
plt.plot('x_values', 'y4_values', data=df,
         marker='', color=colors.RED, linewidth=2, label='OpenMP')
plt.legend()
plt.xlabel('Number of processors')
plt.ylabel('Execution time')
plt.title('Speed-up over processors of parallel processing of a large image')
plt.show()
plt.close()
