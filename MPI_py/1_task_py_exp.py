from mpi4py import MPI
import math

start, ftol, res, cur_num, n = MPI.Wtime(), 1e-50, 0, 1, 1
while cur_num > ftol:
	res, n, cur_num = res + cur_num, n + 1, cur_num/n
end = MPI.Wtime()
print(res, end - start)

