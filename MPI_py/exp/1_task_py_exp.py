from mpi4py import MPI
import math
from decimal import *
getcontext().prec = 50

start, ftol, res, cur_num, n = MPI.Wtime(), 1e-50, Decimal(0.0), Decimal(1.0), Decimal(1.0)
while cur_num > ftol:
	res, n, cur_num = res + cur_num, n + 1, cur_num/n
end = MPI.Wtime()
print(res, end - start)

