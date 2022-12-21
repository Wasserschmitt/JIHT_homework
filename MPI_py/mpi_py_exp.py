from mpi4py import MPI
import numpy
from decimal import *

def fac(n):
    if n == 0:
        return 1
    return fac(n-1) * Decimal(n)

comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()

getcontext().prec = 50

start = MPI.Wtime()

if rank==0:
	res = Decimal("0")
	data = Decimal("1")
	check_sum = [Decimal("1")]*(size-1)
	while True:
		comm.bcast(data, root=0)
		comm.gather(check_sum, root=0)
		res += sum(check_sum) - check_sum[0]
		if check_sum[-1] < 1e-50:
			end = MPI.Wtime()
			print("Calculation time is ", end - start)
			print("Calculated e is ", res)
			data = "End"
			comm.bcast(data, root=0)
			break
		else:
			data = check_nums[-1]


else:
	data = Decimal("0")
	count = Decimal(rank)
	while True:
		comm.bcast(data, root=0)
		if data == "End":
			break
		else:
			check_sum = data/fac(count)
			count = count + size-1
			comm.gather(data, root=0)		
