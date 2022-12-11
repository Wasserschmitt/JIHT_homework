from mpi4py import MPI
import numpy
from decimal import *

comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()

getcontext().prec = 50

start = MPI.Wtime()

if rank==0:
	res = 0
	while True:
		data = 1
		comm.bcast(data, root=0)
		check_nums = [comm.recv(source=i, tag=i) for i in range(1, size)]
		res += sum(check_nums)
		if check_sum[-1] < 1e-50:
			end = MPI.Wtime()
			print("Calculation time is ", end - start)
			print("Calculated e is ", res)
			data = "End"
			comm.bcast(data, root=0)
			break
else:
	data=0
	while True:
		comm.bcast(data, root=0)
		if data == "End":
			break
		


exit
			
		
