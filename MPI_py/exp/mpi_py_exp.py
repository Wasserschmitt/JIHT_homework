from mpi4py import MPI
from decimal import *

comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()

getcontext().prec = 50

tolerance = Decimal(1e-50)
s = Decimal(0.0)
num = Decimal(1.0)
count = 0

start = MPI.Wtime()

for i in range(rank):
    num = num/(i+1)
s = s + num

while num > tolerance:
    for i in range(2,size+2):
	    num = num/(rank+i+count*size)
    s = s + num
    count += 1

exp = comm.reduce(s, op=MPI.SUM, root=0) + 1
end = MPI.Wtime()
if rank == 0:
    print("Calculation time is ", end - start)
    print("Calculated e is ", exp)
        

#Я пытался доработать версию скрипта с bcast() и gather(), но у меня не вышло
#Поэтому я разработал новый скрипт, скопировав структуру скрипта для pi под MPI с reduce()
#В итоге скрипт вышел проще и намного понятнее.
#Результаты тестирования видны на графиках, все они подписаны.
