from mpi4py import MPI

N = 1000000
h = 1.0 / N
s = 0.0
start = MPI.Wtime()
for i in range(N):
    x = h * (i + 0.5)
    s += 4.0 / (1.0 + x**2)
pi = h * s
end = MPI.Wtime()
print("Calculation time is ", end - start)
print("Calculated pi is ", pi)
