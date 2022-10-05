using LinearAlgebra
using TridiagonalSolving
using TrianglesSolving


function solution_a()
	A = [8 9 -4 1; 0 4 1 0; 0 0 -1 6; 0 0 0 11]
	b = [9, 3, -1, 2]
	x = backwardsub(A, b)
	return x, b - A*x
end



function solution_b()
	B = [-2 1 0 0 0; 1 -2 1 0 0; 0 1 -2 1 0; 0 0 1 -2 1; 0 0 0 1 -2]
	b = [1, 1, 1, 1, 1]
	x = tridiagsolve(B, b)
	return x, b - B*x
end



function solution_c()
	C = [1 8 -3 9; 0 4 10 -2; 8 2 -5 1; 3 1 6 12]
	b = [3, 6, 1, 4]
	L, U, p = lu(C)
	x = backwardsub(U, forwardsub(L, b))
	return x, b - C*x
end
	
	
