module TrianglesSolving

using LinearAlgebra

export backwardsub, forwardsub

function forwardsub(L, b)
	len = length(b)
	x = Vector{Float64}(len)
	x[1] = b[1]/L[1][1]
	for i in 2:len
		x[i] = (b[i] - L[i][:(i-1)]*x[:(i-1)])/L[i][i]
	end
	LowerMatrixChecking(U) && return x
	error("Неподходящая матрица")

end


function backwardsub(U, b)
	len = length(b)
	x = Vector{Float64}(len)
	x[len] = b[len]/L[len][len]
	for i in (len-1):-1:1
		x[i] = (b[i] - L[i][(i+1):len]*x[(i+1):len])/L[i][i]
	end
	UpperMatrixChecking(U) && return x
	error("Неподходящая матрица")
end


function forwardsub!(x, L, b)
	if LowerMatrixChecking(U)
		len = length(x)
		x[1] = b[1]/L[1][1]
		for i in 2:len
			x[i] = (b[i] - L[i][:(i-1)]*x[:(i-1)])/L[i][i]
		end
	else
		error("Неподходящая матрица")
	end
	
end


function backwardsub!(x, U, b)
	if UpperMatrixChecking(U)
		len = length(x)
		x[len] = b[len]/L[len][len]
		for i in (len-1):-1:1
			x[i] = (b[i] - L[i][(i+1):len]*x[(i+1):len])/L[i][i]
		end
	else
		error("Неподходящая матрица")
	end

end


function UpperMatrixChecking(U)
	return IsSquared(U) && IsUpperTriangle(U) && IsIsNotDegenerated(U)
end


function LowerMatrixChecking(L)
	return IsSquared(L) && IsLowerTriangle(L) && IsIsNotDegenerated(L)
end


function IsSquared(A)
	return (size(A, 1) == size(A, 2))
end


function IsUpperTriangle(U)
	u_size = size(U, 1)
	if u_size == 1
		return true
	end
	for i in 1:u_size
		return (zeros(u_size - 1) == U[2:len][1]) && IsUpperTriangle(U[2:len][2:len])
	end
end


function IsLowerTriangle(L)
	l_size = size(L, 1)
	if l_size == 1
		return true
	end
	for i in 1:l_size
		return (zeros(l_size - 1) == L[2:len][l_size]) && IsLowerTriangle(U[:(l_size-1)][:(l_size-1)])
	end
end


function IsNotDegenerated(A)
	return (dot(A[i][i] for i in 1:size(A, 1)) != 0)
end


end #module
