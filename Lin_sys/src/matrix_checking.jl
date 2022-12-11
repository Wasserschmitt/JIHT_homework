function UpperMatrixChecking(U)
	return (IsSquared(U) && IsUpperTriangle(U) && IsNotDegenerated(U))
end


function LowerMatrixChecking(L)
	return (IsSquared(L) && IsLowerTriangle(L) && IsNotDegenerated(L))
end


function IsSquared(A)
	return (size(A, 1) == size(A, 2))
end


function IsUpperTriangle(U)
	u_size = size(U, 1)
	if u_size == 1
		return true
	else
		return (zeros(u_size - 1) == U[2:end, 1]) && IsUpperTriangle(U[2:end, 2:end])
	end
end


function IsLowerTriangle(L)
	l_size = size(L, 1)
	if l_size == 1
		return true
	else
		return (zeros(l_size - 1) == L[1:(l_size-1), end]) && IsLowerTriangle(L[1:(l_size-1), 1:(l_size-1)])
	end
end


function IsNotDegenerated(A)
	return (prod(diag(A)) != 0)
end
