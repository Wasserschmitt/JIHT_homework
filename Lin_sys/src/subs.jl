function forwardsub(L::Matrix, b::Vector)
	len = length(b)
	x = Vector{Float64}(undef, len)
	x[1] = b[1]/L[1, 1]
	for i in 2:len
		pre_index = i-1
		pre_res = dot(L[i, 1:pre_index], x[1:pre_index])
		x[i] = (b[i] - pre_res)/L[i, i]
	end
	LowerMatrixChecking(L) && return x
	error("Неподходящая матрица")
end


function backwardsub(U::Matrix, b::Vector)
	len = length(b)
	x = Vector{Float64}(undef, len)
	x[end] = b[end]/U[end, end]
	for i in (len-1):-1:1
		pre_index = i+1
		pre_res = dot(U[i, pre_index:len], x[pre_index:len])
		x[i] = (b[i] - pre_res)/U[i, i]
	end
	UpperMatrixChecking(U) && return x
	error("Неподходящая матрица")
end


function forwardsub!(x::Vector{Real}, L::Matrix{Real}, b::Vector{Real})
	if LowerMatrixChecking(L)
		len = length(x)
		x[1] = b[1]/L[1, 1]
		for i in 2:len
			pre_index = i-1
			pre_res = dot(L[i, 1:pre_index], x[1:pre_index])
			x[i] = (b[i] - pre_res)/L[i, i]
		end
	else
		error("Неподходящая матрица")
	end
	
end


function backwardsub!(x::Vector, U::Matrix, b::Vector)
	if UpperMatrixChecking(U)
		len = length(x)
		x[end] = b[end]/U[end][end]
		for i in (len-1):-1:1
			pre_index = i+1
			pre_res = dot(U[i, pre_index:len], x[pre_index:len])
			x[i] = (b[i] - pre_res)/U[i, i]
		end
	else
		error("Неподходящая матрица")
	end

end
