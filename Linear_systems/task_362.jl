module TridiagonalSolving

using LinearAlgebra


function tridiagsolve(a, b, c, f)
	x = zeros(length(a))
	alpha = Vector{Float64}[]
	beta = Vector{Float64}[]
	push!(alpha, -c[1]/b[1])
	push!(beta, -f[1]/b[1])
	for i in 2:(length(a)-1)
		push!(alpha, -c[i]/(b[i] + a[i]*alpha[i-1]))
		push!(beta, f[i]-a[i]*beta[i-1])/(b[i] + a[i]*alpha[i-1]))
	end
	x[end] = (f[end] - a[end]*beta[end])/(b[end] + a[end]*alpha[end])
	for i in (size(x)-1):-1:1
		x[i] = alpha[i]*x[i+1] + beta[i]
	end
	return x
end


function tridiagsolve(A::Tridiagonal, f)
	b = diag(A)
	a = diag(A, -1)
	c = diag(A, 1)
	return tridiagsolve(a, b, c, f)
end 


end #module
