function tridiagsolve(a::Vector, b::Vector, c::Vector, f::Vector)
	x = zeros(length(b))
	alpha = zeros(length(a))
	beta = zeros(length(a))
	alpha[1] = -c[1]/b[1]
	beta[1] = -f[1]/b[1]
	for i in 2:(length(a)-1)
		alpha[i] = -c[i]/(b[i] + a[i]*alpha[i-1])
		beta[i] = (f[i]-a[i]*beta[i-1])/(b[i] + a[i]*alpha[i-1])
	end
	x[end] = (f[end] - a[end]*beta[end])/(b[end] + a[end]*alpha[end])
	for i in (length(x)-1):-1:1
		x[i] = alpha[i]*x[i+1] + beta[i]
	end
	return x
end


function tridiagsolve(A::Matrix, f::Vector)
	b = diag(A)
	a = diag(A, -1)
	c = diag(A, 1)
	return tridiagsolve(a, b, c, f)
end 
