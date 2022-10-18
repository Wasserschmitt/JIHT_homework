function intadapt(f, a, b, tol, xtol=eps(), fa=f(a), fb=f(b), m=(b-a)/2, fm=f(m))
    if a > b; a, b = b, a; end
    
    xl = (a + m)/2; fl = f(xl)  # расположение:
    xr = (m + b)/2; fr = f(xr)  # a -- xl -- m -- xr -- b
    
    T = Vector{Float64}(undef, 3)
    h = b - a
    T[1] = h * (fa + fb)/2
    T[2] = T[1]/2 + h/2 * fm
    T[3] = T[2]/2 + h/4 * (fl + fr)
    S = (4*T[2:end] - T[1:2]) / 3

    err = (S[2] - S[1]) / 15
    
    if abs(err) < tol * (1 + tol * abs(S[2]))
        Q = S[2]
        nodes = [a, xl, m, xr, b]
    else
        b - a ≤ xtol && error("Достигнут предел точности отрезка интегрирования `xtol`.")
        Ql, nodesl = intadapt(f, a, m, tol, xtol, fa, fm, xl, fl)
        Qr, nodesr = intadapt(f, m, b, tol, xtol, fm, fb, xr, fr)
        Q = Ql + Qr
        nodes = [nodesl; nodesr[2:end]]
    end
    return (Q, nodes)
end

foo1(x) = sin(100*x*exp(-x^2))

int_1, nodes_1 = intadapt(foo1, -1/3, 3, 1e-6)
println(int_1)

file = open("52.txt", "w")
foo2(x) = exp(x)*cos(x)
value = (exp(pi/2)-1)/2
for i in 0:10
	integr, nodes_2 = intadapt(foo2, 0, pi/2, 10.0^(-2-i))
	println(file, join([length(nodes_2) abs(integr-value)], ' '))
end
close(file)
