using LinearAlgebra

#функция
function f1(x, temp)
    x₁, x₂ = x
    return [
        8*temp/(3*x₁-1) - 3/(x₁^2) - 8*temp/(3*x₂-1) + 3/(x₂^2),
        -temp*log((3*x₁-1)/(2*exp(-1/2))) + temp/(3*x₁-1) - 9/(4*x₁) + temp*log((3*x₂-1)/(2*exp(-1/2))) - temp/(3*x₂-1) + 9/(4*x₂),
    ]
end

#якобиан
function J1(x, temp)
    x₁, x₂ = x
    return [
        -24*temp/(3*x₁-1)^2+6/(x₁^3)                    24*temp/(3*x₂-1)^2-6/(x₂^3)                  ;
        -temp*3/(3*x₁-1)-3*temp/(3*x₁-1)^2+36/(4*x₁)^2  temp*3/(3*x₂-1)+3*temp/(3*x₂-1)^2-36/(4*x₂)^2;
    ]
end


function newtonsys(f, x, J, temp; maxiter=50, xtol=1e-6, ftol=1e-6)
    x = float(copy(x))
    δx, y = similar.((x, x))
    for i in 1:maxiter
        y .= f(x, temp)
        δx .= .- (J(x, temp) \ y)
        x .+= δx
	print(x)
        norm(δx) < xtol && return x
        norm(y) < ftol && return x
    end
    error("Превышено число итераций.")
end


#


file = open("dump.txt", "w")

for i in 1:8
    x = [0.5, 2]
    cur_temp = 0.85+0.2*(i-1)
    res = newtonsys(f1, x, J1, cur_temp)
    pushfirst!(res, cur_temp)
    println(file, join(res, ' '))
end
close(file)

