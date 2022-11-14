#условия

z_1 = [0.9, 0.1]
K_1 = [1.5, 0.01]

z_2 = [0.2463, 0.2208, 0.2208, 0.3121]
K_2 = [40, 25, 0.6, 0.005]

#описание функции


rach_rice_eq_one(z::Float64, K::Float64, x::Float64) = z*(K-1)/(x*(K-1)+1)

rach_rice_eq(z::Vector{Float64}, K::Vector{Float64}, x::Float64) = sum([rach_rice_eq_one(z[i], K[i], x) for i in 1:length(z)])


#нахождение корня уравнения


function bisection(f, x₁, x₂; xtol=eps(), ftol=eps())
    if x₁ > x₂; x₁, x₂ = x₂, x₁; end
    y₁, y₂ = f(x₁), f(x₂)

    sign(y₁) == sign(y₂) && error("Функция должна иметь разные знаки в концах отрезка")
    abs(y₁) < ftol && return x₁
    abs(y₂) < ftol && return x₂
    
    maxiter = ceil(Int, log2((x₂-x₁)/xtol))
    
    for i in 1:maxiter
        xnew = (x₂ + x₁) / 2
        ynew = f(xnew)
        
        if sign(y₂) == sign(ynew)
            x₂, y₂ = xnew, ynew
        elseif sign(y₁) == sign(ynew)
            x₁, y₁ = xnew, ynew
        else
            return xnew
        end
        abs(ynew) < ftol && return xnew
    end
    return (x₂ + x₁)/2
end
 
function root_search(z, K) 
    cur_func(x) = rach_rice_eq(z, K, x)
    return bisection(cur_func, 1/(1-maximum(K)), 1/(1-minimum(K)))
end


#вывод данных для построения графика


function dumping(z, K, file)
    for i in 1:999
	arg = (1/(1-minimum(K)) - 1/(1-maximum(K)))/1000*i + 1/(1-maximum(K))
	println(file, join([arg rach_rice_eq(z, K, arg)], ' '))
	end
    close(file)
    return nothing
end
    

#Функция полной обработки данных


function task_solver(z, K, file)
    dumping(z, K, file)
    return root_search(z, K)
end


#
file_1 = open("dump1.txt", "w")
file_2 = open("dump2.txt", "w")

println(task_solver(z_1, K_1, file_1))
println(task_solver(z_2, K_2, file_2))
    

