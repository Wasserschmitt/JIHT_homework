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
