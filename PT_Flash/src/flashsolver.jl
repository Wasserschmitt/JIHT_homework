struct SplitResult{T}
    converged::Bool     # Найдено ли решение
    iters::Int          # Число итераций в методе простой итерации
    x::Vector{T}        # Состав жидкости (мольные доли)
    y::Vector{T}        # Состав газа (мольные доли)
    V::T                # Мольная доля газовой фазы
end

flashsolve_success(iters, x, y, V) = SplitResult(true, Int(iters), x, y, V)
flashsolve_fail(n::Int) = SplitResult(false, NaN, fill(NaN, n), fill(NaN, n), NaN)

rach_rice_eq_one(z::Float64, K::Float64, x::Float64) = z*(1-K)/(x*(K-1)+1)

rach_rice_eq(z::Vector{Float64}, K::Vector{Float64}, x::Float64) = sum([rach_rice_eq_one(z[i], K[i], x) for i in 1:length(z)])


function ptsplit(
    mixture::Mixture, # mixture data
    molfrac,    # known as z_i in Конспект_PTSplit
    pressure,
    RT,
    Kinit;   # initial K{0}
    maxiter=100,
    tolerance=1e-6     # tolerance
)
    K = Kinit   # initial 
    
    for i in 1:maxiter
        try
            # Rashford-Rice solving
            cur_func(arg) = rach_rice_eq(molfrac, K, arg)
            V = bisection(cur_func, 1/(1-maximum(K)), 1/(1-minimum(K)))

            # molfracs in gas and liquid calculation
            x = [molfrac[j]/(1+(K[j]-1)*V) for j in 1:length(molfrac)]
            y = [K[j]*x[j] for j in 1:length(molfrac)]
            
            # fugacities calculation
            ln_φ_gas = log_fugacity_coef(mixture, y, pressure, RT, :gas)
            ln_φ_liq = log_fugacity_coef(mixture, x, pressure, RT, :liq)
            
            # checking
            K_fug = [exp(ln_φ_liq[j] - ln_φ_gas[j]) for j in 1:length(molfrac)]
            diff = maximum([abs(K_fug[j] - K[j]) for j in 1:length(molfrac)])
            if diff < tolerance
                return flashsolve_success(i, x, y, V)
            else
                K = K_fug
            end
        catch e
            @warn e
            return flashsolve_fail(length(molfrac))
        end
    end
    
    return flashsolve_fail(length(molfrac))
end
