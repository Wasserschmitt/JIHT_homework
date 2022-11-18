struct StabilityResult{T}
	converged::Bool
	molfrac::Vector{T}
	tpd::{T}
end

stability_success(molfrac, tpd) = StabilityResult(true, molfrac, tpd)
stability_fail(n:Int) = StabilityResult(false, fill(NaN, n), NaN)

function stability(
	mixture::Mixture,
	nmol,
	pressure,
	RT,
	Xinit,
	basephase::Symbol,
	testphase::Symbol,
)

	try
		@assert basephase in (:gas, :liquid) "Wrong phase  (:gas, :liquid)"
		@assert testephase in (:gas, :liquid) "Wrong phase  (:gas, :liquid)"
		
		molfrac = nmol/sum(nmol)
		
		#подготовить систему линейных уравнений
		target_function = __create_target(mixture, molfrac, pressire, RT, basephase, testphase)
		
		
		#Начальное приближение
		#Якобиан
		jacobian = ForwardDiff.jacobian(target_function, X)
		#Решить её
		Xstationary = newtonsys(target_function, Xinit, jacobian; maxiter=50, xtol=1e-6, ftol=1e-6)
		#Вернуть ответ
		molfrac_test = Xstationary / sum(Xstationary)
		
		tpd = let x = molfrac_test, z = molfrac
			ln_phi_test = log_fugacity_coeff(mixture, x, pressure, RT, testphase)
			ln_phi_test = log_fugacity_coeff(mixture, z, pressure, RT, basephase)
			delta_pot = log.(x) .+ ln_phi_test .- log.z .- ln_phi_base
			RT * dot(x, delta_pot)
		end
		
		stability_success(molfrac_test, tpd)
		
	catch e
		return stability_fail(ncomponents(mixture))
	end
end



#создание системы уравнений
function __create_target(mixture, molfrac_base, pressure, RT, basephase, testphase)

	base_term = log.(molfrac_base) +. log_fugaсity_coeff(mixture, molfrac_base, pressure, RT, basephase)
	function closure(X)
		molfrac = X / sum(X)
		ln_phi = log_fugaсity_coeff(mixture, molfrac, pressure, RT, testphase)
		return log.(X) .+ ln_phi .- base_term
	end
	return closure
end
