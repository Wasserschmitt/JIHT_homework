using CTSolverStability


ch4 = CTSolverStability.VanDerWaalsComponent(0.2303, 0.0000431, 4.5592e6, 190.56 * CTSolverStability.GAS_CONSTANT_SI)
c4h10 = CTSolverStability.VanDerWaalsComponent(1.389, 0.0001164, 3.796e6 , 425.1 * CTSolverStability.GAS_CONSTANT_SI)

mixture = CTSolverStability.VanDerWaalsMixture([ch4, c4h10])

molfracs_1 = range(0, 1;step=0.01)

pressures = range(1e5, 100e5; step = 1e5)

temperature = 294

file = open("npt.txt", "w")

for p in pressures, molfrac_1 in molfracs_1
	RT = temperature * CTSolverStability.GAS_CONSTANT_SI
	Kinit = CTSolverStability.michelsen.(mixture.components, p, RT)
	molfrac = [molfrac_1, 1 - molfrac_1]
	results = [
		CTSolverStability.stability(mixture, molfrac, p, RT, molfrac ./ Kinit, :gas, :liquid)
		CTSolverStability.stability(mixture, molfrac, p, RT, molfrac .* Kinit, :liquid, :gas)
	]
	
	println(
		file,		
		join(Any[
			molfrac_1,
			p,			
			Int(results[1].converged),
			results[1].tpd,
			Int(results[2].converged),
			results[2].tpd
		], '\t')
	)
end

close(file)
