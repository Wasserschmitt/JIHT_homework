using CTSolverStability

c1 = VanDerWaals(2.283 /10, 0.04278/1000, 4. 5592e6,190.56 * CTSolverStability.GAS_CONSTANT_SI)
c4 = VanDerWaals(14.66/10, 0.1226/1000, 3.796e6 , 425.1*CTSolverStability.GAS_CONSTANT_SI)

mixture = VanDerWaalsMixtire([c1, c4])

molfrac = [0.8, 0.2]

pressures = range(1e5, 100e5; length = 5)

temperatures = range(150, 500; length = 5)

for p in pressures, t in temperatures
	RT = t * CTSolverStability.GAS_CONSTANT_SI
	Kinit = CTSolverStability.michelsen.(mixture.components, p, RT)
	
	results = [
		stability(mixture, molfrac, p, RT, molfrac ./ Kinit, :gas, :liquid)
		stability(mixture, molfrac, p, RT, molfrac .* Kinit, :liquid, :gas)
	]
	
	println(
		io,
		join(Any[
			Int(results[1].converged),
			results[1].tpd,
			Int(results[2].converged),
			results[2].tpd,
		], '\t')
	)
end
