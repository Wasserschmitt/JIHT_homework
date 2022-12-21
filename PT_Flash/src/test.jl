using PT_Flash


ch4 = PT_Flash.VanDerWaalsComponent(0.2303, 0.0000431, 4.5592e6, 190.56 * PT_Flash.GAS_CONSTANT_SI)
c4h10 = PT_Flash.VanDerWaalsComponent(1.389, 0.0001164, 3.796e6 , 425.1 * PT_Flash.GAS_CONSTANT_SI)

mixture = PT_Flash.VanDerWaalsMixture([ch4, c4h10])

molfracs_1 = range(0, 1;step=0.01)

pressures = range(1e5, 100e5; step = 1e5)

temperature = 294

file = open("npt.txt", "w")

for p in pressures, molfrac_1 in molfracs_1
	RT = temperature * PT_Flash.GAS_CONSTANT_SI
	Kinit = PT_Flash.michelsen.(mixture.components, p, RT)
	molfrac = [molfrac_1, 1 - molfrac_1]
	results = PT_Flash.ptsplit(mixture, molfrac, p, RT, Kinit)
	println(
		file,		
		join(Any[
			molfrac_1,
			p,			
			join(results.x),
			join(results.y),
			results.V
		], '\t')
	)
end

close(file)
