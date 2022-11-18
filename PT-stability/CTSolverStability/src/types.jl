abstract type Component end

abstract type Mixture end

ncomponents(m:Mixture) = error("Not implemented: ncomponents(::Mixture)")

function log_fugacity_coeff(m::Mixture, molfrac, pressure, RT, phase::Symbol)
	return error("Not implemented: log_fugacity_coeff(::Mixture)")
end
