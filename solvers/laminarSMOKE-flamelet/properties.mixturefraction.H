
// Returns the mixture density [kg/m3]
double Density(const double Z)
{
	const double rhoA = 1.;
	const double rhoB = 1.;

	return 1. / ( (1.-Z)/rhoA + Z/rhoB ); 
}

// Returns the diffusion coefficient of the mixture fraction [m2/s]
double MixtureFractionDiffusivity(const double Z)
{
	const double gammaZmix = 1.6e-5;

	return gammaZmix; 
}

// Returns the mixture thermal conductivity [W/m/K]
double ThermalConductivity(const double Z)
{
	const double lambda = 0.02;

	return lambda; 
}

// Returns the mixture dynamic viscosity [kg/m/s]
double DynamicViscosity(const double Z)
{
	const double mu = 1.6e-5;

	return mu; 
}
