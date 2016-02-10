{
	Info<< "Properties evaluation based on mixture fraction... " ;

	double tStart = OpenSMOKE::OpenSMOKEGetCpuTime();

	// Internal fields
	{
		const scalarField& zMixCells = zMix.internalField(); 
		const scalarField& pCells    = p.internalField(); 

		scalarField& rhoCells       =  rho.internalField();
		scalarField& psiCells       =  psi.internalField();
		scalarField& gammazMixCells =  gammazMix.internalField();
		scalarField& lambdaCells    =  lambda.internalField();
		scalarField& muCells        =  mu.internalField();

		forAll(zMixCells, celli)
		{
			// psi = rho/p [kg/m3/Pa]
			rhoCells[celli] = Density(zMixCells[celli]);
			psiCells[celli] = rhoCells[celli]/pCells[celli];

			// Mixture fraction diffusion coefficient [m2/s]
			gammazMixCells[celli] = MixtureFractionDiffusivity(zMixCells[celli]);

			// Thermal conductivity [W/m/K]
			lambdaCells[celli] = ThermalConductivity(zMixCells[celli]);

			// Dynamic viscosity [kg/m/s]
			muCells[celli] = DynamicViscosity(zMixCells[celli]);
		}
    	}

	// Boundaries
	forAll(zMix.boundaryField(), patchi)
	{
		fvPatchScalarField& pp    = p.boundaryField()[patchi];
		fvPatchScalarField& pzmix = zMix.boundaryField()[patchi];
        
		fvPatchScalarField& prho          = rho.boundaryField()[patchi];
		fvPatchScalarField& ppsi          = psi.boundaryField()[patchi];
		fvPatchScalarField& pgammazmix    = gammazMix.boundaryField()[patchi];
		fvPatchScalarField& plambda       = lambda.boundaryField()[patchi];
		fvPatchScalarField& pmu           = mu.boundaryField()[patchi];

		forAll(pzmix, facei)
		{
			// Density [kg/m3]
			prho[facei] = Density(pzmix[facei]);
			ppsi[facei] = prho[facei]/pp[facei];

			// Mixture fraction diffusion coefficient [m2/s]
			pgammazmix[facei] = MixtureFractionDiffusivity(pzmix[facei]);

			// Thermal conductivity [W/m/K]
			plambda[facei] = ThermalConductivity(pzmix[facei]);

			// Dynamic viscosity [kg/m/s]
			pmu[facei] = DynamicViscosity(pzmix[facei]);
		}
	}

	double tEnd = OpenSMOKE::OpenSMOKEGetCpuTime();

	Info << "done in " << tEnd - tStart << " s  (" << (tEnd-tStart)/double(mesh.nCells())*1000. << " ms per cell)" << endl;
}
