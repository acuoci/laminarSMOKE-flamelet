{
	if(tauEquation == true)
	{
		fvScalarMatrix tauEqn
		(
		    	fvm::div(phi, tau) 
		        == 
                        rho + 
		     	fvOptions(rho, tau)
		);

		// Solve
		tauEqn.relax();
		fvOptions.constrain(tauEqn);
		#if OPENFOAM_VERSION >=60
		tauEqn.solve("tau");
		#else
		tauEqn.solve(mesh.solver("tau"));
		#endif
	}
}
