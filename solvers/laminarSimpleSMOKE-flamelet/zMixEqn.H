{
	if(zMixEquation == true)
	{
		fvScalarMatrix zMixEqn
		(
		    	fvm::div(phi, zMix) 
		        == 
			fvm::laplacian(rho*gammazMix, zMix) +
		     	fvOptions(rho, zMix)
		);

		// Solve
		zMixEqn.relax();
		fvOptions.constrain(zMixEqn);
		zMixEqn.solve();
		fvOptions.correct(zMix);	
	}
}
