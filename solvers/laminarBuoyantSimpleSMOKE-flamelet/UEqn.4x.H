
	MRF.correctBoundaryVelocity(U);

	tmp<fvVectorMatrix> tUEqn
	(
		fvm::div(phi, U)
		+ MRF.DDt(rho, U)
		- fvm::laplacian(mu, U)
      		- fvc::div(mu*dev2(fvc::grad(U)().T()))
		==
		fvOptions(rho, U)
	);
	fvVectorMatrix& UEqn = tUEqn.ref();

	UEqn.relax();

	fvOptions.constrain(UEqn);

	if (simple.momentumPredictor())
	{
		solve
		(
			UEqn
			==
			fvc::reconstruct
			(
				(
					- ghf*fvc::snGrad(rho)
					- fvc::snGrad(p_rgh)
				) * mesh.magSf()
			)
		);

		fvOptions.correct(U);
	}

	Info << " Umin: " << min(U).value() << " Umax: " << max(U).value() << endl;

