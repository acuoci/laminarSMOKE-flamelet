// Solve the Momentum equation

fvVectorMatrix UEqn
(
  fvm::ddt(rho, U)
+ fvm::div(phi, U)
- fvm::laplacian(mu, U)
- fvc::div(mu*dev2(fvc::grad(U)().T()))
     ==
  fvOptions(rho, U)
);

UEqn.relax();

fvOptions.constrain(UEqn);

if (pimple.momentumPredictor())
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
	K = 0.5*magSqr(U);
}

Info << " Umin: " << min(U).value() << " Umax: " << max(U).value() << endl;
