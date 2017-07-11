// Solve the Momentum equation

MRF.correctBoundaryVelocity(U);

tmp<fvVectorMatrix> tUEqn
(
	fvm::ddt(rho, U) 
	+ fvm::div(phi, U)
	+ MRF.DDt(rho, U)
	- fvm::laplacian(mu, U)
	- fvc::div(mu*dev2(fvc::grad(U)().T()))
	==
          rho*g
	+ fvOptions(rho, U)
);
fvVectorMatrix& UEqn = tUEqn.ref();

UEqn.relax();

fvOptions.constrain(UEqn);

if (pimple.momentumPredictor())
{
	solve(UEqn == -fvc::grad(p));

	fvOptions.correct(U);
	K = 0.5*magSqr(U);
}

Info << " Umin: " << min(U).value() << " Umax: " << max(U).value() << endl;

