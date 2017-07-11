   // Solve the Momentum equation

MRF.correctBoundaryVelocity(U);

tmp<fvVectorMatrix> tUEqn
(
	  fvm::div(phi, U)
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

solve(UEqn == -fvc::grad(p));

fvOptions.correct(U);

Info << " Umin: " << min(U).value() << " Umax: " << max(U).value() << endl;

