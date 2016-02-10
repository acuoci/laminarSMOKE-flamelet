// Solve the Momentum equation

tmp<fvVectorMatrix> UEqn
(
    fvm::ddt(rho, U)
  + fvm::div(phi, U)
  - fvm::laplacian(mu, U)
  - fvc::div(mu*dev2(fvc::grad(U)().T()))
 ==
    fvOptions(rho, U)
);

UEqn().relax();

fvOptions.constrain(UEqn());

if (pimple.momentumPredictor())
{
    solve(UEqn() == -fvc::grad(p));

    fvOptions.correct(U);
    K = 0.5*magSqr(U);
}

Info << " Umin: " << min(U).value() << " Umax: " << max(U).value() << endl;
