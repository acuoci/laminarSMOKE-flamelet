{
    volScalarField rAU(1.0/UEqn().A());
    volVectorField HbyA("HbyA", U);
    HbyA = rAU*UEqn().H();

    UEqn.clear();

    bool closedVolume = false;

    if (simple.transonic())
    {
        surfaceScalarField phid
        (
            "phid",
            fvc::interpolate(psi)*(fvc::interpolate(HbyA) & mesh.Sf())
        );

        fvOptions.relativeFlux(fvc::interpolate(psi), phid);

        while (simple.correctNonOrthogonal())
        {
            fvScalarMatrix pEqn
            (
                fvm::div(phid, p)
              - fvm::laplacian(rho*rAU, p)
              ==
                fvOptions(psi, p, rho.name())
            );

            // Relax the pressure equation to ensure diagonal-dominance
            pEqn.relax();

            fvOptions.constrain(pEqn);

            pEqn.setReference(pRefCell, pRefValue);

            pEqn.solve();

            if (simple.finalNonOrthogonalIter())
            {
                phi == pEqn.flux();
            }
        }
    }
    else
    {
        surfaceScalarField phiHbyA
        (
            "phiHbyA",
            fvc::interpolate(rho)*(fvc::interpolate(HbyA) & mesh.Sf())
        );

        fvOptions.relativeFlux(fvc::interpolate(rho), phiHbyA);

        closedVolume = adjustPhi(phiHbyA, U, p);

        while (simple.correctNonOrthogonal())
        {
            fvScalarMatrix pEqn
            (
                fvc::div(phiHbyA)
              - fvm::laplacian(rho*rAU, p)
              ==
                fvOptions(psi, p, rho.name())
            );

            pEqn.setReference(pRefCell, pRefValue);

            fvOptions.constrain(pEqn);

            pEqn.solve();

            if (simple.finalNonOrthogonalIter())
            {
                phi = phiHbyA + pEqn.flux();
            }
        }
    }

    // Continuity errors
    {
	    volScalarField contErr(fvc::div(phi));

	    scalar sumLocalContErr = runTime.deltaTValue()*
		mag(contErr)().weightedAverage(mesh.V()).value();

	    scalar globalContErr = runTime.deltaTValue()*
		contErr.weightedAverage(mesh.V()).value();
	    cumulativeContErr += globalContErr;

	    Info<< "time step continuity errors : sum local = " << sumLocalContErr
		<< ", global = " << globalContErr
		<< ", cumulative = " << cumulativeContErr
		<< endl;
    }

    // Explicitly relax pressure for momentum corrector
    p.relax();

    U = HbyA - rAU*fvc::grad(p);
    U.correctBoundaryConditions();
    fvOptions.correct(U);

    // For closed-volume cases adjust the pressure and density levels
    // to obey overall mass continuity
    if (closedVolume)
    {
        p += (initialMass - fvc::domainIntegrate(psi*p))
            /fvc::domainIntegrate(psi);
    }

    rho = Rho(p,psi);
    rho = max(rho, rhoMin);
    rho = min(rho, rhoMax);

    if (!simple.transonic())
    {
        rho.relax();
    }

    Info<< "rho max/min : "
        << max(rho).value() << " "
        << min(rho).value() << endl;

    // Be careful
    // In case of energy balance equation here you can calculate the DpDt contribution
    // The present solver does not need it
}

