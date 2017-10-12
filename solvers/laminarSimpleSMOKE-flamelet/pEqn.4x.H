{
	volScalarField rAU(1.0/UEqn.A());
	surfaceScalarField rhorAUf("rhorAUf", fvc::interpolate(rho*rAU));
	volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p));
	tUEqn.clear();

	bool closedVolume = false;

	surfaceScalarField phiHbyA("phiHbyA", fvc::interpolate(rho)*fvc::flux(HbyA));
	MRF.makeRelative(fvc::interpolate(rho), phiHbyA);

	// Update the pressure BCs to ensure flux consistency
	constrainPressure(p, rho, U, phiHbyA, rhorAUf, MRF);

	if (simple.transonic())
	{
	    surfaceScalarField phid
	    (
		"phid",
		(fvc::interpolate(psi)/fvc::interpolate(rho))*phiHbyA
	    );

	    phiHbyA -= fvc::interpolate(psi*p)*phiHbyA/fvc::interpolate(rho);

	    while (simple.correctNonOrthogonal())
	    {
		fvScalarMatrix pEqn
		(
		    fvc::div(phiHbyA)
		  + fvm::div(phid, p)
		  - fvm::laplacian(rhorAUf, p)
		  ==
		    fvOptions(psi, p, rho.name())
		);

		// Relax the pressure equation to ensure diagonal-dominance
		pEqn.relax();

		#if DEVVERSION == 1
		pEqn.setReference
		(
		    pressureControl.refCell(),
		    pressureControl.refValue()
		);
		#else
		pEqn.setReference(pRefCell, pRefValue);
		#endif

		pEqn.solve();

		if (simple.finalNonOrthogonalIter())
		{
		    phi = phiHbyA + pEqn.flux();
		}
	    }
	}
	else
	{
	    closedVolume = adjustPhi(phiHbyA, U, p);

	    while (simple.correctNonOrthogonal())
	    {
		fvScalarMatrix pEqn
		(
		    fvc::div(phiHbyA)
		  - fvm::laplacian(rhorAUf, p)
		  ==
		    fvOptions(psi, p, rho.name())
		);

		#if DEVVERSION == 1
		pEqn.setReference
		(
		    pressureControl.refCell(),
		    pressureControl.refValue()
		);
		#else
		pEqn.setReference(pRefCell, pRefValue);
		#endif

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

		scalar sumLocalContErr = runTime.deltaTValue()*mag(contErr)().weightedAverage(mesh.V()).value();

		scalar globalContErr = runTime.deltaTValue()*contErr.weightedAverage(mesh.V()).value();
		cumulativeContErr += globalContErr;

		Info	<< "time step continuity errors : sum local = " << sumLocalContErr
			<< ", global = " << globalContErr
			<< ", cumulative = " << cumulativeContErr
			<< endl;
	}

	// Explicitly relax pressure for momentum corrector
	p.relax();

	U = HbyA - rAU*fvc::grad(p);
	U.correctBoundaryConditions();
	fvOptions.correct(U);

	#if DEVVERSION == 1
	bool pLimited = pressureControl.limit(p);
	#else
	bool pLimited = false;
	#endif

	// For closed-volume cases adjust the pressure and density levels
	// to obey overall mass continuity
	if (closedVolume)
	{
	    p += (initialMass - fvc::domainIntegrate(psi*p))
		/fvc::domainIntegrate(psi);
	}

	if (pLimited || closedVolume)
	{
	    p.correctBoundaryConditions();
	}

	rho = Rho(p,psi);

	if (!simple.transonic())
	{
	    rho.relax();
	}
	Info<< "rho max/min : " << max(rho).value() << " " << min(rho).value() << endl;

	// Be careful
	// In case of energy balance equation here you can calculate the DpDt contribution
	// The present solver does not need it
}
