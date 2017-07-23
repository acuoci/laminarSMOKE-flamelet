#include "compressibleCreatePhi.H"

#if OPENFOAM_VERSION >= 40 
	#if STEADYSTATE != 1
		pressureControl pressureControl(p, rho, pimple.dict(), false);
	#else
		pressureControl pressureControl(p, rho, simple.dict());
	#endif
	mesh.setFluxRequired(p.name());
#endif


// Kinetic energy
Info<< "Creating field kinetic energy K\n" << endl;
volScalarField K("K", 0.5*magSqr(U));

#if   OPENFOAM_VERSION >= 40

	#include "readhRef.H"
	#include "gh.H"
	#if STEADYSTATE != 1
	#include "readpRef.H"
	#endif

	volScalarField p_rgh
	(
	    IOobject
	    (
		"p_rgh",
		runTime.timeName(),
		mesh,
		IOobject::MUST_READ,
		IOobject::AUTO_WRITE
	    ),
	    mesh
	);

	mesh.setFluxRequired(p_rgh.name());

	// Force p_rgh to be consistent with p
	p_rgh = p - rho*gh;

	#if STEADYSTATE == 1

		label pRefCell = 0;
		scalar pRefValue = 0.0;

		setRefCell
		(
			p,
			p_rgh,
			simple.dict(),
			pRefCell,
			pRefValue
		);

	#endif

#elif   OPENFOAM_VERSION == 30

	#include "readhRef.H"
	#include "gh.H"

	volScalarField p_rgh
	(
	    IOobject
	    (
		"p_rgh",
		runTime.timeName(),
		mesh,
		IOobject::MUST_READ,
		IOobject::AUTO_WRITE
	    ),
	    mesh
	);

	// Force p_rgh to be consistent with p
	p_rgh = p - rho*gh;

	mesh.setFluxRequired(p_rgh.name());

	#if STEADYSTATE == 1

		label pRefCell = 0;
		scalar pRefValue = 0.0;

		setRefCell
		(
			p,
			p_rgh,
			simple.dict(),
			pRefCell,
			pRefValue
		);

	#endif

#else

	Info<< "Calculating field g.h\n" << endl;
	volScalarField gh("gh", g & mesh.C());

	surfaceScalarField ghf("ghf", g & mesh.Cf());

	volScalarField p_rgh
	(
		IOobject
		(
		    "p_rgh",
		    runTime.timeName(),
		    mesh,
		    IOobject::MUST_READ,
		    IOobject::AUTO_WRITE
		),
		mesh
	);

	// Force p_rgh to be consistent with p
	p_rgh = p - rho*gh;

	#if STEADYSTATE == 1

		label pRefCell = 0;
		scalar pRefValue = 0.0;

		setRefCell
		(
			p,
			p_rgh,
			mesh.solutionDict().subDict("SIMPLE"),
			pRefCell,
			pRefValue
		);

	#endif

#endif

// Check the gravity
if ( (mag(g)).value() == 0.)
{
	Info << endl;
	Info << "------------------------------------------------------------------" << endl;
	Info << " ERROR:   gravity is turned off                                   " << endl;
	Info << "          Please consider that this solver is conceived           " << endl;
	Info << "          for flows in which gravity play a major role.           " << endl;
	Info << "          If you think gravity plays a negligible role, consider  " << endl;
	Info << "          the possibility to use laminarPimpleSMOKE-flamelet      " << endl;
	Info << "------------------------------------------------------------------" << endl;
	Info << endl;
	abort();
}
