#include "compressibleCreatePhi.H"

#if OPENFOAM_VERSION >= 40 
	#if DEVVERSION == 1
		#if STEADYSTATE != 1
			pressureControl pressureControl(p, rho, pimple.dict(), false);
		#else
			pressureControl pressureControl(p, rho, simple.dict());
		#endif
	#endif
	mesh.setFluxRequired(p.name());
#endif


// Kinetic energy
Info<< "Creating field kinetic energy K\n" << endl;
volScalarField K("K", 0.5*magSqr(U));

// Check the gravity
if ( (mag(g)).value() != 0.)
{
	Info << endl;
	Info << "------------------------------------------------------------------" << endl;
	Info << " ERROR:   gravity is turned on                                    " << endl;
	Info << "          Please consider that this solver is not conceived       " << endl;
	Info << "          for bouyancy driven flows. If you think gravity         " << endl;
	Info << "          plays a non-negligible role, consider the possibility   " << endl;
	Info << "          to use the laminarBouyantPimpleSMOKE-flamelet solver    " << endl;
	Info << "          (unsteady flows) or laminarBouyantSimpleSMOKE-flamelet  " << endl;
	Info << "          solver (for steady-state flows)                         " << endl;
	Info << "------------------------------------------------------------------" << endl;
	Info << endl;
	abort();
}

