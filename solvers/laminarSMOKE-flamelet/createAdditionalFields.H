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

