/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    laminarPimpleSMOKE-flamelet

Author
    Alberto Cuoci
    Politecnico di Milano
    Department of Chemistry, Materials, and Chemical Engineering
    P.zza Leonardo da Vinci 32, 20133 Milano (Italy)
    email: alberto.cuoci@polimi.it

\*---------------------------------------------------------------------------*/

// This is not a steady state simulation
#define STEADYSTATE 0

// OpenFOAM
#include "fvCFD.H"
#include "pimpleControl.H"
#if OPENFOAM_VERSION >= 40
	#include "fvOptions.H"
	#if DEVVERSION == 1
		#include "pressureControl.H"
	#endif
#else
#include "fvIOoptionList.H"
#endif

// Mixture fraction properties
#include "properties.mixturefraction.H"	

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Returns the density from psi
tmp<volScalarField> Rho(volScalarField& p, volScalarField& psi)
{
	return p*psi;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
	#include "setRootCase.H"
	#include "createTime.H"
	#include "createMesh.H"

	#if OPENFOAM_VERSION >= 40
		#include "createControl.H"
		#include "createTimeControls.H"
		#include "createMRF.H"
	#else
		pimpleControl pimple(mesh);	
	#endif

	#include "readGravitationalAcceleration.H"
	#include "createBasicFields.H"
	#include "readOptions.H"
	#include "createFvOptions.H"
	#include "properties.H"
	#include "createAdditionalFields.H"
	#include "initContinuityErrs.H"

	// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

	Info<< "\nStarting time loop\n" << endl;

	while (runTime.run())
	{
		// Chooses the best time step
		#include "readTimeControls.H"
		#include "compressibleCourantNo.H"
		#include "setDeltaT.H"

		runTime++;
		Info<< "Time = " << runTime.timeName() << nl << endl;

		// Continuity equation
		if (pimple.nCorrPIMPLE() <= 1)
		{
			#include "rhoEqn.H"
		}

		// Pimple loop
		while (pimple.loop())
		{
			#include "UEqn.H"
			#include "zMixEqn.H"

			// --- Pressure corrector loop
			while (pimple.correct())
			{
				#if OPENFOAM_VERSION >= 40
				if (pimple.consistent())
				{
					#include "pcEqn.H"
				}
				else
				#endif
				{
					#include "pEqn.H"
				}
			}

			#include "properties.H"
		}

		// Passive scalars
		#include "tauEqn.H"

		// Write on file
		runTime.write();

		Info 	<< "ExecutionTime = " 	<< runTime.elapsedCpuTime() 	<< " s"
			<< "  ClockTime = " 	<< runTime.elapsedClockTime()	<< " s"
			<< nl << endl;
	}

	Info<< "End\n" << endl;

	return 0;
}

// ************************************************************************* //
