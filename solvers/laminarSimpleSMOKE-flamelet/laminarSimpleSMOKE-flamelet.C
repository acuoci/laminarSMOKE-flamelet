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
    laminarSimpleSMOKE-flamelet

Author
    Alberto Cuoci
    Politecnico di Milano
    Department of Chemistry, Materials, and Chemical Engineering
    P.zza Leonardo da Vinci 32, 20133 Milano (Italy)
    email: alberto.cuoci@polimi.it

\*---------------------------------------------------------------------------*/

// This is a steady state simulation
#define STEADYSTATE 1

// OpenFOAM
#include "fvCFD.H"
#include "simpleControl.H"
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
	#include "readGravitationalAcceleration.H"

	#if OPENFOAM_VERSION >= 40
		#include "createControl.H"
		#include "createMRF.H"
	#else
		simpleControl simple(mesh);	
	#endif

	#include "createBasicFields.H"
	#include "readOptions.H"
	#include "createFvOptions.H"
	#include "properties.H"
	#include "createAdditionalFields.H"
	#include "initContinuityErrs.H"

	dimensionedScalar initialMass = fvc::domainIntegrate(rho);

	// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

	Info<< "\nStarting time loop\n" << endl;

	while (simple.loop())
	{
		Info<< "Time = " << runTime.timeName() << nl << endl;

		// Pressure-velocity SIMPLE corrector
		{
			#include "UEqn.H"
			
			#if OPENFOAM_VERSION >= 40
			if (simple.consistent())
			{
				#include "pcEqn.4x.H"
			}
			else
			#endif
			{
				#include "pEqn.H"
			}

			#include "zMixEqn.H"

			#include "properties.H"
		}	

		// Passive scalars
		#include "tauEqn.H"

		runTime.write();

		Info 	<< "ExecutionTime = " 	<< runTime.elapsedCpuTime() 	<< " s"
			<< "  ClockTime = " 	<< runTime.elapsedClockTime()	<< " s"
			<< nl << endl;
	}

	Info<< "End\n" << endl;

	return 0;
}

// ************************************************************************* //

