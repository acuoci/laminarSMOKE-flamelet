Info<< "Reading Solver options\n" << endl;
IOdictionary solverOptions
(
	IOobject
	(
		"solverOptions",
		U.time().constant(),
		U.db(),
		IOobject::MUST_READ,
		IOobject::NO_WRITE
	)
);



// Physical model
Switch zMixEquation;
Switch tauEquation;
scalar thermodynamicPressure = 101325.;

const dictionary& physicalModelDictionary = solverOptions.subDict("PhysicalModel");
{
	tauEquation = Switch(physicalModelDictionary.lookup(word("tauEquation")));
	zMixEquation = Switch(physicalModelDictionary.lookup(word("zMixEquation")));
	thermodynamicPressure = readScalar(physicalModelDictionary.lookup("thermodynamicPressure"));
}

Info<< "Solver options correctly read\n" << endl;
