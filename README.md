laminarSMOKE-flamelet
=====================
CFD solver for laminar flows based on OpenFOAM and flamelets

Supported OpenFOAM versions
---------------------------
- OpenFOAM-2.x
- OpenFOAM-4.x
- OpenFOAM-dev

Compulsory libraries
--------------------
- Eigen (http://eigen.tuxfamily.org/index.php?title=Main_Page)
- RapidXML (http://rapidxml.sourceforge.net/)
- Boost C++ (http://www.boost.org/)

Optional libraries
------------------
- Intel MKL (https://software.intel.com/en-us/intel-mkl)

Compilation
-----------
Two different options are available to compile the code
1. Minimalist: no external Intel MKL libraries are needed
2. Minimalist + Intel MKL: better performances thanks to the linking to the Intel MKL libraries is needed


Instructions to compile the Minimalist version
----------------------------------------------
1. Open the `mybashrc.minimalist`, select the OpenFOAM version you are working with and update the required paths
2. Type: `source mybashrc.minimalist`
3. Compile the steady-state solver: from the `solver/laminarSimpleSMOKE-flamelet` folder type `wclean` and then `wmake`
4. Compile the unsteady solver: from the `solver/laminarPimpleSMOKE-flamelet` folder type `wclean` and then `wmake`


Instructions to compile the Minimalist+MKL version
--------------------------------------------------
1. Open the `mybashrc.minimalist.mkl`select the OpenFOAM version you are working with and update the required paths (in particular those pointing to the Intel MKL libraries)
2. Type: `source mybashrc.minimalist.mkl`
3. Compile the steady-state solver: from the `solver/laminarSimpleSMOKE-flamelet` folder type `wclean` and then `wmake`
4. Compile the unsteady solver: from the `solver/laminarPimpleSMOKE-flamelet` folder type `wclean` and then `wmake`


Run your first case
-------------------
The folder `cases/CounterFlow/Coarse` contains a simple test case (counterflow diffusion flame).

1. Unsteady simulation: Open the `laminarPimpleSMOKE-flamelet` folder, build the mesh using the `blockMesh` utility, and run the case using the `laminarPimpleSMOKE-flamelet` solver. Even if you are interested in steady state conditions, we strongly suggest to always start with unsteady calculations to create a reasonable first-guess solution for the application of the steady state solver.

2. Steady state simulation: you can now move to the `laminarSimpleSMOKE-flamelet` folder. Copy the last time folder calculated by the unsteady solver (point 1 above), build the mesh using the `blockMesh` utility, and run the case using the `laminarSimpleSMOKE-flamelet` solver. In order to reach the steady state conditions, 5000-6000 iterations are enough.
