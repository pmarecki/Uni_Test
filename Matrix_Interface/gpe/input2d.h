// %%%%%%%%%%%%%%%% General Preferences %%%%%%%%%%%%%%%%%%%%%%%%%%%
// %%% ALL PARAMETERS ARE IN DIMENSIONLESS UNITS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

//**************** Initial state ****************

#define CALCULATE 0
#define PREVIOUS 1
#define MATLAB 2
#define GROUND 3
#define INITIAL_STATE PREVIOUS

//**************** Time settings ****************

// total time of propagation
//#define DURATION 10.//
// number of time steps
//#define T_STEPS 10000//

#define TOLERANCE 1.e-7//1.e-4
//#define MAXIT 1000000

//**************** Spatial settings ****************

// size of the area in x- and y- directions
//#define LENGTH 10.

// number of spatial points in each direction
// at which the speckle pattern will be created
// must be power of 2!
#define N_POINTS  512//1024//16384//8192//4096//2048     // POWER OF 2  !!

//**************** Atomic parameters ****************

#define SCATTERING_LENGTH 0.1
//#define RINT 0.3
#define ATOM_NUMBER 62
#define DENSITY 0.1

//**************** Potential parameters ****************

#define CORRELATION_LENGTH 1.
//#define VS 0.18
#define VSLIST {32}
//#define RR {42270,16233,15149,738,15282} //RHO=0.1 N=10 NP=128
//#define RR {28834}
//#define RR {32474,2276,15007,24303,28834} //RHO=0.1 N=16 NP=128
//#define RR {32805,8950,32161,47059,23317} //RHO=0.1 N=28 NP=256
//#define RR {5923,1142,12310,11446,24348} //RHO=0.1 N=40 NP=256
//#define RR {3734,9900,27231,45390,24107} //RHO=0.1 N=62 NP=512
#define RR {24107} //RHO=0.1 N=62 NP=512
//#define RR {2613,9908,44183,15037,28946} //RHO=0.1 N=250 NP=1024

#define KS 1.e-5
#define DIRECTION 1