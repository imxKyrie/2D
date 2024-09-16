#ifndef GAS_CONSTANTS_CUH
#define GAS_CONSTANTS_CUH

#include "Common.cuh"

#ifndef INIT_GAS_CU

// Constants
extern __constant__ double d_pi;
extern __constant__ double d_boltz;
extern __constant__ double d_intError;
extern __constant__ double d_initCellError;
extern __constant__ double d_crossError;

// Gas Properties
extern __constant__ double d_viscosity;
extern __constant__ double d_refTem;
extern __constant__ double d_rmass;
extern __constant__ double d_gamma;

extern __constant__ double d_initFnd;
extern __constant__ double d_initTem;

extern __constant__ double d_PrBGK;

extern __constant__ double d_omega;
extern __constant__ double d_refDia;
extern __constant__ double d_visomk;

extern __constant__ double d_refCxs;
extern __constant__ double d_initVm;
extern __constant__ double d_lambda;
extern __constant__ double d_refVis;

// Constants
extern double h_pi;
extern double h_boltz;

// Gas Properties
extern double h_rmass;
extern double h_gamma;

extern double h_initFnd;
extern double h_initTem;

#endif

void InitializeGas();

void verifyConstantsOnHost();

#endif // GAS_CONSTANTS_CUH