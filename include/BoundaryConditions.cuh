#ifndef BOUNDARY_CONDITIONS_CUH
#define BOUNDARY_CONDITIONS_CUH

#include "Common.cuh"

struct SurfaceValues
{
    double pressure;
    double stress;
    double heatFlux;
};

#ifndef INIT_BOUNDARY_CU

extern __constant__ double d_xMeshLength, d_yMeshLength, d_zMeshLength;

extern __constant__ unsigned int d_initParticles;
extern __constant__ double d_Mach, d_CFL;
extern __constant__ double d_wallTem;

extern __constant__ bool d_USP_BGK;
extern __constant__ double d_critFit, d_critFit_FP;

extern __constant__ double d_inflowFnd, d_inflowTem, d_inflowPressure;
extern __constant__ double d_inflowVx, d_inflowVy, d_inflowVz;
extern __constant__ double d_dtm;

extern __constant__ double d_ghostCells;

extern __constant__ unsigned int d_numSurfaces;

extern __device__ SurfaceValues *d_surfacesValues;

extern unsigned int h_initParticles;

#endif

void InitializeBoundary();

#endif // BOUNDARY_CONDITIONS_CUH