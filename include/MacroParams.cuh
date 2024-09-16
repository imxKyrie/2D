#ifndef MACRO_PARAMS_CUH
#define MACRO_PARAMS_CUH

#include "Common.cuh"

struct CellParams
{
    double3 averageVel;
    double averageTem;
    double averageFnd;

    double3 vel;
    double tem;
    double fnd;
    
    unsigned int numParticles;

    double pxy[6];
    double2 q_heat;

    double currentTem;
    double equRange;
    double equRange_FP;
    double tau_FP;
    double equRange_FP_x;
    double localLambda;
    double localKn;

    double ccg1;
    double ccg2;
    double remainder_BGK;

    unsigned int index;

    bool useDSMC;
};

struct CellSamples
{
    unsigned int index;
    double vel[3];
    double tem;
    double fnd;
    double equ;
    double equ_FP;
    double separate_regime;
    double localLambda;
    double localKn;
};

#ifndef INIT_CELLS_MACRO_CU

extern __device__ CellParams *d_cellsParams;
extern __device__ CellSamples *d_cellsSamples;

extern CellParams *h_cellsParams;
extern CellSamples *h_cellsSamples;

#endif

void InitializeCellsMacro();

#endif // MACRO_PARAMS_CUH