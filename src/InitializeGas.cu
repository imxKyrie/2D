#define INIT_GAS_CU

#include "GasConstants.cuh"

// Constants
__constant__ double d_pi;
__constant__ double d_boltz;
__constant__ double d_intError;
__constant__ double d_initCellError;
__constant__ double d_crossError;

// Gas Properties
__constant__ double d_viscosity;
__constant__ double d_refTem;
__constant__ double d_rmass;
__constant__ double d_gamma;

__constant__ double d_initFnd;
__constant__ double d_initTem;

__constant__ double d_PrBGK;

__constant__ double d_omega;
__constant__ double d_refDia;
__constant__ double d_visomk;

__constant__ double d_refCxs;
__constant__ double d_initVm;
__constant__ double d_lambda;
__constant__ double d_refVis;

// Constants
double h_pi = 3.14159265358979323846;
double h_boltz = 1.38064852e-23;

// Gas Properties
double h_rmass = 66.3e-27;
double h_gamma = 5.0/3.0;

double h_initFnd = 2.454e25;
double h_initTem = 298.15;

void InitializeGas() {
    // Constants
    double h_intError = 1.0e-7;
    double h_initCellError = 0.999999;
    double h_crossError = 1.0e-7;

    // Gas Properties
    double h_viscosity = 2.117e-5;
    double h_refTem = 273.0;

    double h_PrBGK = 2.0/3.0;

    double h_omega = 0.81;
    double h_refDia = 4.17e-10;
    double h_visomk = 0.81;
    
    double h_refCxs = h_pi * pow(h_refDia, 2.0);
    double h_initVm = sqrt(2.0 * h_boltz * h_initTem / h_rmass);
    double h_lambda = h_viscosity * 16.0 / (5.0 * sqrt(h_pi) * h_rmass * h_initFnd * h_initVm);
    double h_refVis = 15.0 * sqrt(h_pi * h_rmass * h_boltz * h_refTem) / (2.0 * h_refCxs * (5.0 - 2.0 * h_omega) * (7.0 - 2.0 * h_omega));

    // Transfer constants to GPU
    CUDA_CHECK(cudaMemcpyToSymbol(d_pi, &h_pi, sizeof(double)));
    CUDA_CHECK(cudaMemcpyToSymbol(d_boltz, &h_boltz, sizeof(double)));
    CUDA_CHECK(cudaMemcpyToSymbol(d_intError, &h_intError, sizeof(double)));
    CUDA_CHECK(cudaMemcpyToSymbol(d_initCellError, &h_initCellError, sizeof(double)));
    CUDA_CHECK(cudaMemcpyToSymbol(d_crossError, &h_crossError, sizeof(double)));

    CUDA_CHECK(cudaMemcpyToSymbol(d_viscosity, &h_viscosity, sizeof(double)));
    CUDA_CHECK(cudaMemcpyToSymbol(d_refTem, &h_refTem, sizeof(double)));
    CUDA_CHECK(cudaMemcpyToSymbol(d_rmass, &h_rmass, sizeof(double)));
    CUDA_CHECK(cudaMemcpyToSymbol(d_gamma, &h_gamma, sizeof(double)));

    CUDA_CHECK(cudaMemcpyToSymbol(d_initFnd, &h_initFnd, sizeof(double)));
    CUDA_CHECK(cudaMemcpyToSymbol(d_initTem, &h_initTem, sizeof(double)));

    CUDA_CHECK(cudaMemcpyToSymbol(d_PrBGK, &h_PrBGK, sizeof(double)));

    CUDA_CHECK(cudaMemcpyToSymbol(d_omega, &h_omega, sizeof(double)));
    CUDA_CHECK(cudaMemcpyToSymbol(d_refDia, &h_refDia, sizeof(double)));
    CUDA_CHECK(cudaMemcpyToSymbol(d_visomk, &h_visomk, sizeof(double)));

    CUDA_CHECK(cudaMemcpyToSymbol(d_refCxs, &h_refCxs, sizeof(double)));
    CUDA_CHECK(cudaMemcpyToSymbol(d_initVm, &h_initVm, sizeof(double)));
    CUDA_CHECK(cudaMemcpyToSymbol(d_lambda, &h_lambda, sizeof(double)));
    CUDA_CHECK(cudaMemcpyToSymbol(d_refVis, &h_refVis, sizeof(double)));
}