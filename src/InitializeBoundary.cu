#define INIT_BOUNDARY_CU

#include "BoundaryConditions.cuh"
#include "GasConstants.cuh"
#include "MemoryManager.cuh"
#include "FileData.h"

__constant__ double d_xMeshLength, d_yMeshLength, d_zMeshLength;

__constant__ unsigned int d_initParticles;
__constant__ double d_Mach, d_CFL;
__constant__ double d_wallTem;

__constant__ bool d_USP_BGK;
__constant__ double d_critFit, d_critFit_FP;

__constant__ double d_inflowFnd, d_inflowTem, d_inflowPressure;
__constant__ double d_inflowVx, d_inflowVy, d_inflowVz;
__constant__ double d_dtm;

__constant__ double d_ghostCells;

__constant__ unsigned int d_numSurfaces;

__device__ SurfaceValues *d_surfacesValues;

unsigned int h_initParticles = 500;

__global__ void SetInitialSurfaceValues(SurfaceValues *surfacesValues) {
    unsigned int idx = blockIdx.x * blockDim.x + threadIdx.x;

    if (idx >= d_numSurfaces) return;

    surfacesValues[idx].pressure = 0.0;
    surfacesValues[idx].stress = 0.0;
    surfacesValues[idx].heatFlux = 0.0;
}

void InitializeBoundary() {
    // Initialize Constants
    double h_Mach = 0.0, h_CFL = 0.005;
    double h_wallTem = 298.15;

    bool h_USP_BGK = true;
    double h_critFit = 0.1, h_critFit_FP = 1.5;

    double h_xMeshLength = FileData::xMax - FileData::xMin;
    double h_yMeshLength = FileData::yMax - FileData::yMin;
    double h_zMeshLength = h_initParticles / (h_initFnd * FileData::xCellAverageLength * FileData::yCellAverageLength);

    double h_inflowFnd = h_initFnd;
    double h_inflowTem = h_initTem;
    double h_inflowPressure = h_initFnd * h_boltz * h_initTem;

    double h_inflowVx = h_Mach * sqrt(h_gamma * h_boltz * h_inflowTem / h_rmass);
    double h_inflowVy = 0.0;
    double h_inflowVz = 0.0;

    double h_initVm = sqrt(2.0 * h_boltz * h_initTem / h_rmass);
    double h_dtm = FileData::xCellAverageLength / (h_inflowVx + h_initVm) * h_CFL;

    double h_ghostCells = 3.0 * FileData::xCellAverageLength;

    unsigned int h_numSurfaces = FileData::numBoundaryFaces;

    CUDA_CHECK(cudaMemcpyToSymbol(d_initParticles, &h_initParticles, sizeof(int)));
    CUDA_CHECK(cudaMemcpyToSymbol(d_Mach, &h_Mach, sizeof(double)));
    CUDA_CHECK(cudaMemcpyToSymbol(d_CFL, &h_CFL, sizeof(double)));
    CUDA_CHECK(cudaMemcpyToSymbol(d_wallTem, &h_wallTem, sizeof(double)));

    CUDA_CHECK(cudaMemcpyToSymbol(d_USP_BGK, &h_USP_BGK, sizeof(bool)));
    CUDA_CHECK(cudaMemcpyToSymbol(d_critFit, &h_critFit, sizeof(double)));
    CUDA_CHECK(cudaMemcpyToSymbol(d_critFit_FP, &h_critFit_FP, sizeof(double)));

    CUDA_CHECK(cudaMemcpyToSymbol(d_xMeshLength, &h_xMeshLength, sizeof(double)));
    CUDA_CHECK(cudaMemcpyToSymbol(d_yMeshLength, &h_yMeshLength, sizeof(double)));
    CUDA_CHECK(cudaMemcpyToSymbol(d_zMeshLength, &h_zMeshLength, sizeof(double)));

    CUDA_CHECK(cudaMemcpyToSymbol(d_inflowFnd, &h_inflowFnd, sizeof(double)));
    CUDA_CHECK(cudaMemcpyToSymbol(d_inflowTem, &h_inflowTem, sizeof(double)));
    CUDA_CHECK(cudaMemcpyToSymbol(d_inflowPressure, &h_inflowPressure, sizeof(double)));

    CUDA_CHECK(cudaMemcpyToSymbol(d_inflowVx, &h_inflowVx, sizeof(double)));
    CUDA_CHECK(cudaMemcpyToSymbol(d_inflowVy, &h_inflowVy, sizeof(double)));
    CUDA_CHECK(cudaMemcpyToSymbol(d_inflowVz, &h_inflowVz, sizeof(double)));

    CUDA_CHECK(cudaMemcpyToSymbol(d_dtm, &h_dtm, sizeof(double)));

    CUDA_CHECK(cudaMemcpyToSymbol(d_ghostCells, &h_ghostCells, sizeof(double)));

    CUDA_CHECK(cudaMemcpyToSymbol(d_numSurfaces, &h_numSurfaces, sizeof(int)));

    // Initialize d_surfaceValues
    SurfaceValues *d_surfacesValues_temp;
    CUDA_CHECK(cudaMalloc((void**)&d_surfacesValues_temp, h_numSurfaces * sizeof(SurfaceValues)));

    unsigned int blockSize = 128;
    unsigned int gridSize = (h_numSurfaces + blockSize - 1) / blockSize;

    dim3 block(blockSize);
    dim3 grid(gridSize);

    SetInitialSurfaceValues<<<grid, block>>>(d_surfacesValues_temp);
    CUDA_CHECK(cudaDeviceSynchronize());

    CUDA_CHECK(cudaMemcpyToSymbol(d_surfacesValues, &d_surfacesValues_temp, sizeof(SurfaceValues*)));

    MemoryManager::AddDevicePointer(d_surfacesValues_temp);
}