#define INIT_CELLS_MACRO_CU

#include "MacroParams.cuh"
#include "GasConstants.cuh"
#include "BoundaryConditions.cuh"
#include "Mesh.cuh"
#include "MemoryManager.cuh"
#include "FileData.h"


__device__ CellParams *d_cellsParams;
__device__ CellSamples *d_cellsSamples;
CellParams *h_cellsParams;
CellSamples *h_cellsSamples;

__global__ void InitializeCellsParams(CellParams *cellsParams) {
    unsigned int idx = blockIdx.x * blockDim.x + threadIdx.x;

    if (idx >= d_numCells) return;

    cellsParams[idx].averageVel.x = d_inflowVx * d_inflowFnd;
    cellsParams[idx].averageVel.y = d_inflowVy * d_inflowFnd;
    cellsParams[idx].averageVel.z = d_inflowVz * d_inflowFnd;
    cellsParams[idx].averageTem = (1.5 * d_boltz * d_inflowTem / d_rmass +
                                       0.5 * (pow(d_inflowVx, 2.0) + pow(d_inflowVy, 2.0) + pow(d_inflowVz, 2.0))) * d_inflowFnd;
    cellsParams[idx].averageFnd = d_inflowFnd;

    cellsParams[idx].vel.x = d_inflowVx * d_inflowFnd;
    cellsParams[idx].vel.y = d_inflowVy * d_inflowFnd;
    cellsParams[idx].vel.z = d_inflowVz * d_inflowFnd;
    cellsParams[idx].tem = (1.5 * d_boltz * d_inflowTem / d_rmass +
                                0.5 * (pow(d_inflowVx, 2.0) + pow(d_inflowVy, 2.0) + pow(d_inflowVz, 2.0))) * d_inflowFnd;
    cellsParams[idx].fnd = d_inflowFnd;

    cellsParams[idx].numParticles = int(d_cells[idx].volume * d_zMeshLength * d_inflowFnd + 0.5 + d_intError);

    cellsParams[idx].pxy[0] = 0.0;
    cellsParams[idx].pxy[1] = 0.0;
    cellsParams[idx].pxy[2] = 0.0;
    cellsParams[idx].pxy[3] = 0.0;
    cellsParams[idx].pxy[4] = 0.0;
    cellsParams[idx].pxy[5] = 0.0;

    cellsParams[idx].q_heat.x = 0.0;
    cellsParams[idx].q_heat.y = 0.0;

    cellsParams[idx].currentTem = 0.0;
    cellsParams[idx].equRange = 0.0;
    cellsParams[idx].equRange_FP = 0.0;
    cellsParams[idx].tau_FP = 0.0;
    cellsParams[idx].equRange_FP_x = 0.0;
    cellsParams[idx].localLambda = 0.0;
    cellsParams[idx].localKn = 0.0;

    cellsParams[idx].ccg1 = 300.0 * d_refCxs * sqrt(d_inflowTem / 300.0);
    cellsParams[idx].ccg2 = 0.0;
    cellsParams[idx].remainder_BGK = 0.0;

    cellsParams[idx].index = idx + 1;

    cellsParams[idx].useDSMC = false;
}

__global__ void InitializeCellsSample(CellSamples *cellsSamples) {
    unsigned int idx = blockIdx.x * blockDim.x + threadIdx.x;

    if (idx >= d_numCells) return;

    cellsSamples[idx].index = idx;
    cellsSamples[idx].vel[0] = 0.0;
    cellsSamples[idx].vel[1] = 0.0;
    cellsSamples[idx].vel[2] = 0.0;
    cellsSamples[idx].tem = 0.0;
    cellsSamples[idx].fnd = 0.0;
    cellsSamples[idx].equ = 0.0;
    cellsSamples[idx].equ_FP = 0.0;
    cellsSamples[idx].separate_regime = 0.0;
    cellsSamples[idx].localLambda = 0.0;
    cellsSamples[idx].localKn = 0.0;
}

void InitializeCellsMacro() {
    CellParams *d_cellsParams_temp;
    CellSamples *d_cellsSamples_temp;

    unsigned int h_numCells = FileData::numCells;

    CUDA_CHECK(cudaMalloc((void**)&d_cellsParams_temp, h_numCells * sizeof(CellParams)));
    CUDA_CHECK(cudaMalloc((void**)&d_cellsSamples_temp, h_numCells * sizeof(CellSamples)));

    unsigned int blockSize = 128;
    unsigned int gridSize = (h_numCells + blockSize - 1) / blockSize;

    dim3 block(blockSize);
    dim3 grid(gridSize);

    InitializeCellsParams<<<grid, block>>>(d_cellsParams_temp);

    InitializeCellsSample<<<grid, block>>>(d_cellsSamples_temp);

    CUDA_CHECK(cudaMallocHost((void**)&h_cellsParams, h_numCells * sizeof(CellParams)));
    CUDA_CHECK(cudaMallocHost((void**)&h_cellsSamples, h_numCells * sizeof(CellSamples)));
    CUDA_CHECK(cudaMemcpy(h_cellsParams, d_cellsParams_temp, h_numCells * sizeof(CellParams), cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(h_cellsSamples, d_cellsSamples_temp, h_numCells * sizeof(CellSamples), cudaMemcpyDeviceToHost));

    CUDA_CHECK(cudaMemcpyToSymbol(d_cellsParams, &d_cellsParams_temp, sizeof(CellParams*)));
    CUDA_CHECK(cudaMemcpyToSymbol(d_cellsSamples, &d_cellsSamples_temp, sizeof(CellSamples*)));

    MemoryManager::AddDevicePointer(d_cellsParams_temp);
    MemoryManager::AddDevicePointer(d_cellsSamples_temp);
}