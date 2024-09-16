#define INIT_PARTICLES_CU

#include "Particles.cuh"
#include "GasConstants.cuh"
#include "BoundaryConditions.cuh"
#include "MacroParams.cuh"
#include "Mesh.cuh"
#include "MathTool.cuh"
#include "MemoryManager.cuh"
#include "FileData.h"

unsigned int ParticlesManager::numParticles;

__device__ unsigned int d_numParticles;

// Define ParticlesManager
thrust::device_vector<Particle>& ParticlesManager::GetParticles() {
    static thrust::device_vector<Particle> d_particles;
    return d_particles;
}

unsigned int ParticlesManager::GetNumParticles() {
    numParticles = GetParticles().size();
    return numParticles;
}

void ParticlesManager::Initialize(unsigned int& numParticles) {
    GetParticles().resize(numParticles);
    ParticlesManager::numParticles = numParticles;
}

void ParticlesManager::Free() {
    GetParticles().clear();
    GetParticles().shrink_to_fit();
}

// Initialize particles
__device__ double2 RandomParticlePosition(unsigned int& idx, unsigned int& cellIdx) {
    curandState localState = d_states[idx];
    double u = curand_uniform_double(&localState);
    double v = curand_uniform_double(&localState);
    if (u + v > 1.0) {
        u = 1.0 - u;
        v = 1.0 - v;
    }

    double2 pos;
    if (d_cells[cellIdx].type == 1) {
        DeviceNode A = d_nodes[d_cells[cellIdx].nodes[0]-1];
        DeviceNode B = d_nodes[d_cells[cellIdx].nodes[1]-1];
        DeviceNode C = d_nodes[d_cells[cellIdx].nodes[2]-1];
        pos.x = A.x + u * (B.x - A.x) + v * (C.x - A.x);
        pos.y = A.y + u * (B.y - A.y) + v * (C.y - A.y);
    }
    else {
        DeviceNode A = d_nodes[d_cells[cellIdx].nodes[0]-1];
        DeviceNode B = d_nodes[d_cells[cellIdx].nodes[1]-1];
        DeviceNode C = d_nodes[d_cells[cellIdx].nodes[2]-1];
        DeviceNode D = d_nodes[d_cells[cellIdx].nodes[3]-1];
        double selector = (TriangleArea(A, B, C)) / (TriangleArea(A, B, C) + TriangleArea(A, C, D));
        double w = curand_uniform_double(&localState);
        if (w < selector)    // Select triangle ABC
        {
            pos.x = A.x + u * (B.x - A.x) + v * (C.x - A.x);
            pos.y = A.y + u * (B.y - A.y) + v * (C.y - A.y);
        }
        else    // Select triangle ACD
        {
            pos.x = A.x + u * (C.x - A.x) + v * (D.x - A.x);
            pos.y = A.y + u * (C.y - A.y) + v * (D.y - A.y);
        }
    }

    d_states[idx] = localState;

    return pos;
}

__device__ double3 RandomParticleVelocity(unsigned int& idx, double3& cellVel, double& cellTem, double& cellFnd) {
    curandState localState = d_states[idx];
    double vm = sqrt(2.0 * d_boltz * cellTem / d_rmass);
    double3 vel = make_double3(0.0, 0.0, 0.0);
    double* vel_ptr = &vel.x;
    double* cellVel_ptr = &cellVel.x;

    double u = 0.0, v = 0.0;

    for (unsigned int i = 0; i < 3; i++) {
        u = sqrt(-log(curand_uniform_double(&localState)));
        v = 2.0 * d_pi * curand_uniform_double(&localState);
        vel_ptr[i] = u * sin(v) * vm + cellVel_ptr[i];
    }

    d_states[idx] = localState;

    return vel;
}

__global__ void RandomParticles(Particle *particles) {
    unsigned int idx = blockIdx.x * blockDim.x + threadIdx.x;

    if (idx >= d_numParticles) return;

    unsigned int cellIdx = idx / d_initParticles;

    double averageFnd = d_cellsParams[cellIdx].averageFnd;

    double3 cellVel_temp = make_double3(d_cellsParams[cellIdx].averageVel.x / averageFnd,
                                        d_cellsParams[cellIdx].averageVel.y / averageFnd,
                                        d_cellsParams[cellIdx].averageVel.z / averageFnd);

    double tem_temp = (d_cellsParams[cellIdx].averageTem/d_cellsParams[cellIdx].averageFnd - 0.5 * (pow(cellVel_temp.x, 2.0) + pow(cellVel_temp.y, 2.0) + pow(cellVel_temp.z, 2.0))) / 1.5 / d_boltz * d_rmass;

    double2 pos = RandomParticlePosition(idx, cellIdx);
    double3 vel = RandomParticleVelocity(idx, cellVel_temp, tem_temp, averageFnd);

    Particle& particle = particles[idx];
    particle.pos = pos;
    particle.vel = vel;
    particle.averageVel = cellVel_temp;
    particle.tem = tem_temp;
    particle.cellId = cellIdx + 1;
    particle.active = true;
}

void InitializeParticles() {
    unsigned int h_numCells = FileData::numCells;
    unsigned int h_numParticles = h_numCells * h_initParticles;

    CUDA_CHECK(cudaMemcpyToSymbol(d_numParticles, &h_numParticles, sizeof(unsigned int)));

    ParticlesManager::Initialize(h_numParticles);
    InitializeRandom(h_numParticles);

    thrust::device_vector<Particle>& d_particles_temp = ParticlesManager::GetParticles();
    Particle *d_particles_ptr = thrust::raw_pointer_cast(d_particles_temp.data());

    unsigned int blockSize = 128;
    unsigned int gridSize = (h_numParticles + blockSize - 1) / blockSize;

    dim3 block(blockSize);
    dim3 grid(gridSize);

    RandomParticles<<<grid, block>>>(d_particles_ptr);
    CUDA_CHECK(cudaDeviceSynchronize());
}