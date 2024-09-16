#define PARTICLES_MOTION_CU

#include "ParticlesMotion.cuh"
#include "GasConstants.cuh"
#include "BoundaryConditions.cuh"
#include "MacroParams.cuh"
#include "Mesh.cuh"
#include "Particles.cuh"
#include "MathTool.cuh"
#include "FileData.h"

__device__ void SetParticleMacro(unsigned int &cellIdx, double3 &cellVel, double &cellTem) {
    if (d_cellsParams[cellIdx].averageFnd > 0.0) {
        double3 averageVel = d_cellsParams[cellIdx].averageVel;
        double averageTem = d_cellsParams[cellIdx].averageTem;
        double averageFnd = d_cellsParams[cellIdx].averageFnd;
        cellVel.x = averageVel.x / averageFnd;
        cellVel.y = averageVel.y / averageFnd;
        cellVel.z = averageVel.z / averageFnd;
        cellTem = abs((averageTem / averageFnd - 0.5 * (pow(cellVel.x, 2.0) + pow(cellVel.y, 2.0) + pow(cellVel.z, 2.0))) / 1.5 / d_boltz * d_rmass);
    }
    else {
        cellVel.x = 0.0;
        cellVel.y = 0.0;
        cellVel.z = 0.0;
        cellTem = d_inflowTem;
    }
}

__device__ void ParticleCross(double2 &pos_temp, double2 &pos, DeviceNode &node_1, DeviceNode &node_2, bool &cross) {
    double error = 1e-7;
    double2 vec_1 = pos - pos_temp;
    double2 vec_2 = node_2 - node_1;
    double2 vec_3 = node_1 - pos_temp;
    double r = vec_1.x * vec_2.y - vec_1.y * vec_2.x;
    if (r == 0.0) return;
    double s = (vec_3.x * vec_1.y - vec_3.y * vec_1.x) / r;
    if (s < 0.0 + error || s > 1.0 - error) return;
    double t = (vec_3.x * vec_2.y - vec_3.y * vec_2.x) / r;
    if (t < 0.0 + error || t > 1.0 - error) return;
    cross = true;
    pos.x = pos_temp.x + t * vec_1.x;
    pos.y = pos_temp.y + t * vec_1.y;
}

__device__ void CorrectPositionInCell(double2 &pos, unsigned int &cellidx) {
    DeviceCell cell = d_cells[cellidx];
    double error = 0.9999;
    DeviceNode nodes[4];
    unsigned int numFaces = (cell.type == 1) ? 3 : 4;
    DeviceNode centroid = cell.centroid;
    nodes[0] = d_nodes[cell.nodes[0]-1];
    nodes[1] = d_nodes[cell.nodes[1]-1];
    nodes[2] = d_nodes[cell.nodes[2]-1];
    nodes[3] = d_nodes[cell.nodes[3]-1];
    while (true) {
        if (IsParticleWithinCell(pos, nodes, numFaces)) break;
        else {
            pos.x = centroid.x + (pos.x - centroid.x) * error;
            pos.y = centroid.y + (pos.y - centroid.y) * error;
        }
    }
}

__device__ double3 VelocityAfterWallReflection(unsigned int &idx) {
    curandState localstate = d_states[idx];
    double vm = sqrt(2.0 * d_boltz * d_wallTem / d_rmass);
    double3 vel = make_double3(0.0, 0.0, 0.0);
    double u = 0.0, v = 0.0;
    u = sqrt(-log(curand_uniform_double(&localstate)));
    v = 2.0 * d_pi * curand_uniform_double(&localstate);
    vel.x = u * sin(v) * vm;
    u = sqrt(-log(curand_uniform_double(&localstate)));
    vel.y = u * vm;
    u = sqrt(-log(curand_uniform_double(&localstate)));
    v = 2.0 * d_pi * curand_uniform_double(&localstate);
    vel.z = u * sin(v) * vm;
    d_states[idx] = localstate;
    return vel;
}

__device__ void SurfaceSample(double3 &vel_temp, double3 &vel_wall, double2 &face_normal, unsigned int &cellIdx) {
    vel_wall.x = vel_temp.x * face_normal.y - vel_temp.y * face_normal.x;
}

__device__ void BoundaryReflect(double2 &pos_temp, double3 &vel_temp, double3 &averageVel_temp, double &tem_temp, unsigned int &cellIdx, bool &active_temp,
                                double2 &pos, bool &cross_before, int &cell_cross_before, double &dtr, unsigned int &idx) {
    DeviceCell cell = d_cells[cellIdx];
    DeviceNode nodes[4];
    DeviceFace faces[4];
    unsigned int numFaces = (cell.type == 1) ? 3 : 4;
    DeviceNode centroid = cell.centroid;
    nodes[0] = d_nodes[cell.nodes[0]-1];
    nodes[1] = d_nodes[cell.nodes[1]-1];
    nodes[2] = d_nodes[cell.nodes[2]-1];
    nodes[3] = d_nodes[cell.nodes[3]-1];
    faces[0] = d_faces[cell.faces[0]-1];
    faces[1] = d_faces[cell.faces[1]-1];
    faces[2] = d_faces[cell.faces[2]-1];
    faces[3] = d_faces[cell.faces[3]-1];
    bool cross = false;
    unsigned int k_out;
    for (unsigned int i = 0; i < numFaces; i++) {
        ParticleCross(pos_temp, pos, nodes[i], nodes[(i+1)%numFaces], cross);
        if (cross) {
            k_out = i;
            break;
        }
    }
    if (cross) {
        cross_before = cross;
        dtr -= (abs(vel_temp.x) > abs(vel_temp.y)) ?
               (pos.x - pos_temp.x) / vel_temp.x : (pos.y - pos_temp.y) / vel_temp.y;
        if (faces[k_out].type == 2)    // Interior
        {
            cell_cross_before = cellIdx;
            cellIdx = (faces[k_out].adjacentCells[0] - 1 == cellIdx) ?
                      (faces[k_out].adjacentCells[1] - 1) : (faces[k_out].adjacentCells[0] - 1);
            CorrectPositionInCell(pos, cellIdx);
        }
        else if (faces[k_out].type == 3)    // Wall
        {
            double2 face_normal = CalculateNormal(nodes[k_out], nodes[(k_out+1)%numFaces], centroid);
            double3 vel_wall = VelocityAfterWallReflection(idx);
            SurfaceSample(vel_temp, vel_wall, face_normal, cellIdx);
            vel_temp.x = vel_wall.x * face_normal.y + vel_wall.y * face_normal.x;
            vel_temp.y = - vel_wall.x * face_normal.x + vel_wall.y * face_normal.y;
            vel_temp.z = vel_wall.z;
            cell_cross_before = cellIdx;
            averageVel_temp = make_double3(0.0, 0.0, 0.0);
            tem_temp = d_wallTem;
        }
        else if (faces[k_out].type == 4 || faces[k_out].type == 5)    // Inlet or Outlet
        {
            dtr = 0.0;
            active_temp = false;
        }
        else if (faces[k_out].type == 7)    //Symmetry
        {
            double2 face_normal = CalculateNormal(nodes[k_out], nodes[(k_out+1)%numFaces], centroid);
            double2 vel_wall;
            vel_wall.x = vel_temp.x * face_normal.y - vel_temp.y * face_normal.x;
            vel_wall.y = - (vel_temp.x * face_normal.x + vel_temp.y * face_normal.y);
            vel_temp.x = vel_wall.x * face_normal.y + vel_wall.y * face_normal.x;
            vel_temp.y = - vel_wall.x * face_normal.x + vel_wall.y * face_normal.y;
        }
    }
    else dtr = 0.0;
}

__global__ void Motion(Particle *particles) {
    unsigned int idx = blockIdx.x * blockDim.x + threadIdx.x;

    if (idx >= d_numParticles) return;

    Particle &particle = particles[idx];
    double2 pos_temp = particle.pos;
    double3 vel_temp = particle.vel;
    double3 averageVel_temp = particle.averageVel;
    double tem_temp = particle.tem;
    unsigned int cellIdx = particle.cellId - 1;
    bool active_temp = particle.active;

    SetParticleMacro(cellIdx, averageVel_temp, tem_temp);

    double dtr = d_dtm;
    double2 pos;
    bool cross_before = false;
    int cell_cross_before = -1;
    while (true) {
        pos.x = pos_temp.x + vel_temp.x * dtr;
        pos.y = pos_temp.y + vel_temp.y * dtr;
        BoundaryReflect(pos_temp, vel_temp, averageVel_temp, tem_temp, cellIdx, active_temp, pos, cross_before, cell_cross_before, dtr, idx);
        if (dtr <= 0.0) break;
        pos_temp = pos;
    }

    particle.pos = pos;
    particle.vel = vel_temp;
    particle.averageVel = averageVel_temp;
    particle.tem = tem_temp;
    particle.cellId = cellIdx + 1;
    particle.active = active_temp;
}

__device__ void InflowParticles(unsigned int &cellIdx )
{}

void Inflow() {
    IOFace &inletFace = h_inlet[idx];
    unsigned int cellIdx = inletFace.cellId - 1;
    DeviceCell cell = FileData::cells[cellIdx];
    DeviceNode node_1 = inletFace.nodes[0];
    double dh = sqrt(pow(inletFace.nodes[0].x - inletFace.nodes[1].x, 2.0) + pow(inletFace.nodes[0].y - inletFace.nodes[1].y, 2.0));
    double2 face_normal = CalculateNormal(inletFace.nodes[0], inletFace.nodes[1], cell.centroid);
    double vx = d_cellsParams[inletFace.cellId - 1].averageVel.x / d_cellsParams[inletFace.cellId - 1].averageFnd;
}

void ParticlesMotion() {
    unsigned int numParticles = ParticlesManager::GetNumParticles();
    thrust::device_vector<Particle>& d_particles_temp = ParticlesManager::GetParticles();

    Particle *d_particles_ptr = thrust::raw_pointer_cast(d_particles_temp.data());

    unsigned int blockSize = 128;
    unsigned int gridSize = (numParticles + blockSize - 1) / blockSize;

    Motion<<<gridSize, blockSize>>>(d_particles_ptr);
    CUDA_CHECK(cudaDeviceSynchronize());

    unsigned int numInletFaces = FileData::numInletFaces;
    unsigned int numOutletFaces = FileData::numOutletFaces;

    gridSize = (numInletFaces + blockSize - 1) / blockSize;
    Inflow();
    CUDA_CHECK(cudaDeviceSynchronize());
}