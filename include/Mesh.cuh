#ifndef MESH_CUH
#define MESH_CUH

#include "Common.cuh"

using DeviceNode = double2;

struct DeviceFace
{
    unsigned int index;
    unsigned int num;
    unsigned int type;
    unsigned int nodes[2];
    unsigned int adjacentCells[2];
};

struct IOFace
{
    unsigned int cellId;
    DeviceNode nodes[2];
};

struct DeviceCell
{
    unsigned int index;
    unsigned int type;
    double2 length;
    double volume;
    double2 centroid;
    unsigned int nodes[4];
    unsigned int faces[4];
};

#ifndef INIT_MESH_CU

extern __constant__ unsigned int d_numNodes;
extern __constant__ unsigned int d_numFaces;
extern __constant__ unsigned int d_numInletFaces;
extern __constant__ unsigned int d_numOutletFaces;
extern __constant__ unsigned int d_numCells;

extern __device__ DeviceNode *d_nodes;
extern __device__ DeviceFace *d_faces;
extern __device__ IOFace *d_inlet;
extern __device__ IOFace *d_outlet;
extern __device__ DeviceCell *d_cells;

extern IOFace *h_inlet;
extern IOFace *h_outlet;

#endif

void InitializeMesh();

#endif // MESH_CUH