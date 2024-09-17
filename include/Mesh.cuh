#ifndef MESH_CUH
#define MESH_CUH

#include "Common.cuh"
#include "Common.h"

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
    DeviceCell() 
        : index(0), type(0), length{0, 0}, volume(0), centroid{0, 0},
          nodes{0, 0, 0, 0}, faces{0, 0, 0, 0} {}
    DeviceCell(const Cell& cell)
        : index(cell.index), type(cell.type), length(cell.length), volume(cell.volume), centroid(cell.centroid),
          nodes{cell.nodes[0], cell.nodes[1], cell.nodes[2], cell.nodes[3]}, faces{cell.faces[0], cell.faces[1], cell.faces[2], cell.faces[3]} {}
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