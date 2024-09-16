#define INIT_MESH_CU

#include "Mesh.cuh"
#include "MemoryManager.cuh"
#include "FileData.h"


__constant__ unsigned int d_numNodes;
__constant__ unsigned int d_numFaces;
__constant__ unsigned int d_numInletFaces;
__constant__ unsigned int d_numOutletFaces;
__constant__ unsigned int d_numCells;

__device__ DeviceNode *d_nodes;
__device__ DeviceFace *d_faces;
__device__ IOFace *d_inlet;
__device__ IOFace *d_outlet;
__device__ DeviceCell *d_cells;

IOFace *h_inlet;
IOFace *h_outlet;

__global__ void CalculateCellVolumesAndCentroids() {
    unsigned int idx = blockIdx.x * blockDim.x + threadIdx.x;

    if (idx >= d_numCells) return;

    unsigned int type = d_cells[idx].type;

    if (type == 1) {
        DeviceNode node_1 = d_nodes[d_cells[idx].nodes[0]-1];
        DeviceNode node_2 = d_nodes[d_cells[idx].nodes[1]-1];
        DeviceNode node_3 = d_nodes[d_cells[idx].nodes[2]-1];

        double volume = 0.5 * ((node_1.x - node_2.x) * (node_1.y + node_2.y) +
                               (node_2.x - node_3.x) * (node_2.y + node_3.y) +
                               (node_3.x - node_1.x) * (node_3.y + node_1.y));

        double2 centroid = make_double2((node_1.x + node_2.x + node_3.x) / 3.0, (node_1.y + node_2.y + node_3.y) / 3.0);

        d_cells[idx].volume = volume;
        d_cells[idx].centroid = centroid;
    }
    else if (type == 3) {
        DeviceNode node_1 = d_nodes[d_cells[idx].nodes[0]-1];
        DeviceNode node_2 = d_nodes[d_cells[idx].nodes[1]-1];
        DeviceNode node_3 = d_nodes[d_cells[idx].nodes[2]-1];
        DeviceNode node_4 = d_nodes[d_cells[idx].nodes[3]-1];

        double volume = 0.5 * ((node_1.x - node_3.x) * (node_2.y - node_4.y) +
                               (node_4.x - node_2.x) * (node_1.y - node_3.y));

        double S1 = 0.5 * ((node_1.x - node_2.x) * (node_1.y + node_2.y) +
                           (node_2.x - node_3.x) * (node_2.y + node_3.y) +
                           (node_3.x - node_1.x) * (node_3.y + node_1.y));

        double S2 = 0.5 * ((node_1.x - node_3.x) * (node_1.y + node_3.y) +
                           (node_3.x - node_4.x) * (node_3.y + node_4.y) +
                           (node_4.x - node_1.x) * (node_4.y + node_1.y));

        double C1x = node_1.x + node_2.x + node_3.x;
        double C1y = node_1.y + node_2.y + node_3.y;
        double C2x = node_1.x + node_3.x + node_4.x;
        double C2y = node_1.y + node_3.y + node_4.y;

        double2 centroid = make_double2((C1x * S1 + C2x * S2) / (3.0 * (S1 + S2)), (C1y * S1 + C2y * S2) / (3.0 * (S1 + S2)));

        d_cells[idx].volume = volume;
        d_cells[idx].centroid = centroid;
    }
}

void InitializeMesh() {
    unsigned int h_numNodes = FileData::numNodes;
    unsigned int h_numFaces = FileData::numFaces;
    unsigned int h_numInletFaces = FileData::numInletFaces;
    unsigned int h_numOutletFaces = FileData::numOutletFaces;
    unsigned int h_numCells = FileData::numCells;

    CUDA_CHECK(cudaMemcpyToSymbol(d_numNodes, &h_numNodes, sizeof(unsigned int)));
    CUDA_CHECK(cudaMemcpyToSymbol(d_numFaces, &h_numFaces, sizeof(unsigned int)));
    CUDA_CHECK(cudaMemcpyToSymbol(d_numInletFaces, &h_numInletFaces, sizeof(unsigned int)));
    CUDA_CHECK(cudaMemcpyToSymbol(d_numOutletFaces, &h_numOutletFaces, sizeof(unsigned int)));
    CUDA_CHECK(cudaMemcpyToSymbol(d_numCells, &h_numCells, sizeof(unsigned int)));

    // Transfer mesh's nodes to GPU
    DeviceNode *h_nodes = new DeviceNode[h_numNodes];

    for (unsigned int i = 0; i < h_numNodes; ++i) {
        h_nodes[i].x = FileData::nodes[i].x;
        h_nodes[i].y = FileData::nodes[i].y;
    }

    DeviceNode *d_nodes_temp;
    CUDA_CHECK(cudaMalloc((void**)&d_nodes_temp, h_numNodes * sizeof(DeviceNode)));
    CUDA_CHECK(cudaMemcpy(d_nodes_temp, h_nodes, h_numNodes * sizeof(DeviceNode), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpyToSymbol(d_nodes, &d_nodes_temp, sizeof(DeviceNode*)));

    // Transfer mesh's faces to GPU
    DeviceFace *h_faces = new DeviceFace[h_numFaces];
    CUDA_CHECK(cudaMallocHost((void**)&h_inlet, h_numInletFaces * sizeof(IOFace)));
    CUDA_CHECK(cudaMallocHost((void**)&h_outlet, h_numOutletFaces * sizeof(IOFace)));
    unsigned int count_inlet = 0, count_outlet = 0;

    for (unsigned int i = 0; i < h_numFaces; ++i) {
        h_faces[i].index = FileData::faces[i].index;
        h_faces[i].num = FileData::faces[i].num;
        h_faces[i].type = FileData::faces[i].type;
        h_faces[i].nodes[0] = FileData::faces[i].nodes[0];
        h_faces[i].nodes[1] = FileData::faces[i].nodes[1];
        h_faces[i].adjacentCells[0] = FileData::faces[i].adjacentCells[0];
        h_faces[i].adjacentCells[1] = FileData::faces[i].adjacentCells[1];
        if (h_faces[i].type == 4) {
            h_inlet[count_inlet].cellId = h_faces[i].adjacentCells[0];
            h_inlet[count_inlet].nodes[0].x = h_nodes[h_faces[i].nodes[0]-1].x;
            h_inlet[count_inlet].nodes[0].y = h_nodes[h_faces[i].nodes[0]-1].y;
            h_inlet[count_inlet].nodes[1].x = h_nodes[h_faces[i].nodes[1]-1].x;
            h_inlet[count_inlet].nodes[1].y = h_nodes[h_faces[i].nodes[1]-1].y;
            count_inlet++;
        }
        if (h_faces[i].type == 5) {
            h_outlet[count_outlet].cellId = h_faces[i].adjacentCells[0];
            h_outlet[count_outlet].nodes[0].x = h_nodes[h_faces[i].nodes[0]-1].x;
            h_outlet[count_outlet].nodes[0].y = h_nodes[h_faces[i].nodes[0]-1].y;
            h_outlet[count_outlet].nodes[1].x = h_nodes[h_faces[i].nodes[1]-1].x;
            h_outlet[count_outlet].nodes[1].y = h_nodes[h_faces[i].nodes[1]-1].y;
            count_outlet++;
        }
    }

    DeviceFace *d_faces_temp;
    CUDA_CHECK(cudaMalloc((void**)&d_faces_temp, h_numFaces * sizeof(DeviceFace)));
    CUDA_CHECK(cudaMemcpy(d_faces_temp, h_faces, h_numFaces * sizeof(DeviceFace), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpyToSymbol(d_faces, &d_faces_temp, sizeof(DeviceFace*)));

    IOFace *d_inlet_temp;
    CUDA_CHECK(cudaMalloc((void**)&d_inlet_temp, h_numInletFaces * sizeof(IOFace)));
    CUDA_CHECK(cudaMemcpy(d_inlet_temp, h_inlet, h_numInletFaces * sizeof(IOFace), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpyToSymbol(d_inlet, &d_inlet_temp, sizeof(IOFace*)));

    IOFace* d_outlet_temp;
    CUDA_CHECK(cudaMalloc((void**)&d_outlet_temp, h_numOutletFaces * sizeof(IOFace)));
    CUDA_CHECK(cudaMemcpy(d_outlet_temp, h_outlet, h_numOutletFaces * sizeof(IOFace), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpyToSymbol(d_outlet, &d_outlet_temp, sizeof(IOFace*)));

    // Transfer mesh's cells to GPU
    DeviceCell *h_cells = new DeviceCell[h_numCells];

    for (unsigned int i = 0; i < h_numCells; ++i) {
        h_cells[i].index = FileData::cells[i].index;
        h_cells[i].type = FileData::cells[i].type;
        h_cells[i].length.x = FileData::cells[i].xLength;
        h_cells[i].length.y = FileData::cells[i].yLength;
        h_cells[i].volume = FileData::cells[i].volume;
        h_cells[i].centroid.x = FileData::cells[i].xCentroid;
        h_cells[i].centroid.y = FileData::cells[i].yCentroid;

        if (h_cells[i].type == 1) {
            h_cells[i].nodes[0] = FileData::cells[i].nodes[0];
            h_cells[i].nodes[1] = FileData::cells[i].nodes[1];
            h_cells[i].nodes[2] = FileData::cells[i].nodes[2];
            h_cells[i].nodes[3] = FileData::cells[i].nodes[2];
            h_cells[i].faces[0] = FileData::cells[i].faces[0];
            h_cells[i].faces[1] = FileData::cells[i].faces[1];
            h_cells[i].faces[2] = FileData::cells[i].faces[2];
            h_cells[i].faces[3] = FileData::cells[i].faces[2];
        }
        else if (h_cells[i].type == 3) {
            h_cells[i].nodes[0] = FileData::cells[i].nodes[0];
            h_cells[i].nodes[1] = FileData::cells[i].nodes[1];
            h_cells[i].nodes[2] = FileData::cells[i].nodes[2];
            h_cells[i].nodes[3] = FileData::cells[i].nodes[3];
            h_cells[i].faces[0] = FileData::cells[i].faces[0];
            h_cells[i].faces[1] = FileData::cells[i].faces[1];
            h_cells[i].faces[2] = FileData::cells[i].faces[2];
            h_cells[i].faces[3] = FileData::cells[i].faces[3];
        }
    }

    DeviceCell *d_cells_temp;
    CUDA_CHECK(cudaMalloc((void**)&d_cells_temp, h_numCells * sizeof(DeviceCell)));
    CUDA_CHECK(cudaMemcpy(d_cells_temp, h_cells, h_numCells * sizeof(DeviceCell), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpyToSymbol(d_cells, &d_cells_temp, sizeof(DeviceCell*)));

    // Calculate cells' volumes and Centroids
    unsigned int blockSize = 128;
    unsigned int gridSize = (h_numCells + blockSize - 1) / blockSize;

    dim3 block(blockSize);
    dim3 grid(gridSize);

    CalculateCellVolumesAndCentroids<<<grid, block>>>();
    CUDA_CHECK(cudaDeviceSynchronize());

    delete[] h_nodes;
    delete[] h_faces;
    delete[] h_cells;

    MemoryManager::AddHostPointer(h_inlet);
    MemoryManager::AddHostPointer(h_outlet);
    MemoryManager::AddDevicePointer(d_nodes_temp);
    MemoryManager::AddDevicePointer(d_faces_temp);
    MemoryManager::AddDevicePointer(d_inlet_temp);
    MemoryManager::AddDevicePointer(d_outlet_temp);
    MemoryManager::AddDevicePointer(d_cells_temp);
}