#ifndef FILE_DATA_H
#define FILE_DATA_H

#include "Common.h"

namespace FileData {
    extern unsigned int numNodes;
    extern unsigned int numFaces;
    extern unsigned int numBoundaryFaces;
    extern unsigned int numInteriorFaces;
    extern unsigned int numCells;
    extern unsigned int numTriCells;
    extern unsigned int numQuadCells;
    extern unsigned int numWallFaces;
    extern unsigned int numInletFaces;
    extern unsigned int numOutletFaces;
    extern unsigned int numSymmetryFaces;
    extern double xMin, xMax, yMin, yMax;    // Computational domain of the mesh
    extern double xCellAverageLength, yCellAverageLength;    // Average length of cells in x and y directions

    extern std::vector<Node> nodes;
    extern std::vector<unsigned int> cellsType;
    extern std::vector<Face> faces;
    extern std::vector<Cell> cells;

    void OutputMeshData();
}

#endif // FILE_DATA_H