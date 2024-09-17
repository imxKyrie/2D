#ifndef COMMON_H
#define COMMON_H

#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <unordered_map>
#include <unordered_set>
#include <vector_types.h>

struct Node
{
    double x, y;
    Node(const double &x, const double &y)
        : x(x), y(y) {}
};

struct Face
{
    unsigned int index;
    unsigned int num;
    unsigned int type;
    std::vector<unsigned int> nodes;
    std::vector<unsigned int> adjacentCells;
    Face(const unsigned int &index, const unsigned int &num, const unsigned int &type, 
         const unsigned int &node1, const unsigned int &node2, const unsigned int &cell1, const unsigned int &cell2)
        : index(index), num(num), type(type), nodes{node1, node2}, adjacentCells{cell1, cell2} {}
};

struct Cell
{
    unsigned int index;
    unsigned int type;
    double2 length;
    double volume;
    double2 centroid;
    std::vector<unsigned int> nodes;
    std::vector<unsigned int> faces;
    Cell(const unsigned int &index, const unsigned int &type, const double &xLength, const double &yLength, const double &volume, 
         const double &xCentroid, const double &yCentroid, const std::vector<unsigned int> &nodes, const std::vector<unsigned int> &faces)
        : index(index), type(type), length{xLength, yLength}, volume(volume), centroid{xCentroid, yCentroid}, nodes(nodes), faces(faces) {}
};

#endif // COMMON_H