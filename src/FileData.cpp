#include "FileData.h"
#include <regex>

namespace FileData
{
    unsigned int numNodes = 0;
    unsigned int numFaces = 0;
    unsigned int numBoundaryFaces = 0;
    unsigned int numInteriorFaces = 0;
    unsigned int numCells = 0;
    unsigned int numTriCells = 0;
    unsigned int numQuadCells = 0;
    unsigned int numWallFaces = 0, numInletFaces = 0, numOutletFaces = 0, numSymmetryFaces = 0;
    double xMin = 0.0, xMax = 0.0, yMin = 0.0, yMax = 0.0;
    double xCellAverageLength = 0.0, yCellAverageLength = 0.0;

    std::vector<Node> nodes;
    std::vector<unsigned int> cellsType;
    std::vector<Face> faces;
    std::vector<Cell> cells;

    bool OpenFileWithPrompt(std::ifstream &file, std::string &filename)
    {
        while (true)
        {
            file.open(filename);
            if (file.is_open()) return true;
            else
            {
                std::cerr << "Fail to open the file, please ensure that the file exists and the filename is correct!" << std::endl;
                std::cout << "Do you want to exit the program? (y/n): ";
                std::string choice;
                std::cin >> choice;
                if (choice == "y" || choice == "Y") return false;
                else
                {
                    std::cout << "Please enter the filename again: ";
                    std::cin >> filename;
                }
            }
        }
    }

    unsigned int JudgeFaceType(const std::string &line)
    {
        std::regex re("\\d+");
        std::sregex_iterator it(line.begin(), line.end(), re);
        std::sregex_iterator end;
        std::vector<unsigned int> nums;
        while (it!= end)
        {
            nums.push_back(std::stoi(it->str()));
            ++it;
        }

        if (nums.size() >= 2) return nums[nums.size() - 2];
        else return -1;
    }

    void ProcessFile(std::string &filename,
                     std::vector<Node> &nodes,
                     std::vector<unsigned int> &cellsType,
                     std::vector<Face> &faces,
                     unsigned int &numNodes,
                     unsigned int &numFaces,
                     unsigned int &numBoundaryFaces,
                     unsigned int &numInteriorFaces,
                     unsigned int &numCells,
                     unsigned int &numTriCells,
                     unsigned int &numQuadCells,
                     unsigned int &numWallFaces,
                     unsigned int &numInletFaces,
                     unsigned int &numOutletFaces,
                     unsigned int &numSymmetryFaces)
    {
        std::ifstream file;

        if (!OpenFileWithPrompt(file, filename))
        {
            std::cerr << "Exiting program." << std::endl;
            std::exit(1);
        }

        std::cout << "Successfully found the file, now start processing..." << std::endl;

        std::string line;
        bool captureNodes, capturecellsType, capturenumInteriorFaces, captureWallFaces,
             captureInletFaces, captureOutletFaces, captureSymmetryFaces;

        while (getline(file, line))
        {
            if (line.find("(0 \"Number of Nodes") != std::string::npos)
            {
                unsigned int start = line.find_first_of("0123456789", line.find("Number of Nodes"));
                if (start != std::string::npos)
                {
                    unsigned int end = line.find_first_not_of("0123456789", start);
                    std::string numberstr = line.substr(start, end - start);
                    std::istringstream(numberstr) >> numNodes;
                    nodes.reserve(numNodes);
                    std::cout << "Number of Nodes: " << numNodes << std::endl;
                }
                continue;
            }

            if (line.find("(0 \"Total Number of Faces") != std::string::npos)
            {
                unsigned int start = line.find_first_of("0123456789", line.find("Total Number of Faces"));
                if (start != std::string::npos)
                {
                    unsigned int end = line.find_first_not_of("0123456789", start);
                    std::string numberstr = line.substr(start, end - start);
                    std::istringstream(numberstr) >> numFaces;
                    std::cout << "Number of Faces: " << numFaces <<std::endl;
                }
                continue;
            }

            if (line.find("Boundary Faces") != std::string::npos)
            {
                unsigned int start = line.find_first_of("0123456789", line.find("Boundary Faces"));
                if (start != std::string::npos)
                {
                    unsigned int end = line.find_first_not_of("0123456789", start);
                    std::string numberstr = line.substr(start, end - start);
                    std::istringstream(numberstr) >> numBoundaryFaces;
                    std::cout << "Boundary Faces: " << numBoundaryFaces <<std::endl;
                }
                continue;
            }

            if (line.find("Interior Faces") != std::string::npos)
            {
                unsigned int start = line.find_first_of("0123456789", line.find("Interior Faces"));
                if (start != std::string::npos)
                {
                    unsigned int end = line.find_first_not_of("0123456789", start);
                    std::string numberstr = line.substr(start, end - start);
                    std::istringstream(numberstr) >> numInteriorFaces;
                    std::cout << "Interior Faces: " << numInteriorFaces <<std::endl;
                }
                continue;
            }

            if (line.find("Total Number of Cells") != std::string::npos)
            {
                unsigned int start = line.find_first_of("0123456789", line.find("Total Number of Cells"));
                if (start != std::string::npos)
                {
                    unsigned int end = line.find_first_not_of("0123456789", start);
                    std::string numberstr = line.substr(start, end - start);
                    std::istringstream(numberstr) >> numCells;
                    cellsType.reserve(numCells);
                    std::cout << "Number of Cells: " << numCells <<std::endl;
                }
                continue;
            }

            if (line.find("Tri cells") != std::string::npos)
            {
                unsigned int start = line.find_first_of("0123456789", line.find("Tri cells"));
                if (start != std::string::npos)
                {
                    unsigned int end = line.find_first_not_of("0123456789", start);
                    std::string numberstr = line.substr(start, end - start);
                    std::istringstream(numberstr) >> numTriCells;
                    std::cout << "Tri cells: " << numTriCells <<std::endl;
                }
                continue;
            }

            if (line.find("Quad cells") != std::string::npos)
            {
                unsigned int start = line.find_first_of("0123456789", line.find("Quad cells"));
                if (start != std::string::npos)
                {
                    unsigned int end = line.find_first_not_of("0123456789", start);
                    std::string numberstr = line.substr(start, end - start);
                    std::istringstream(numberstr) >> numQuadCells;
                    std::cout << "Quad cells: " << numQuadCells <<std::endl;
                }
                continue;
            }

            if (line.find("(10 (1 1") != std::string::npos)
            {
                captureNodes = true;
                continue;
            }

            if (captureNodes)
            {
                std::istringstream ss(line);
                double x, y;
                while(ss >> x >> y) nodes.emplace_back(x, y);
                if (nodes.size() >= numNodes) captureNodes = false;
                continue;
            }

            if (line.find("(12 (2 1") !=std::string::npos)
            {
                capturecellsType = true;
                continue;
            }

            if (capturecellsType)
            {
                std::istringstream ss(line);
                unsigned int type;
                while (ss >> type) cellsType.emplace_back(type);
                if (cellsType.size() >= numCells) capturecellsType = false;
                continue;
            }

            if (line.find("(13 (") != std::string::npos && line.find("(13 (0") == std::string::npos)
            {
                if (JudgeFaceType(line) == 2)
                {
                    capturenumInteriorFaces = true;
                    captureWallFaces = false;
                    captureInletFaces = false;
                    captureOutletFaces = false;
                    captureSymmetryFaces = false;
                }
                else if (JudgeFaceType(line) == 3)
                {
                    capturenumInteriorFaces = false;
                    captureWallFaces = true;
                    captureInletFaces = false;
                    captureOutletFaces = false;
                    captureSymmetryFaces = false;
                }
                else if (JudgeFaceType(line) == 4)
                {
                    capturenumInteriorFaces = false;
                    captureWallFaces = false;
                    captureInletFaces = true;
                    captureOutletFaces = false;
                    captureSymmetryFaces = false;
                }
                else if (JudgeFaceType(line) == 5)
                {
                    capturenumInteriorFaces = false;
                    captureWallFaces = false;
                    captureInletFaces = false;
                    captureOutletFaces = true;
                    captureSymmetryFaces = false;
                }
                else if (JudgeFaceType(line) == 7)
                {
                    capturenumInteriorFaces = false;
                    captureWallFaces = false;
                    captureInletFaces = false;
                    captureOutletFaces = false;
                    captureSymmetryFaces = true;
                }
                continue;
            }

            if (capturenumInteriorFaces)
            {
                std::istringstream ss(line);
                unsigned int i = faces.size() + 1;
                unsigned int num;
                std::string node1, node2, cell1, cell2;
                while (ss >> num >> node1 >> node2 >> cell1 >> cell2)
                {
                    faces.emplace_back(i, num, 2, std::stoi(node1, nullptr, 16), std::stoi(node2, nullptr, 16), std::stoi(cell1, nullptr, 16), std::stoi(cell2, nullptr, 16));
                }
                if (faces.size() >= numInteriorFaces) capturenumInteriorFaces = false;
                continue;
            }

            if (captureWallFaces)
            {
                std::istringstream ss(line);
                unsigned int i = faces.size() + 1;
                unsigned int num;
                std::string node1, node2, cell1, cell2;
                while (ss >> num >> node1 >> node2 >> cell1 >> cell2)
                {
                    faces.emplace_back(i, num, 3, std::stoi(node1, nullptr, 16), std::stoi(node2, nullptr, 16), std::stoi(cell1, nullptr, 16), std::stoi(cell2, nullptr, 16));
                    numWallFaces++;
                }
                continue;
            }

            if (captureInletFaces)
            {
                std::istringstream ss(line);
                unsigned int i = faces.size() + 1;
                unsigned int num;
                std::string node1, node2, cell1, cell2;
                while (ss >> num >> node1 >> node2 >> cell1 >> cell2)
                {
                    faces.emplace_back(i, num, 4, std::stoi(node1, nullptr, 16), std::stoi(node2, nullptr, 16), std::stoi(cell1, nullptr, 16), std::stoi(cell2, nullptr, 16));
                    numInletFaces++;
                }
                continue;
            }

            if (captureOutletFaces)
            {
                std::istringstream ss(line);
                unsigned int i = faces.size() + 1;
                unsigned int num;
                std::string node1, node2, cell1, cell2;
                while (ss >> num >> node1 >> node2 >> cell1 >> cell2)
                {
                    faces.emplace_back(i, num, 5, std::stoi(node1, nullptr, 16), std::stoi(node2, nullptr, 16), std::stoi(cell1, nullptr, 16), std::stoi(cell2, nullptr, 16));
                    numOutletFaces++;
                }
                continue;
            }

            if (captureSymmetryFaces)
            {
                std::istringstream ss(line);
                unsigned int i = faces.size() + 1;
                unsigned int num;
                std::string node1, node2, cell1, cell2;
                while (ss >> num >> node1 >> node2 >> cell1 >> cell2)
                {
                    faces.emplace_back(i, num, 7, std::stoi(node1, nullptr, 16), std::stoi(node2, nullptr, 16), std::stoi(cell1, nullptr, 16), std::stoi(cell2, nullptr, 16));
                    numSymmetryFaces++;
                }
                continue;
            }
        }
    
        file.close();
    }

    void FindMaxMinNodes(const std::vector<Node> &nodes, double &xMin, double &xMax, double &yMin, double &yMax)
    {
        auto xMinNode = std::min_element(nodes.begin(), nodes.end(), [](const Node &a, const Node &b) { return a.x < b.x; });
        auto xMaxNode = std::max_element(nodes.begin(), nodes.end(), [](const Node &a, const Node &b) { return a.x < b.x; });
        auto yMinNode = std::min_element(nodes.begin(), nodes.end(), [](const Node &a, const Node &b) { return a.y < b.y; });
        auto yMaxNode = std::max_element(nodes.begin(), nodes.end(), [](const Node &a, const Node &b) { return a.y < b.y; });

        xMin = xMinNode->x;
        xMax = xMaxNode->x;
        yMin = yMinNode->y;
        yMax = yMaxNode->y;
    }

    void InputDataToCells(const unsigned int &numCells,
                          const std::vector<Face> &faces,
                          const std::vector<unsigned int> &cellsType,
                          const std::vector<Node> &nodes,
                          double &xCellAverageLength,
                          double &yCellAverageLength,
                          std::vector<Cell> &cells)
    {
        std::unordered_map<unsigned int, std::vector<Face>> FacesInCells;

        double xLengthSum = 0.0, yLengthSum = 0.0;

        cells.reserve(numCells);
        FacesInCells.reserve(numCells + 1);
    
        for (const auto &face : faces)
        {
            for (const auto &cell : face.adjacentCells)
            {
                FacesInCells[cell].emplace_back(face);
            }
        }

        for (unsigned int i = 1; i < FacesInCells.size(); ++i)
        {
            std::unordered_set<unsigned int> uniqueNodesSet, visited;
            std::unordered_map<unsigned int, std::vector<unsigned int>> adjacentNodesMap;
            std::vector<unsigned int> uniqueNodes;
            std::vector<unsigned int> uniqueFaces;

            unsigned int current = 0;
            double xLength = 0.0, yLength = 0.0;

            if (FacesInCells.at(i)[0].adjacentCells[0] == i) current = FacesInCells.at(i)[0].nodes[0];
            else current = FacesInCells.at(i)[0].nodes[1];

            for (const auto &face : FacesInCells.at(i))
            {
                uniqueNodesSet.insert(face.nodes[0]);
                uniqueNodesSet.insert(face.nodes[1]);
                adjacentNodesMap[face.nodes[0]].push_back(face.nodes[1]);
                adjacentNodesMap[face.nodes[1]].push_back(face.nodes[0]);
            }

            uniqueNodes.push_back(current);
            uniqueFaces.push_back(FacesInCells.at(i)[0].index);
            visited.insert(current);

            while (uniqueNodes.size() < uniqueNodesSet.size())
            {
                for (const unsigned int &neighbor : adjacentNodesMap[current])
                {
                    if (!visited.count(neighbor))
                    {
                        visited.emplace(neighbor);
                        uniqueNodes.push_back(neighbor);
                        current = neighbor;
                        break;
                    }
                }
            }

            unsigned int j = 1;
            while (uniqueFaces.size() < uniqueNodesSet.size())
            {
                for (const auto &face : FacesInCells.at(i))
                {
                    if (face.nodes[0] == uniqueNodes[j] || face.nodes[1] == uniqueNodes[j])
                    {
                        if (std::find(uniqueFaces.begin(), uniqueFaces.end(), face.index) == uniqueFaces.end())
                        {
                            uniqueFaces.push_back(face.index);
                            j++;
                            break;
                        }
                    }
                }
            }

            if (uniqueNodes.size() == 3)
            {
                xLength = std::max({nodes[uniqueNodes[0]-1].x, nodes[uniqueNodes[1]-1].x, nodes[uniqueNodes[2]-1].x}) - std::min({nodes[uniqueNodes[0]-1].x, nodes[uniqueNodes[1]-1].x, nodes[uniqueNodes[2]-1].x});
                yLength = std::max({nodes[uniqueNodes[0]-1].y, nodes[uniqueNodes[1]-1].y, nodes[uniqueNodes[2]-1].y}) - std::min({nodes[uniqueNodes[0]-1].y, nodes[uniqueNodes[1]-1].y, nodes[uniqueNodes[2]-1].y});
                xLengthSum += xLength;
                yLengthSum += yLength;
                uniqueNodes.push_back(uniqueNodes[2]);
                uniqueFaces.push_back(uniqueFaces[2]);
            }
            else if (uniqueNodes.size() == 4)
            {
                xLength = std::max({nodes[uniqueNodes[0]-1].x, nodes[uniqueNodes[1]-1].x, nodes[uniqueNodes[2]-1].x, nodes[uniqueNodes[3]-1].x}) - std::min({nodes[uniqueNodes[0]-1].x, nodes[uniqueNodes[1]-1].x, nodes[uniqueNodes[2]-1].x, nodes[uniqueNodes[3]-1].x});
                yLength = std::max({nodes[uniqueNodes[0]-1].y, nodes[uniqueNodes[1]-1].y, nodes[uniqueNodes[2]-1].y, nodes[uniqueNodes[3]-1].y}) - std::min({nodes[uniqueNodes[0]-1].y, nodes[uniqueNodes[1]-1].y, nodes[uniqueNodes[2]-1].y, nodes[uniqueNodes[3]-1].y});
                xLengthSum += xLength;
                yLengthSum += yLength;
            }

            // uniqueFaces.assign(FacesofCells[i].begin(), FacesofCells[i].end());

            cells.emplace_back(i, cellsType[i-1], xLength, yLength, 0.0, 0.0, 0.0, uniqueNodes, uniqueFaces);
        }

        xCellAverageLength = xLengthSum / numCells;
        yCellAverageLength = yLengthSum / numCells;
    }

    void OutputMeshData()
    {
        std::string filename = "circle.cas";
        // std::cout << "Now begin to reading mesh...\n";
        // std::cout << "Please enter the name of CAS file: ";
        // std::cin >> filename;

        ProcessFile(filename, nodes, cellsType, faces,
                    numNodes, numFaces, numBoundaryFaces,
                    numInteriorFaces, numCells, numTriCells, numQuadCells,
                    numWallFaces, numInletFaces, numOutletFaces, numSymmetryFaces);

        FindMaxMinNodes(nodes, xMin, xMax, yMin, yMax);

        cells.reserve(numCells);
        
        InputDataToCells(numCells, faces, cellsType, nodes, xCellAverageLength, yCellAverageLength, cells);
    }
}