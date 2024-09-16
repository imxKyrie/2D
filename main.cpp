#include "FileData.h"
#include "GasConstants.cuh"
#include "BoundaryConditions.cuh"
#include "Mesh.cuh"
#include "MacroParams.cuh"
#include "Particles.cuh"
#include "ParticlesMotion.cuh"
#include "MemoryManager.cuh"

int main() {
    FileData::OutputMeshData();

    InitializeGas();

    InitializeBoundary();

    InitializeMesh();

    InitializeCellsMacro();

    InitializeParticles();

    for (unsigned int i = 0; i < 100; i++)
    {
        ParticlesMotion();
        std::cout << "Time step " << i << " completed" << std::endl;
    }
    verifyConstantsOnHost();
    std::cout << "End of simulation" << std::endl;

    // Free all allocated memory
    ParticlesManager::Free();
    MemoryManager::FreeAll();
}