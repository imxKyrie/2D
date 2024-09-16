#ifndef MEMORY_MANAGER_CUH
#define MEMORY_MANAGER_CUH

#include "Common.cuh"
#include <vector>

class MemoryManager
{
public:
    static std::vector<void *>& GetHostPointers();
    static std::vector<void *>& GetDevicePointers();
    static void AddHostPointer(void *ptr);
    static void AddDevicePointer(void *ptr);
    static void FreeAll();

private:
    MemoryManager() = default;
    ~MemoryManager() = default;
    MemoryManager(const MemoryManager&) = delete;
    MemoryManager& operator=(const MemoryManager&) = delete;
};

#endif // MEMORY_MANAGER_CUH