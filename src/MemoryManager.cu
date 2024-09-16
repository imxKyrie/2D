#include "MemoryManager.cuh"

std::vector<void *>& MemoryManager::GetHostPointers() {
    static std::vector<void *> h_pointers;
    return h_pointers;
}

std::vector<void *>& MemoryManager::GetDevicePointers() {
    static std::vector<void *> d_pointers;
    return d_pointers;
}

void MemoryManager::AddHostPointer(void *ptr) {
    GetHostPointers().push_back(ptr);
}

void MemoryManager::AddDevicePointer(void *ptr) {
    GetDevicePointers().push_back(ptr);
}

void MemoryManager::FreeAll() {
    for (void *ptr : GetHostPointers()) CUDA_CHECK(cudaFreeHost(ptr));
    for (void *ptr : GetDevicePointers()) CUDA_CHECK(cudaFree(ptr));
    GetHostPointers().clear();
    GetDevicePointers().clear();
}