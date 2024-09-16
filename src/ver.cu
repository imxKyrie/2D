#include "GasConstants.cuh"
#include "Common.h"
#include "BoundaryConditions.cuh"
#include "Mesh.cuh"
#include "MacroParams.cuh"
#include "Particles.cuh"

__global__ void verifyConstants(double *results,Particle* particles) {
    results[0] = d_pi;
    results[1] = d_boltz;
    results[2] = d_intError;
    results[3] = d_initCellError;
    results[4] = d_gamma;
    results[5] = particles[1].cellId;
    results[6] = d_faces[2444].index;
    results[7] = particles[1].pos.y;
    printf("id: %f\n", results[5]);
    printf("x: %f\n", results[6]);
    printf("y: %f\n", results[7]);
    // Add more as needed
}

/*
__global__ void priProperties(int *result) {

    result[0] = d_numSurfaces;
    result[1] = d_initParticles;
    printf("Number of boundary faces: ");
    printf("%f", d_numSurfaces);
    // Add more as needed
}
*/
// 在主机端验证常量内存
void verifyConstantsOnHost() {
    double* h_results = new double[8];
    double *d_results;
    CUDA_CHECK(cudaMalloc(&d_results, 8 * sizeof(double)));
    thrust::device_vector<Particle>& d_particles_temp = ParticlesManager::GetParticles();
    unsigned int num_particles = d_particles_temp.size();
    Particle* d_particles_ptr = thrust::raw_pointer_cast(d_particles_temp.data());
    // printGasProperties<<<1, 1>>>();
    verifyConstants<<<1, 1>>>(d_results, d_particles_ptr);
    CUDA_CHECK(cudaDeviceSynchronize());
    CUDA_CHECK(cudaMemcpy(h_results, d_results,8 * sizeof(double), cudaMemcpyDeviceToHost));

    std::cout << "d_pi: " << num_particles << std::endl;
    std::cout << "d_boltz: " << h_results[1] << std::endl;
    std::cout << "d_intError: " << h_results[2] << std::endl;
    std::cout << "d_initCellError: " << h_results[3] << std::endl;
    std::cout << "d_crossError: " << h_results[4] << std::endl;

    CUDA_CHECK(cudaFree(d_results));
    // CUDA_CHECK(cudaFree(d_particles_ptr));
}
/*
void verboundary() {
            int *d_result;
        int* h_result = new int[2];
    CUDA_CHECK(cudaMalloc(&d_result, 2 * sizeof(int)));

    
        // priProperties<<<1, 1>>>(d_result);
    CUDA_CHECK(cudaDeviceSynchronize());
    CUDA_CHECK(cudaMemcpy(h_result, d_result, 2*sizeof(double), cudaMemcpyDeviceToHost));
    std::cout << "Number of boundary faces: " << h_result[0] << std::endl;
}*/