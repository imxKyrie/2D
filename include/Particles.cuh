#ifndef PARTICLES_CUH
#define PARTICLES_CUH

#include "Common.cuh"
#include <thrust/device_vector.h>
#include <thrust/sort.h>

struct Particle
{
    double2 pos;
    double3 vel;
    double3 averageVel;
    double tem;
    unsigned int cellId;
    bool active;
};

class ParticlesManager
{
public:
    static thrust::device_vector<Particle>& GetParticles();
    static unsigned int GetNumParticles();
    static void Initialize(unsigned int &numParticles);
    static void Free();

private:
    static unsigned int numParticles;

    ParticlesManager() = default;
    ~ParticlesManager() = default;
    ParticlesManager(const ParticlesManager&) = delete;
    ParticlesManager& operator=(const ParticlesManager&) = delete;
};

#ifndef INIT_PARTICLES_CU

extern __device__ unsigned int d_numParticles;

#endif

void InitializeParticles();

#endif // PARTICLES_CUH