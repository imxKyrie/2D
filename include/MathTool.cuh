#ifndef MATH_TOOL_CUH
#define MATH_TOOL_CUH

#include "Common.cuh"
#include "Mesh.cuh"
#include <curand_kernel.h>

#ifndef MATH_TOOL_CU

extern __device__ curandState *d_states;

#endif

void InitializeRandom(unsigned int &numParticles);

__device__ double TriangleArea(DeviceNode &A, DeviceNode &B, DeviceNode &C);

__device__ __host__ double2 operator+(const double2 &a, const double2 &b);
__device__ __host__ double2 operator-(const double2 &a, const double2 &b);
__device__ __host__ double3 operator+(const double3 &a, const double3 &b);
__device__ __host__ double3 operator-(const double3 &a, const double3 &b);

__device__ __host__ bool IsParticleWithinCell(double2 &P, DeviceNode (&nodes)[4], unsigned int &numFace);

__device__ __host__ double2 CalculateNormal(double2 &A, double2 &B, double2 &P);

__device__ __host__ double Erf(double x);
#endif // MATH_TOOL_CUH