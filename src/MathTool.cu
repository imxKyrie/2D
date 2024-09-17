#define MATH_TOOL_CU

#include "MathTool.cuh"
#include "MemoryManager.cuh"

__device__ curandState *d_states;

// Random number generation
__global__ void SetKernel(curandState *state, unsigned long seed) {
    unsigned int idx = blockIdx.x * blockDim.x + threadIdx.x;
    curand_init(seed, idx, 0, &state[idx]);
}

void InitializeRandom(unsigned int& numParticles) {
    unsigned int blockSize = 128;
    unsigned int gridSize = (3 * numParticles + blockSize - 1) / blockSize;

    curandState *d_states_temp;
    CUDA_CHECK(cudaMalloc((void**)&d_states_temp, blockSize * gridSize * sizeof(curandState)));
    
    dim3 block(blockSize);
    dim3 grid(gridSize);

    SetKernel<<<grid, block>>>(d_states_temp, time(nullptr));
    CUDA_CHECK(cudaDeviceSynchronize());

    CUDA_CHECK(cudaMemcpyToSymbol(d_states, &d_states_temp, sizeof(curandState*)));

    MemoryManager::AddDevicePointer(d_states_temp);
}

// Calculate area of a triangle
__device__ double TriangleArea(DeviceNode &A, DeviceNode &B, DeviceNode &C) {
    return 0.5 * fabs(A.x * (B.y - C.y) + B.x * (C.y - A.y) + C.x * (A.y - B.y));
}

// Calculate vector
__device__ __host__ double2 operator+(const double2 &A, const double2 &B) {
    double2 result;
    result.x = A.x + B.x;
    result.y = A.y + B.y;
    return result;
}

__device__ __host__ double2 operator-(const double2 &A, const double2 &B) {
    double2 result;
    result.x = A.x - B.x;
    result.y = A.y - B.y;
    return result;
}

__device__ __host__ double3 operator+(const double3 &A, const double3 &B) {
    double3 result;
    result.x = A.x + B.x;
    result.y = A.y + B.y;
    result.z = A.z + B.z;
    return result;
}

__device__ __host__ double3 operator-(const double3 &A, const double3 &B) {
    double3 result;
    result.x = A.x - B.x;
    result.y = A.y - B.y;
    result.z = A.z - B.z;
    return result;
}

// Determine whether it is within cell
__device__ __host__ double CrossProduct(const double2 &A, const double2 &B, const double2 &P) {
    return (B.x - A.x) * (P.y - A.y) - (B.y - A.y) * (P.x - A.x);
}

__device__ __host__ bool IsParticleWithinCell(double2 &P, DeviceNode (&nodes)[4], unsigned int &numFace) {
    double cross[4];
    for (unsigned int i = 0; i < numFace; i++) {
        cross[i] = CrossProduct(nodes[i], nodes[(i + 1) % numFace], P);
    }

    if (numFace == 3) return (cross[0] > 0 && cross[1] > 0 && cross[2] > 0) || (cross[0] < 0 && cross[1] < 0 && cross[2] < 0);
    else return (cross[0] > 0 && cross[1] > 0 && cross[2] > 0 && cross[3] > 0) || (cross[0] < 0 && cross[1] < 0 && cross[2] < 0 && cross[3] < 0);
}

// Calculate normal of face
__device__ __host__ double2 CalculateNormal(double2 &A, double2 &B, double2 &P) {
    double2 AB = B - A;
    double2 AP = P - A;
    double2 normal = make_double2(-AB.y, AB.x);
    if (normal.x * AP.x + normal.y * AP.y < 0) {
        normal.x = -normal.x;
        normal.y = -normal.y;
    }
    normal.x = normal.x / sqrt(normal.x * normal.x + normal.y * normal.y);
    normal.y = normal.y / sqrt(normal.x * normal.x + normal.y * normal.y);

    return normal;
}

// Error Function and Gamma Function
__device__ __host__ double Gamma_1(double x) {
    if (x <= 0.0) {
        printf("Gamma_1 Error");
        return -1.0;
    }
    double data[11] = {0.0000677106, -0.0003442342, 0.0015397681, -0.0024467480, 0.0109736958,
                       -0.0002109075,0.0742379071,0.0815782188,0.4118402518,0.4227843370,1.0};
    double gamma = data[0];
    double t = 0.0;
    if (x <= 1.0) {
        t = 1.0 / (x * (x + 1.0));
        x += 2.0;
    }
    else if (x <= 2.0) {
        t = 1.0 / x;
        x += 1.0;
    }
    else if (x <= 3.0) {
        t = 1.0;
    }
    else {
        t = 1.0;
        while (x > 3.0) {
            x -= 1.0;
            t *= x;
        }
    }
    double y = x - 2.0;
    for (unsigned int i = 0; i < 10; i++) {
        gamma = gamma * y + data[i + 1];
    }
    gamma *= t;

    return gamma;
}

__device__ __host__ double Gamma_2(double a, double x) {
    if (a <= 0.0 || x <= 0.0) {
        printf("Gamma_2 Error");
        return -1.0;
    }
    if (fabs(x) < 1.0e-10) return 0.0;
    if (x > 1.0e35) return 1.0;
    double gamma = 0.0;
    double q = exp(log(x) * a);
    if (x < 1.0 + a) {
        double d = 1.0 / a, p = a;
        gamma = d;
        for (unsigned int i = 1; i <= 100; i++) {
            p = 1.0 + p;
            d = d * x / p;
            gamma += d;
            if (fabs(d) < 1.0e-7 * fabs(gamma)) {
                gamma = gamma * exp(-x) * q / Gamma_1(a);
                return gamma;
            }
        }
    }
    else {
        double p0 = 0.0, p1 = 1.0, q0 = 1.0, q1 = x;
        gamma = 1.0 / x;
        for (unsigned int i = 1; i <= 100; i++) {
            p0 = p1 + (i - a) * p0;
            q0 = q1 + (i - a) * q0;
            p1 = x * p0 + i * p1;
            q1 = x * q0 + i * q1;
            if (fabs(q1) > 1.0e-10) {
                double gamma0 = p1 / q1;
                if (fabs((gamma0 - gamma)/gamma0) < 1.0e-7) {
                    gamma = gamma0 * exp(-x) * q / Gamma_1(a);
                    return 1.0 - gamma;
                }
                gamma = gamma0;
            }
        }
    }
    printf("A too large!\n");
    gamma = 1.0 - gamma * exp(-x) * q / Gamma_1(a);

    return gamma;
}

__device__ __host__ double Erf(double x) {
    return (x >= 0.0) ? Gamma_2(0.5, x * x) : -Gamma_2(0.5, x * x);
}