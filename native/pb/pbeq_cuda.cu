// Restricted GPU implementation of CHARMM PBEQ1's nonperiodic odd-grid branch.
// Original reference: CHARMM 51 source/misc/pbeq.F90, SUBROUTINE PBEQ1.
// Compile with --fmad=false and without --use_fast_math. The explicitly rounded
// operations preserve the original mixed float/double expression evaluation.
#include "pbeq_cuda.h"
#include <cuda_runtime.h>
#include <cub/block/block_reduce.cuh>
#include <cub/device/device_reduce.cuh>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <cstring>
#include <limits>
#include <memory>
#include <new>

namespace {
constexpr int kThreads = 256;
constexpr size_t kReserveBytes = size_t(2) << 30;

struct DeviceStorage {
    float *grids[6]{};
    double *partials = nullptr;
    double *norm = nullptr;
    void *scratch = nullptr;
    ~DeviceStorage() {
        for (auto p : grids) if (p) cudaFree(p);
        if (partials) cudaFree(partials);
        if (norm) cudaFree(norm);
        if (scratch) cudaFree(scratch);
    }
};

int fail(char *error, int len, const char *stage, cudaError_t status) {
    if (error && len > 0)
        std::snprintf(error, size_t(len), "%s: %s", stage, cudaGetErrorString(status));
    return 1;
}
int invalid(char *error, int len, const char *message) {
    if (error && len > 0) std::snprintf(error, size_t(len), "%s", message);
    return 2;
}

__device__ __forceinline__ float sum6(float a, float b, float c,
                                     float d, float e, float f) {
    // Fortran evaluates this all-float subexpression before promoting it.
    float v = __fadd_rn(a, b);
    v = __fadd_rn(v, c);
    v = __fadd_rn(v, d);
    v = __fadd_rn(v, e);
    return __fadd_rn(v, f);
}

__global__ void update_color(
    int nx, int ny, int nz, int xmin, int xmax, int ymin, int ymax,
    int zbase, int zmaxbase, int iteration, double omega, double epsw,
    float *phi, const float *epsx, const float *epsy, const float *epsz,
    const float *cden, const float *fkapa, size_t count,
    double *block_norms) {
    using Reduce = cub::BlockReduce<double, kThreads>;
    __shared__ typename Reduce::TempStorage storage;
    const size_t yz = size_t(ny) * size_t(nz);
    const int zslots = (nz - 1) / 2;
    const int yrows = ny - 2;
    double thread_norm = 0.0;
    for (size_t q = size_t(blockIdx.x) * blockDim.x + threadIdx.x;
         q < count; q += size_t(blockDim.x) * gridDim.x) {
        const size_t row = q / size_t(zslots);
        // x/y/z are zero-based; i/j/kg below use the original Fortran indices.
        const int x = int(row / size_t(yrows)) + 1;
        const int y = int(row % size_t(yrows)) + 1;
        const int ksw = 1 + ((iteration + x + y - 3) & 1);
        const int z = ksw + 2 * int(q % size_t(zslots));
        if (z >= nz - 1) continue;
        const int i = x + 1, j = y + 1, kg = z + 1;
        const size_t p = size_t(x) * yz + size_t(y) * size_t(nz) + size_t(z);
        const bool left = i < xmin;
        const int zlo = zbase + ksw - 1;
        const int zhi = zmaxbase - ksw + 1;
        const bool heterogeneous = !left && i <= xmax && j >= ymin &&
                                   j <= ymax && kg >= zlo && kg <= zhi;
        double e123, residual;
        if (heterogeneous) {
            float ef = sum6(epsx[p], epsy[p], epsz[p],
                            epsx[p-yz], epsy[p-nz], epsz[p-1]);
            e123 = double(__fadd_rn(ef, fkapa[p]));
            float weighted = sum6(
                __fmul_rn(phi[p-yz], epsx[p-yz]),
                __fmul_rn(phi[p+yz], epsx[p]),
                __fmul_rn(phi[p-nz], epsy[p-nz]),
                __fmul_rn(phi[p+nz], epsy[p]),
                __fmul_rn(phi[p-1], epsz[p-1]),
                __fmul_rn(phi[p+1], epsz[p]));
            residual = __dsub_rn(double(weighted), __dmul_rn(double(phi[p]), e123));
            residual = __dadd_rn(residual, double(cden[p]));
        } else {
            e123 = __dadd_rn(__dmul_rn(6.0, epsw), double(fkapa[p]));
            float neighbors = sum6(phi[p-yz], phi[p+yz], phi[p-nz],
                                   phi[p+nz], phi[p-1], phi[p+1]);
            residual = __dsub_rn(__dmul_rn(double(neighbors), epsw),
                                __dmul_rn(double(phi[p]), e123));
            // The original left-X region includes CDEN; all other uniform
            // regions omit it. This asymmetry must not be "corrected" here.
            if (left) residual = __dadd_rn(residual, double(cden[p]));
        }
        thread_norm = __dadd_rn(thread_norm, fabs(residual));
        phi[p] = __double2float_rn(__dadd_rn(double(phi[p]),
            __ddiv_rn(__dmul_rn(omega, residual), e123)));
    }
    const double sum = Reduce(storage).Sum(thread_norm);
    if (threadIdx.x == 0) block_norms[blockIdx.x] = sum;
}
}  // namespace

static int pbeq_cuda_core(
    int nx, int ny, int nz,
    int ixmin, int ixmax, int iymin, int iymax, int izmin, int izmax,
    int maxits, double omega, double tol, double epsw, double rjac,
    float *phi, const float *epsx, const float *epsy, const float *epsz,
    const float *cden, const float *fkapa,
    int *iterations, double *anorm, char *error, int error_len,
    double *norm_trace, double *omega_trace, float *center_trace,
    int trace_capacity) {
    if (error && error_len > 0) error[0] = '\0';
    if (!phi || !epsx || !epsy || !epsz || !cden || !fkapa ||
        !iterations || !anorm)
        return invalid(error, error_len, "PBEQ CUDA received a null array/output pointer");
    if (nx < 5 || ny < 5 || nz < 5 || !(nx&1) || !(ny&1) || !(nz&1) ||
        maxits < 1 || !std::isfinite(omega) || !std::isfinite(tol) ||
        !std::isfinite(epsw) || !std::isfinite(rjac) || tol < 0 || epsw <= 0)
        return invalid(error, error_len, "PBEQ CUDA requires odd grids >=5 and finite valid scalars");
    const bool tracing = norm_trace || omega_trace || center_trace || trace_capacity;
    if (tracing && (!norm_trace || !omega_trace || !center_trace || trace_capacity < maxits))
        return invalid(error,error_len,"PBEQ CUDA trace requires three output arrays with capacity >= MAXITS");
    const bool vacuum = ixmin == nx && ixmax == nx-1;
    if (!vacuum && !(ixmin >= 2 && ixmin <= ixmax && ixmax <= nx-1 &&
        iymin >= 2 && iymin <= iymax && iymax <= ny-1 &&
        ((ixmax-ixmin)&1) == 0 && ((iymax-iymin)&1) == 0 &&
        ((ixmin+iymin)&1) == 0 && izmin >= 2 && izmin < izmax &&
        izmax <= nz-1 && (izmin&1) == 0 && (izmax&1) == 1))
        return invalid(error, error_len, "PBEQ CUDA received unsupported adjusted bounds");
    const size_t max_points = size_t(std::numeric_limits<int>::max());
    if (size_t(nx) > max_points/size_t(ny) ||
        size_t(nx)*size_t(ny) > max_points/size_t(nz))
        return invalid(error, error_len, "PBEQ CUDA grid exceeds CHARMM's 32-bit indexing");
    const size_t points = size_t(nx) * size_t(ny) * size_t(nz);
    // Original PBEQ1 NC3 and grid indexing are 32-bit integers.
    if (points > size_t(std::numeric_limits<int>::max()))
        return invalid(error, error_len, "PBEQ CUDA grid exceeds CHARMM's 32-bit indexing");
    const size_t bytes = points * sizeof(float);
    const size_t center = size_t(ny)*size_t(nz)*size_t(nx/2) +
                          size_t(nz)*size_t(ny/2) + size_t(nz/2);
    const size_t work = size_t(nx-2) * size_t(ny-2) * size_t((nz-1)/2);
    int device = 0;
    cudaError_t status = cudaGetDevice(&device);
    if (status != cudaSuccess) return fail(error, error_len, "cudaGetDevice", status);
    cudaDeviceProp props{};
    status = cudaGetDeviceProperties(&props, device);
    if (status != cudaSuccess) return fail(error, error_len, "cudaGetDeviceProperties", status);
    const int blocks = int(std::min<size_t>((work+kThreads-1)/kThreads,
                           size_t(std::max(1, props.multiProcessorCount)) * 32));
    size_t scratch_bytes = 0;
    status = cub::DeviceReduce::Sum(nullptr, scratch_bytes,
        static_cast<double *>(nullptr), static_cast<double *>(nullptr), blocks);
    if (status != cudaSuccess) return fail(error, error_len, "CUB workspace query", status);
    size_t free_bytes=0, total_bytes=0;
    status = cudaMemGetInfo(&free_bytes, &total_bytes);
    if (status != cudaSuccess) return fail(error, error_len, "cudaMemGetInfo", status);
    const size_t allocation = 6*bytes + size_t(blocks+1)*sizeof(double) + scratch_bytes;
    // Also reserve a small allowance for runtime/module allocations after this
    // query. Keep at least 2 GiB free for the display and other user workloads.
    if (free_bytes < allocation + kReserveBytes + (size_t(64)<<20)) {
        if (error && error_len > 0) std::snprintf(error, size_t(error_len),
            "Insufficient GPU memory: need %.3f GiB plus 2 GiB reserve, free %.3f GiB",
            double(allocation)/(1u<<30), double(free_bytes)/(1u<<30));
        return 3;
    }
    DeviceStorage storage;
    const auto started = std::chrono::steady_clock::now();
    const float *host[6] = {phi, epsx, epsy, epsz, cden, fkapa};
    for (int i=0; i<6; ++i) {
        status = cudaMalloc(reinterpret_cast<void **>(&storage.grids[i]), bytes);
        if (status != cudaSuccess) return fail(error, error_len, "cudaMalloc grid", status);
        status = cudaMemcpy(storage.grids[i], host[i], bytes, cudaMemcpyHostToDevice);
        if (status != cudaSuccess) return fail(error, error_len, "copy grid to GPU", status);
    }
    status = cudaMalloc(reinterpret_cast<void **>(&storage.partials), size_t(blocks)*sizeof(double));
    if (status != cudaSuccess) return fail(error, error_len, "cudaMalloc partial norms", status);
    status = cudaMalloc(reinterpret_cast<void **>(&storage.norm), sizeof(double));
    if (status != cudaSuccess) return fail(error, error_len, "cudaMalloc norm", status);
    status = cudaMalloc(&storage.scratch, scratch_bytes);
    if (status != cudaSuccess) return fail(error, error_len, "cudaMalloc reduction workspace", status);
    std::fprintf(stdout, " PBEQ1 CUDA: %s, grid %dx%dx%d, allocation %.3f GiB, reserved 2 GiB\n",
        props.name, nx, ny, nz, double(allocation)/(1u<<30));
    std::fflush(stdout);
    int n = 1;
    double last_norm = 0.0;
    const double rjac2 = rjac*rjac;
    for (; n <= maxits; ++n) {
        update_color<<<blocks,kThreads>>>(nx,ny,nz,ixmin,ixmax,iymin,iymax,
            izmin,izmax,n,omega,epsw,storage.grids[0],storage.grids[1],
            storage.grids[2],storage.grids[3],storage.grids[4],storage.grids[5],
            work,storage.partials);
        status = cudaGetLastError();
        if (status != cudaSuccess) return fail(error, error_len, "launch PBEQ color", status);
        status = cub::DeviceReduce::Sum(storage.scratch, scratch_bytes,
                                      storage.partials,storage.norm,blocks);
        if (status != cudaSuccess) return fail(error, error_len, "reduce PBEQ norm", status);
        // A synchronous scalar transfer completes this color before the next
        // launches. All six large grids remain resident throughout the solve.
        double norm_sum = 0.0;
        status = cudaMemcpy(&norm_sum,storage.norm,sizeof(double),cudaMemcpyDeviceToHost);
        if (status != cudaSuccess) return fail(error, error_len, "copy PBEQ norm from GPU", status);
        last_norm = norm_sum/double(points);
        if (!std::isfinite(last_norm))
            return invalid(error,error_len,"PBEQ CUDA produced a nonfinite residual norm");
        if (tracing) {
            float center_value=0.0f;
            status=cudaMemcpy(&center_value,storage.grids[0]+center,sizeof(float),cudaMemcpyDeviceToHost);
            if (status!=cudaSuccess) return fail(error,error_len,"copy PBEQ trace center from GPU",status);
            norm_trace[n-1]=last_norm;
            omega_trace[n-1]=omega;
            center_trace[n-1]=center_value;
        }
        if (n == 1) omega = 1.0/(1.0-0.5*rjac2);
        else omega = 1.0/(1.0-0.25*rjac2*omega);
        if (n > 1 && last_norm < tol) break;
    }
    // Commit host PHI only after the entire device-to-host transfer succeeds.
    // This also keeps the caller's input intact after a transfer/driver error.
    std::unique_ptr<float[]> output(new (std::nothrow) float[points]);
    if (!output) return invalid(error,error_len,"Cannot allocate host potential output staging buffer");
    status = cudaMemcpy(output.get(),storage.grids[0],bytes,cudaMemcpyDeviceToHost);
    if (status != cudaSuccess) return fail(error,error_len,"copy potential from GPU",status);
    std::memcpy(phi,output.get(),bytes);
    *iterations = n;
    *anorm = last_norm;
    const double seconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now()-started).count();
    std::fprintf(stdout," PBEQ1 CUDA: iterations=%d ANORM=%.17g elapsed=%.6f s\n",
                 n,last_norm,seconds);
    std::fflush(stdout);
    return 0;
}

extern "C" int charmm_pbeq1_cuda(
    int nx,int ny,int nz,int ixmin,int ixmax,int iymin,int iymax,int izmin,int izmax,
    int maxits,double omega,double tol,double epsw,double rjac,
    float *phi,const float *epsx,const float *epsy,const float *epsz,
    const float *cden,const float *fkapa,int *iterations,double *anorm,
    char *error,int error_len) {
    return pbeq_cuda_core(nx,ny,nz,ixmin,ixmax,iymin,iymax,izmin,izmax,
        maxits,omega,tol,epsw,rjac,phi,epsx,epsy,epsz,cden,fkapa,
        iterations,anorm,error,error_len,nullptr,nullptr,nullptr,0);
}

extern "C" int charmm_pbeq1_cuda_trace(
    int nx,int ny,int nz,int ixmin,int ixmax,int iymin,int iymax,int izmin,int izmax,
    int maxits,double omega,double tol,double epsw,double rjac,
    float *phi,const float *epsx,const float *epsy,const float *epsz,
    const float *cden,const float *fkapa,int *iterations,double *anorm,
    char *error,int error_len,double *norm_trace,double *omega_trace,
    float *center_trace,int trace_capacity) {
    return pbeq_cuda_core(nx,ny,nz,ixmin,ixmax,iymin,iymax,izmin,izmax,
        maxits,omega,tol,epsw,rjac,phi,epsx,epsy,epsz,cden,fkapa,
        iterations,anorm,error,error_len,norm_trace,omega_trace,center_trace,trace_capacity);
}
