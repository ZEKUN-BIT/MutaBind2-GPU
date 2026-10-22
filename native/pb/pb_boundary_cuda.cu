// Exact-order GPU coarse boundary sum for CHARMM MAYER_BP_INT.
// Each thread owns one coarse boundary node. Sources are accumulated in the
// original atom order and PHI is rounded back to float after every atom.
// Only zero-salt, zero-membrane, no-image odd-grid calls are supported; the
// Fortran bridge retains the original interpolation and all other branches.
#include "pb_boundary_cuda.h"
#include <cuda_runtime.h>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <cstring>
#include <limits>
#include <memory>
#include <new>

namespace {
constexpr int kBoundaryThreads = 256;
constexpr size_t kBoundaryReserve = size_t(2) << 30;

struct BoundaryDeviceStorage {
    float *phi = nullptr;
    double *sources = nullptr;
    int *invalid = nullptr;
    ~BoundaryDeviceStorage() {
        if (phi) cudaFree(phi);
        if (sources) cudaFree(sources);
        if (invalid) cudaFree(invalid);
    }
};

int boundary_error(char *error, int len, const char *message) {
    if (error && len > 0) std::snprintf(error,size_t(len),"%s",message);
    return 2;
}
int boundary_cuda_error(char *error, int len, const char *stage, cudaError_t status) {
    if (error && len > 0)
        std::snprintf(error,size_t(len),"%s: %s",stage,cudaGetErrorString(status));
    return 1;
}

__global__ void boundary_coarse_sum(
    int nx, int ny, int nz, int ncharges, double dcel,
    const double *xs, const double *ys, const double *zs,
    const double *charges, const double *q_over_eps,
    size_t xy_count, size_t yz_count, size_t total_count,
    float *phi, int *invalid) {
    const size_t node = size_t(blockIdx.x)*blockDim.x + threadIdx.x;
    if (node >= total_count) return;
    int gx, gy, gz;  // zero-based grid coordinates
    // These three groups are disjoint and cover exactly the original coarse
    // nodes. XY owns all z-end nodes. YZ excludes z ends. XZ also excludes x ends.
    if (node < xy_count) {
        const size_t pair = node/2;
        gx = 2*int(pair/size_t((ny+1)/2));
        gy = 2*int(pair%size_t((ny+1)/2));
        gz = (node&1) ? nz-1 : 0;
    } else if (node < xy_count+yz_count) {
        const size_t local = node-xy_count;
        const size_t pair = local/2;
        gx = (local&1) ? nx-1 : 0;
        gy = 2*int(pair/size_t((nz-3)/2));
        gz = 2+2*int(pair%size_t((nz-3)/2));
    } else {
        const size_t local = node-xy_count-yz_count;
        const size_t pair = local/2;
        gx = 2+2*int(pair/size_t((nz-3)/2));
        gy = (local&1) ? ny-1 : 0;
        gz = 2+2*int(pair%size_t((nz-3)/2));
    }
    const size_t p = (size_t(gx)*size_t(ny)+size_t(gy))*size_t(nz)+size_t(gz);
    const double xgrid = __dmul_rn(double(gx),dcel);
    const double ygrid = __dmul_rn(double(gy),dcel);
    const double zgrid = __dmul_rn(double(gz),dcel);
    float value = phi[p];
    // Do not parallelize or reorder this sum over atoms: every addition is
    // rounded to float by the original Fortran PHI assignment.
    for (int atom=0; atom<ncharges; ++atom) {
        if (charges[atom] == 0.0) continue;
        const double dx = __dsub_rn(xgrid,xs[atom]);
        const double dy = __dsub_rn(ygrid,ys[atom]);
        const double dz = __dsub_rn(zgrid,zs[atom]);
        const double r2 = __dadd_rn(__dadd_rn(__dmul_rn(dx,dx),
                                            __dmul_rn(dy,dy)),__dmul_rn(dz,dz));
        const double distance = __dsqrt_rn(r2);
        const double increment = __ddiv_rn(q_over_eps[atom],distance);
        value = __double2float_rn(__dadd_rn(double(value),increment));
    }
    if (!isfinite(value)) atomicExch(invalid,1);
    phi[p] = value;
}
}  // namespace

extern "C" int charmm_pb_boundary_cuda(
    int nx, int ny, int nz, int ncharges, double dcel, double epsw,
    const double *x, const double *y, const double *z, const double *q,
    float *phi, char *error, int error_len) {
    if (error && error_len > 0) error[0] = '\0';
    if (nx<5 || ny<5 || nz<5 || !(nx&1) || !(ny&1) || !(nz&1) ||
        ncharges<0 || !std::isfinite(dcel) || dcel<=0 ||
        !std::isfinite(epsw) || epsw<=0)
        return boundary_error(error,error_len,"Boundary CUDA requires odd grids >=5, nonnegative charge count and finite positive spacing/dielectric");
    if (!phi || (ncharges && (!x || !y || !z || !q)))
        return boundary_error(error,error_len,"Boundary CUDA received a null input pointer");
    const size_t max_points = size_t(std::numeric_limits<int>::max());
    if (size_t(nx)>max_points/size_t(ny) ||
        size_t(nx)*size_t(ny)>max_points/size_t(nz))
        return boundary_error(error,error_len,"Boundary CUDA grid exceeds CHARMM's 32-bit indexing");
    const size_t points = size_t(nx)*size_t(ny)*size_t(nz);
    const size_t grid_bytes = points*sizeof(float);
    const size_t source_bytes = size_t(ncharges)*sizeof(double);
    const size_t xy_count = 2*size_t((nx+1)/2)*size_t((ny+1)/2);
    const size_t yz_count = 2*size_t((ny+1)/2)*size_t((nz-3)/2);
    const size_t xz_count = 2*size_t((nx-3)/2)*size_t((nz-3)/2);
    const size_t coarse_count = xy_count+yz_count+xz_count;
    if (ncharges==0) return 0;  // original atom loop is empty

    std::unique_ptr<double[]> q_over_eps(new (std::nothrow) double[size_t(ncharges)]);
    if (!q_over_eps)
        return boundary_error(error,error_len,"Cannot allocate boundary charge preprocessing buffer");
    for (int atom=0; atom<ncharges; ++atom) {
        if (!std::isfinite(x[atom]) || !std::isfinite(y[atom]) ||
            !std::isfinite(z[atom]) || !std::isfinite(q[atom]))
            return boundary_error(error,error_len,"Boundary CUDA received a nonfinite source coordinate/charge");
        // Original expression is CHI/EPSW/DIST. The first division is constant
        // per atom, and is evaluated once in double without changing its value.
        q_over_eps[atom] = q[atom]/epsw;
        if (!std::isfinite(q_over_eps[atom]))
            return boundary_error(error,error_len,"Boundary charge/dielectric division is nonfinite");
    }
    int device=0;
    cudaError_t status=cudaGetDevice(&device);
    if (status!=cudaSuccess) return boundary_cuda_error(error,error_len,"boundary cudaGetDevice",status);
    cudaDeviceProp props{};
    status=cudaGetDeviceProperties(&props,device);
    if (status!=cudaSuccess) return boundary_cuda_error(error,error_len,"boundary cudaGetDeviceProperties",status);
    size_t free_bytes=0,total_bytes=0;
    status=cudaMemGetInfo(&free_bytes,&total_bytes);
    if (status!=cudaSuccess) return boundary_cuda_error(error,error_len,"boundary cudaMemGetInfo",status);
    const size_t allocation=grid_bytes+5*source_bytes+sizeof(int);
    if (free_bytes<allocation+kBoundaryReserve+(size_t(64)<<20)) {
        if (error && error_len>0) std::snprintf(error,size_t(error_len),
            "Boundary CUDA needs %.3f GiB plus 2 GiB reserve, free %.3f GiB",
            double(allocation)/(1u<<30),double(free_bytes)/(1u<<30));
        return 3;
    }
    BoundaryDeviceStorage storage;
    const auto started=std::chrono::steady_clock::now();
    status=cudaMalloc(reinterpret_cast<void **>(&storage.phi),grid_bytes);
    if (status!=cudaSuccess) return boundary_cuda_error(error,error_len,"boundary allocate PHI",status);
    status=cudaMalloc(reinterpret_cast<void **>(&storage.sources),5*source_bytes);
    if (status!=cudaSuccess) return boundary_cuda_error(error,error_len,"boundary allocate sources",status);
    status=cudaMalloc(reinterpret_cast<void **>(&storage.invalid),sizeof(int));
    if (status!=cudaSuccess) return boundary_cuda_error(error,error_len,"boundary allocate error flag",status);
    status=cudaMemcpy(storage.phi,phi,grid_bytes,cudaMemcpyHostToDevice);
    if (status!=cudaSuccess) return boundary_cuda_error(error,error_len,"boundary copy PHI to GPU",status);
    const double *host_sources[5]={x,y,z,q,q_over_eps.get()};
    for (int i=0;i<5;++i) {
        status=cudaMemcpy(storage.sources+size_t(i)*size_t(ncharges),
                          host_sources[i],source_bytes,cudaMemcpyHostToDevice);
        if (status!=cudaSuccess) return boundary_cuda_error(error,error_len,"boundary copy sources to GPU",status);
    }
    status=cudaMemset(storage.invalid,0,sizeof(int));
    if (status!=cudaSuccess) return boundary_cuda_error(error,error_len,"boundary clear error flag",status);
    const int blocks=int((coarse_count+kBoundaryThreads-1)/kBoundaryThreads);
    boundary_coarse_sum<<<blocks,kBoundaryThreads>>>(nx,ny,nz,ncharges,dcel,
        storage.sources,storage.sources+size_t(ncharges),storage.sources+2*size_t(ncharges),
        storage.sources+3*size_t(ncharges),storage.sources+4*size_t(ncharges),
        xy_count,yz_count,coarse_count,storage.phi,storage.invalid);
    status=cudaGetLastError();
    if (status!=cudaSuccess) return boundary_cuda_error(error,error_len,"boundary launch coarse sum",status);
    int invalid=0;
    status=cudaMemcpy(&invalid,storage.invalid,sizeof(int),cudaMemcpyDeviceToHost);
    if (status!=cudaSuccess) return boundary_cuda_error(error,error_len,"boundary synchronize coarse sum",status);
    if (invalid)
        return boundary_error(error,error_len,"Boundary CUDA produced a nonfinite coarse potential (check source/grid coincidence)");
    // Never expose a partial result if a device-to-host transfer fails.
    std::unique_ptr<float[]> output(new (std::nothrow) float[points]);
    if (!output) return boundary_error(error,error_len,"Cannot allocate boundary output staging buffer");
    status=cudaMemcpy(output.get(),storage.phi,grid_bytes,cudaMemcpyDeviceToHost);
    if (status!=cudaSuccess) return boundary_cuda_error(error,error_len,"boundary copy PHI from GPU",status);
    std::memcpy(phi,output.get(),grid_bytes);
    const double seconds=std::chrono::duration<double>(
        std::chrono::steady_clock::now()-started).count();
    std::fprintf(stdout," MAYER_BP_INT CUDA: %s, sources=%d coarse_nodes=%zu allocation=%.3f GiB elapsed=%.6f s\n",
                 props.name,ncharges,coarse_count,double(allocation)/(1u<<30),seconds);
    std::fflush(stdout);
    return 0;
}
