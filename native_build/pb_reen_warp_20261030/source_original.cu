// CHARMM MAYER_REEN probe evaluation with original double-precision geometry.
// CPU supplies the original random points, sorted neighbor tables and ML/NL.
// Only negative-to-positive grid restoration is parallelized. The operation
// is idempotent and preserves each grid cell's original magnitude.
#include "pb_reen_cuda.h"
#include <cuda_runtime.h>
#include <cuda/atomic>
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
constexpr int kProbeThreads=128;
constexpr size_t kReserve=size_t(2)<<30;

struct ReenDeviceStorage {
    float *grids[4]{};
    int *selected=nullptr,*ml=nullptr,*nl=nullptr,*neighbors=nullptr;
    double *atoms=nullptr,*probes=nullptr;
    ~ReenDeviceStorage() {
        for (auto p:grids) if(p) cudaFree(p);
        if(selected) cudaFree(selected);
        if(ml) cudaFree(ml);
        if(nl) cudaFree(nl);
        if(neighbors) cudaFree(neighbors);
        if(atoms) cudaFree(atoms);
        if(probes) cudaFree(probes);
    }
};

int reen_error(char *error,int length,const char *message) {
    if(error && length>0) std::snprintf(error,size_t(length),"%s",message);
    return 2;
}
int reen_cuda_error(char *error,int length,const char *stage,cudaError_t status) {
    if(error && length>0)
        std::snprintf(error,size_t(length),"%s: %s",stage,cudaGetErrorString(status));
    return 1;
}

__device__ __forceinline__ double square_rn(double v) {
    return __dmul_rn(v,v);
}

__device__ __forceinline__ void restore_negative(float *cell) {
    // An atomic load avoids racing another probe's restore. A successful
    // restore is the only possible state transition, so a failed CAS means
    // another thread already performed precisely the required operation.
    cuda::atomic_ref<float,cuda::thread_scope_device> atom(*cell);
    float old=atom.load(cuda::memory_order_relaxed);
    if(old<0.0f)
        atom.compare_exchange_strong(old,-old,
            cuda::memory_order_relaxed,cuda::memory_order_relaxed);
}

__global__ void restore_accessible_probes(
    int nx,int ny,int nz,int nselected,int mapt,int neighbor_stride,int qfkap,
    double dcel,double radw,double rsmall,
    double tranx,double trany,double tranz,double xbcen,double ybcen,double zbcen,
    const int *selected,const int *ml,const int *nl,const int *neighbors,
    const double *x,const double *y,const double *z,const double *pbrad,
    const double *pox,const double *poy,const double *poz,
    float *rmay,float *rmayx,float *rmayy,float *rmayz,size_t work) {
    const size_t ncyz=size_t(ny)*size_t(nz);
    const int nfil=__double2int_rz(__ddiv_rn(radw,dcel))+2;
    const double bisq=__dadd_rn(square_rn(radw),rsmall);
    const double halfcell=__dmul_rn(0.5,dcel);
    for(size_t task=size_t(blockIdx.x)*blockDim.x+threadIdx.x;
        task<work;task+=size_t(gridDim.x)*blockDim.x) {
        const int il=int(task/size_t(mapt));
        const int jl=int(task%size_t(mapt));
        if(jl>=ml[il]) continue;
        const int i=selected[il]-1;
        const double raneti=__dadd_rn(pbrad[i],radw);
        double xp=__dadd_rn(x[i],__dmul_rn(raneti,pox[jl]));
        double yp=__dadd_rn(y[i],__dmul_rn(raneti,poy[jl]));
        double zp=__dadd_rn(z[i],__dmul_rn(raneti,poz[jl]));
        bool accessible=true;
        for(int kl=0;kl<nl[il];++kl) {
            const int ll=neighbors[size_t(il)*size_t(neighbor_stride)+size_t(kl)]-1;
            const double aij=__dadd_rn(
                __dadd_rn(square_rn(__dsub_rn(xp,x[ll])),
                          square_rn(__dsub_rn(yp,y[ll]))),
                square_rn(__dsub_rn(zp,z[ll])));
            const double radius=__dadd_rn(pbrad[ll],radw);
            if(aij<square_rn(radius)) {accessible=false;break;}
        }
        if(!accessible) continue;

        xp=__dsub_rn(__dadd_rn(xp,tranx),xbcen);
        yp=__dsub_rn(__dadd_rn(yp,trany),ybcen);
        zp=__dsub_rn(__dadd_rn(zp,tranz),zbcen);
        const int ix=__double2int_rz(__ddiv_rn(xp,dcel))+1;
        const int iy=__double2int_rz(__ddiv_rn(yp,dcel))+1;
        const int iz=__double2int_rz(__ddiv_rn(zp,dcel))+1;
        const int jx1=max(1,ix-nfil+1),jx2=min(nx,ix+nfil);
        const int jy1=max(1,iy-nfil+1),jy2=min(ny,iy+nfil);
        const int jz1=max(1,iz-nfil+1),jz2=min(nz,iz+nfil);
        for(int k=jx1;k<=jx2;++k) {
            const size_t ipx=size_t(k-1)*ncyz;
            const double xc=__dmul_rn(double(k-1),dcel);
            const double xsq=square_rn(__dsub_rn(xc,xp));
            for(int l=jy1;l<=jy2;++l) {
                const size_t ipy=size_t(l-1)*size_t(nz);
                const double yc=__dmul_rn(double(l-1),dcel);
                const double ysq=square_rn(__dsub_rn(yc,yp));
                for(int m=jz1;m<=jz2;++m) {
                    const size_t p=ipx+ipy+size_t(m-1);
                    const double zc=__dmul_rn(double(m-1),dcel);
                    const double zsq=square_rn(__dsub_rn(zc,zp));
                    const double dsq=__dadd_rn(__dadd_rn(xsq,ysq),zsq);
                    if(!qfkap && dsq<=bisq) restore_negative(rmay+p);
                    // Calculating the predicates before loading the signs
                    // avoids atomic reads outside the probe sphere. Geometry
                    // is independent of mutable grid values in the reference.
                    const double xc1=__dadd_rn(xc,halfcell);
                    const double dx=__dadd_rn(
                        __dadd_rn(square_rn(__dsub_rn(xc1,xp)),ysq),zsq);
                    if(dx<=bisq) restore_negative(rmayx+p);
                    const double yc1=__dadd_rn(yc,halfcell);
                    const double dy=__dadd_rn(
                        __dadd_rn(square_rn(__dsub_rn(yc1,yp)),xsq),zsq);
                    if(dy<=bisq) restore_negative(rmayy+p);
                    const double zc1=__dadd_rn(zc,halfcell);
                    // Preserve the original Z expression's Y-before-X order.
                    const double dz=__dadd_rn(
                        __dadd_rn(square_rn(__dsub_rn(zc1,zp)),ysq),xsq);
                    if(dz<=bisq) restore_negative(rmayz+p);
                }
            }
        }
    }
}
}  // namespace

extern "C" int charmm_pb_reen_cuda(
    int nx,int ny,int nz,int natom,int nselected,int mapt,int neighbor_stride,int qfkap,
    double dcel,double radw,double epsp,double rsmall,
    double tranx,double trany,double tranz,double xbcen,double ybcen,double zbcen,
    const int *selected,const int *ml,const int *nl,const int *neighbors,
    const double *x,const double *y,const double *z,const double *pbrad,
    const double *pox,const double *poy,const double *poz,
    float *rmay,float *rmayx,float *rmayy,float *rmayz,char *error,int error_len) {
    if(error && error_len>0) error[0]='\0';
    if(nx<1 || ny<1 || nz<1 || natom<1 || nselected<0 || mapt<1 || neighbor_stride<1 || neighbor_stride>512 ||
       (qfkap!=0 && qfkap!=1) || !std::isfinite(dcel) || dcel<=0 ||
       !std::isfinite(radw) || radw<0 || !std::isfinite(epsp) ||
       !std::isfinite(rsmall) || rsmall<0 || !std::isfinite(tranx) ||
       !std::isfinite(trany) || !std::isfinite(tranz) ||
       !std::isfinite(xbcen) || !std::isfinite(ybcen) || !std::isfinite(zbcen))
        return reen_error(error,error_len,"REEN CUDA received invalid dimensions or geometry scalars");
    if(!rmay || !rmayx || !rmayy || !rmayz || !selected || !ml || !nl ||
       !neighbors || !x || !y || !z || !pbrad || !pox || !poy || !poz)
        return reen_error(error,error_len,"REEN CUDA received a null array pointer");
    const size_t max_points=size_t(std::numeric_limits<int>::max());
    if(size_t(nx)>max_points/size_t(ny) || size_t(nx)*size_t(ny)>max_points/size_t(nz))
        return reen_error(error,error_len,"REEN CUDA grid exceeds CHARMM's 32-bit indexing");
    if(radw/dcel>1e6)
        return reen_error(error,error_len,"REEN CUDA probe-to-grid ratio exceeds safe index bounds");
    if(nselected==0) return 0;
    bool any_probes=false;
    for(int i=0;i<natom;++i) {
        if(!std::isfinite(x[i]) || !std::isfinite(y[i]) || !std::isfinite(z[i]) || !std::isfinite(pbrad[i]))
            return reen_error(error,error_len,"REEN CUDA received nonfinite atom coordinates/radius");
    }
    double max_point=0.0;
    for(int p=0;p<mapt;++p) {
        if(!std::isfinite(pox[p]) || !std::isfinite(poy[p]) || !std::isfinite(poz[p]))
            return reen_error(error,error_len,"REEN CUDA received nonfinite probe coordinates");
        max_point=std::max(max_point,std::max(std::abs(pox[p]),std::max(std::abs(poy[p]),std::abs(poz[p]))));
    }
    const double max_shift=std::max(std::abs(tranx)+std::abs(xbcen),
        std::max(std::abs(trany)+std::abs(ybcen),std::abs(tranz)+std::abs(zbcen)));
    for(int atom=0;atom<nselected;++atom) {
        if(selected[atom]<1 || selected[atom]>natom || ml[atom]<0 || ml[atom]>mapt || nl[atom]<-1 ||
           nl[atom]>neighbor_stride || (nl[atom]<0 && ml[atom]!=0))
            return reen_error(error,error_len,"REEN CUDA received an invalid atom identifier or ML/NL table");
        const int i=selected[atom]-1;
        const double bound=std::max(std::abs(x[i]),std::max(std::abs(y[i]),std::abs(z[i])))+
                           std::abs(pbrad[i]+radw)*max_point+max_shift;
        if(!std::isfinite(bound) || bound/dcel>1e9)
            return reen_error(error,error_len,"REEN CUDA transformed probe can exceed safe grid index bounds");
        for(int j=0;j<nl[atom];++j) {
            const int neighbor=neighbors[size_t(atom)*size_t(neighbor_stride)+size_t(j)];
            if(neighbor<1 || neighbor>natom)
                return reen_error(error,error_len,"REEN CUDA received an invalid neighbor atom identifier");
        }
        any_probes=any_probes || ml[atom]>0;
    }
    if(!any_probes) return 0;
    const size_t points=size_t(nx)*size_t(ny)*size_t(nz);
    const size_t grid_bytes=points*sizeof(float);
    const size_t selected_bytes=size_t(nselected)*sizeof(int);
    const size_t neighbor_bytes=size_t(nselected)*size_t(neighbor_stride)*sizeof(int);
    const size_t atom_bytes=size_t(natom)*sizeof(double);
    const size_t probe_bytes=size_t(mapt)*sizeof(double);
    const size_t allocation=4*grid_bytes+3*selected_bytes+neighbor_bytes+4*atom_bytes+3*probe_bytes;
    int device=0;
    cudaError_t status=cudaGetDevice(&device);
    if(status!=cudaSuccess) return reen_cuda_error(error,error_len,"REEN cudaGetDevice",status);
    cudaDeviceProp props{};
    status=cudaGetDeviceProperties(&props,device);
    if(status!=cudaSuccess) return reen_cuda_error(error,error_len,"REEN cudaGetDeviceProperties",status);
    size_t free_bytes=0,total_bytes=0;
    status=cudaMemGetInfo(&free_bytes,&total_bytes);
    if(status!=cudaSuccess) return reen_cuda_error(error,error_len,"REEN cudaMemGetInfo",status);
    if(free_bytes<allocation+kReserve+(size_t(64)<<20)) {
        if(error && error_len>0) std::snprintf(error,size_t(error_len),
            "REEN CUDA needs %.3f GiB plus 2 GiB reserve, free %.3f GiB",
            double(allocation)/(1u<<30),double(free_bytes)/(1u<<30));
        return 3;
    }
    const auto started=std::chrono::steady_clock::now();
    ReenDeviceStorage storage;
    float *host_grids[4]={rmay,rmayx,rmayy,rmayz};
    for(int g=0;g<4;++g) {
        status=cudaMalloc(reinterpret_cast<void **>(&storage.grids[g]),grid_bytes);
        if(status!=cudaSuccess) return reen_cuda_error(error,error_len,"REEN allocate grid",status);
        status=cudaMemcpy(storage.grids[g],host_grids[g],grid_bytes,cudaMemcpyHostToDevice);
        if(status!=cudaSuccess) return reen_cuda_error(error,error_len,"REEN copy grid to GPU",status);
    }
    int **device_integer_arrays[4]={&storage.selected,&storage.ml,&storage.nl,&storage.neighbors};
    const int *host_integer_arrays[4]={selected,ml,nl,neighbors};
    for(int t=0;t<4;++t) {
        const size_t bytes=t==3?neighbor_bytes:selected_bytes;
        status=cudaMalloc(reinterpret_cast<void **>(device_integer_arrays[t]),bytes);
        if(status!=cudaSuccess) return reen_cuda_error(error,error_len,"REEN allocate atom tables",status);
        status=cudaMemcpy(*device_integer_arrays[t],host_integer_arrays[t],bytes,cudaMemcpyHostToDevice);
        if(status!=cudaSuccess) return reen_cuda_error(error,error_len,"REEN copy atom tables",status);
    }
    status=cudaMalloc(reinterpret_cast<void **>(&storage.atoms),4*atom_bytes);
    if(status!=cudaSuccess) return reen_cuda_error(error,error_len,"REEN allocate atom coordinates",status);
    const double *host_atoms[4]={x,y,z,pbrad};
    for(int a=0;a<4;++a) {
        status=cudaMemcpy(storage.atoms+size_t(a)*size_t(natom),host_atoms[a],atom_bytes,cudaMemcpyHostToDevice);
        if(status!=cudaSuccess) return reen_cuda_error(error,error_len,"REEN copy atom coordinates",status);
    }
    status=cudaMalloc(reinterpret_cast<void **>(&storage.probes),3*probe_bytes);
    if(status!=cudaSuccess) return reen_cuda_error(error,error_len,"REEN allocate probe coordinates",status);
    const double *host_probes[3]={pox,poy,poz};
    for(int p=0;p<3;++p) {
        status=cudaMemcpy(storage.probes+size_t(p)*size_t(mapt),host_probes[p],probe_bytes,cudaMemcpyHostToDevice);
        if(status!=cudaSuccess) return reen_cuda_error(error,error_len,"REEN copy probe coordinates",status);
    }
    const size_t work=size_t(nselected)*size_t(mapt);
    const int blocks=int(std::min<size_t>((work+kProbeThreads-1)/kProbeThreads,
        size_t(std::max(1,props.multiProcessorCount))*32));
    restore_accessible_probes<<<blocks,kProbeThreads>>>(nx,ny,nz,nselected,mapt,neighbor_stride,qfkap,
        dcel,radw,rsmall,tranx,trany,tranz,xbcen,ybcen,zbcen,
        storage.selected,storage.ml,storage.nl,storage.neighbors,
        storage.atoms,storage.atoms+size_t(natom),storage.atoms+2*size_t(natom),storage.atoms+3*size_t(natom),
        storage.probes,storage.probes+size_t(mapt),storage.probes+2*size_t(mapt),
        storage.grids[0],storage.grids[1],storage.grids[2],storage.grids[3],work);
    status=cudaGetLastError();
    if(status!=cudaSuccess) return reen_cuda_error(error,error_len,"REEN launch probe kernel",status);
    status=cudaDeviceSynchronize();
    if(status!=cudaSuccess) return reen_cuda_error(error,error_len,"REEN synchronize probe kernel",status);
    // Acquire all staging storage before transferring or modifying any output.
    std::unique_ptr<float[]> staging[4];
    for(auto &buffer:staging) {
        buffer.reset(new(std::nothrow) float[points]);
        if(!buffer) return reen_error(error,error_len,"Cannot allocate REEN transactional host output buffers");
    }
    for(int g=0;g<4;++g) {
        status=cudaMemcpy(staging[g].get(),storage.grids[g],grid_bytes,cudaMemcpyDeviceToHost);
        if(status!=cudaSuccess) return reen_cuda_error(error,error_len,"REEN copy output grid from GPU",status);
    }
    for(int g=0;g<4;++g) std::memcpy(host_grids[g],staging[g].get(),grid_bytes);
    const double seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-started).count();
    std::fprintf(stdout," MAYER_REEN CUDA: %s, selected=%d points=%d allocation=%.3f GiB elapsed=%.6f s\n",
        props.name,nselected,mapt,double(allocation)/(1u<<30),seconds);
    std::fflush(stdout);
    return 0;
}
