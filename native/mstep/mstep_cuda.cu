// MAYER_MSTEP: preserve double expression order and float grid values.
// Each atom contributes only shell/terminal flags. OR is exact and commutative;
// terminal dominates shell. Apply once to the original grids, after all atoms.
#include <cuda_runtime.h>
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <climits>
#include <memory>
#include <new>
#include <vector>

struct Atom { double x,y,z,r,rw,ri; int x0,x1,y0,y1,z0,z1; };
namespace {
struct Device {
    Atom* atoms=nullptr; unsigned* marks=nullptr;
    ~Device(){if(atoms)cudaFree(atoms);if(marks)cudaFree(marks);}
};
__device__ double square(double x){return __dmul_rn(x,x);}
__device__ unsigned region(double d,double r,double rw,bool reen){
    return d<=(reen?r:rw)?2u:(reen&&d>r&&d<=rw?1u:0u);
}
__global__ void mark(const Atom* atoms,int ny,int nz,double h,bool reen,bool fkap,
                     bool stern,int ox,int oy,int oz,unsigned* marks){
    const Atom a=atoms[blockIdx.x];
    const int yl=a.y1-a.y0,zl=a.z1-a.z0;
    const size_t count=size_t(a.x1-a.x0)*yl*zl;
    for(size_t q=threadIdx.x;q<count;q+=blockDim.x){
        int x=a.x0+q/(size_t(yl)*zl), y=a.y0+(q/zl)%yl, z=a.z0+q%zl;
        double xc=__dmul_rn(double(x),h),yc=__dmul_rn(double(y),h),zc=__dmul_rn(double(z),h);
        double xs=square(__dsub_rn(xc,a.x)),ys=square(__dsub_rn(yc,a.y)),zs=square(__dsub_rn(zc,a.z));
        double half=__dmul_rn(0.5,h);
        unsigned mask=0;
        if(!fkap){
            double ds=__dadd_rn(__dadd_rn(xs,ys),zs);
            mask=(!reen||stern)?(ds<=a.ri?2u:0u):region(ds,a.r,a.rw,true);
        }
        double dx=square(__dsub_rn(__dadd_rn(xc,half),a.x));
        double dy=square(__dsub_rn(__dadd_rn(yc,half),a.y));
        double dz=square(__dsub_rn(__dadd_rn(zc,half),a.z));
        mask|=region(__dadd_rn(__dadd_rn(dx,ys),zs),a.r,a.rw,reen)<<2;
        mask|=region(__dadd_rn(__dadd_rn(dy,xs),zs),a.r,a.rw,reen)<<4;
        // Original Z expression adds Y before X.
        mask|=region(__dadd_rn(__dadd_rn(dz,ys),xs),a.r,a.rw,reen)<<6;
        if(mask)atomicOr(marks+(size_t(x-ox)*ny+(y-oy))*nz+(z-oz),mask);
    }
}

}
extern "C" int charmm_pb_mstep_cuda(int nx,int ny,int nz,int na,int ns,
    double h,double rw,double ri,double ep,double tx,double ty,double tz,
    double bx,double by,double bz,int reen,int fkap,const int* sel,
    const double* x,const double* y,const double* z,const double* rad,
    float* g,float* gx,float* gy,float* gz,char* err,int errlen){
    auto fail=[&](int code,const char* msg){if(err&&errlen>0)std::snprintf(err,errlen,"%s",msg);return code;};
    if(err&&errlen>0)err[0]=0;
    if(nx<1||ny<1||nz<1||na<0||ns<0||!sel||!x||!y||!z||!rad||!g||!gx||!gy||!gz)
        return fail(2,"Invalid MSTEP dimensions or arrays");
    if(size_t(nx)>size_t(INT_MAX)/ny||size_t(nx)*ny>size_t(INT_MAX)/nz)
        return fail(2,"MSTEP grid exceeds 32-bit indexing");
    if(!(h>0)||!std::isfinite(h)||!std::isfinite(rw)||!std::isfinite(ri)||rw<0||ri<0||
       !std::isfinite(ep)||double(float(ep))!=ep||ep<0||
       !std::isfinite(tx)||!std::isfinite(ty)||!std::isfinite(tz)||
       !std::isfinite(bx)||!std::isfinite(by)||!std::isfinite(bz)||
       (reen!=0&&reen!=1)||(fkap!=0&&fkap!=1))
        return fail(2,"Unsupported MSTEP scalars; use original CPU routine");
    try {
        std::vector<Atom> atoms;atoms.reserve(ns);
        for(int s=0;s<ns;++s){
            int i=sel[s]-1;if(i<0||i>=na)return fail(2,"MSTEP selection out of range");
            double r=rad[i];if(!std::isfinite(r))return fail(2,"Nonfinite MSTEP radius");
            if(r<=0)continue;
            Atom a{};a.x=x[i]+tx-bx;a.y=y[i]+ty-by;a.z=z[i]+tz-bz;
            double rwater=r+rw,rion=r+ri, extent=(rwater+ri)/h;
            double xx=a.x/h,yy=a.y/h,zz=a.z/h;
            // Also bound integer additions in the original Fortran loop setup.
            constexpr double lim=double(INT_MAX/4);
            if(!std::isfinite(extent)||extent>lim||!std::isfinite(xx)||!std::isfinite(yy)||
               !std::isfinite(zz)||std::abs(xx)>lim||std::abs(yy)>lim||std::abs(zz)>lim)
                return fail(2,"Unsupported MSTEP coordinate or extent");
            int nf=int(extent)+2,ix=int(xx)+1,iy=int(yy)+1,iz=int(zz)+1;
            a.x0=std::max(1,ix-nf+1)-1;a.x1=std::min(nx,ix+nf-1);
            a.y0=std::max(1,iy-nf+1)-1;a.y1=std::min(ny,iy+nf-1);
            a.z0=std::max(1,iz-nf+1)-1;a.z1=std::min(nz,iz+nf-1);
            a.r=r*r;a.rw=rwater*rwater;a.ri=rion*rion;
            if(!std::isfinite(a.r)||!std::isfinite(a.rw)||!std::isfinite(a.ri))
                return fail(2,"MSTEP squared radius overflow");
            if(a.x1>a.x0&&a.y1>a.y0&&a.z1>a.z0)atoms.push_back(a);
        }
        if(atoms.empty())return 0;
        int ox=nx,oy=ny,oz=nz,xend=0,yend=0,zend=0;
        for(const Atom& a:atoms){
            ox=std::min(ox,a.x0);oy=std::min(oy,a.y0);oz=std::min(oz,a.z0);
            xend=std::max(xend,a.x1);yend=std::max(yend,a.y1);zend=std::max(zend,a.z1);
        }
        const int sx=xend-ox,sy=yend-oy,sz=zend-oz;
        const size_t points=size_t(sx)*sy*sz;
        Device d;
        #define CUDA_OK(expr) do{cudaError_t status=(expr);if(status!=cudaSuccess)return fail(1,cudaGetErrorString(status));}while(0)
        size_t free=0,total=0;CUDA_OK(cudaMemGetInfo(&free,&total));
        size_t required=points*sizeof(unsigned)+atoms.size()*sizeof(Atom);
        if(free<required+(size_t(2)<<30))return fail(1,"MSTEP requires 2 GiB GPU reserve");
        std::unique_ptr<unsigned[]> staged(new unsigned[points]);
        CUDA_OK(cudaMalloc(&d.atoms,atoms.size()*sizeof(Atom)));
        CUDA_OK(cudaMalloc(&d.marks,points*sizeof(unsigned)));
        CUDA_OK(cudaMemcpy(d.atoms,atoms.data(),atoms.size()*sizeof(Atom),cudaMemcpyHostToDevice));
        CUDA_OK(cudaMemset(d.marks,0,points*sizeof(unsigned)));
        mark<<<unsigned(atoms.size()),128>>>(d.atoms,sy,sz,h,reen,fkap,ri>0,ox,oy,oz,d.marks);
        CUDA_OK(cudaGetLastError());
        CUDA_OK(cudaMemcpy(staged.get(),d.marks,points*sizeof(unsigned),cudaMemcpyDeviceToHost));
        // Only flags are transferred. The downstream CPU preparation needs these
        // host grids, so a single bounded host application avoids four grid uploads
        // and downloads, without repeating any distance or atom-neighbor work.
        float* host[4]={g,gx,gy,gz};
        size_t q=0;
        for(int xx=ox;xx<xend;++xx)for(int yy=oy;yy<yend;++yy){
            size_t p=(size_t(xx)*ny+yy)*nz+oz;
            for(int zz=oz;zz<zend;++zz,++q,++p){
                unsigned mask=staged[q];if(!mask)continue;
                for(int k=0;k<4;++k){
                    unsigned action=(mask>>(2*k))&3u;if(!action)continue;
                    float v=host[k][p],terminal=k?float(ep):0.0f;
                    if(v!=terminal){
                        if(action&2u)host[k][p]=terminal;
                        else if(v>0.0f)host[k][p]=-v;
                    }
                }
            }
        }
        return 0;
        #undef CUDA_OK
    }catch(const std::bad_alloc&){return fail(1,"MSTEP host allocation failed");}
}
