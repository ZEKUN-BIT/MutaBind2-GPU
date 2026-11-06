// Exact integer nonbond pair generation and stable interaction selection.
// No energy, force, floating-point reduction, or CHARMM exclusion is changed.
#include <cuda_runtime.h>
#include <cub/block/block_reduce.cuh>
#include <cub/block/block_scan.cuh>
#include <algorithm>
#include <climits>
#include <cstdio>
#include <cstring>
#include <new>
#include <vector>
namespace {
constexpr size_t chunk_pairs=size_t(16)<<20;
struct Device {
    std::vector<void*> allocations;
    ~Device(){for(void* p:allocations)cudaFree(p);}
    cudaError_t alloc(int** p,size_t n){auto s=cudaMalloc(p,std::max<size_t>(n,1)*sizeof(int));if(s==cudaSuccess)allocations.push_back(*p);return s;}
};
__global__ void dense_pairs(int n,const int* offsets,const int* exrow,const int* ex,
                            const int* posprefix,int first,int base,int* out){
    int i=blockIdx.x+first+1,begin=exrow[i-1],end=exrow[i];
    for(int j=i+1+threadIdx.x;j<=n;j+=blockDim.x){
        int l=begin,r=end;
        while(l<r){int m=l+(r-l)/2;if(abs(ex[m])<j)l=m+1;else r=m;}
        if(l<end&&ex[l]==j)continue;
        int p=offsets[i-1]-base+j-i-1-(posprefix[l]-posprefix[begin]);
        out[p]=(l<end&&ex[l]==-j)?-j:j;
    }
}
__device__ bool keep(int a,int b,int j,const int* sel,const int* other){
    int k=abs(j)-1;return (a==1&&other[k]==1)||(b==1&&sel[k]==1);
}
__global__ void count_selected(int n,const int* rows,const int* pairs,const int* sel,
                               const int* other,int first,int base,int* counts){
    using Reduce=cub::BlockReduce<int,256>;__shared__ typename Reduce::TempStorage tmp;
    int i=blockIdx.x+first,begin=(i?rows[i-1]:0)-base,end=rows[i]-base,sum=0;
    for(int p=begin+threadIdx.x;p<end;p+=blockDim.x)sum+=keep(sel[i],other[i],pairs[p],sel,other);
    int total=Reduce(tmp).Sum(sum);if(threadIdx.x==0)counts[blockIdx.x]=total;
}
__global__ void select_pairs(int n,const int* rows,const int* pairs,const int* sel,
                             const int* other,const int* destrows,int first,int base,int* out){
    using Scan=cub::BlockScan<int,256>;__shared__ typename Scan::TempStorage tmp;
    int i=blockIdx.x+first,begin=(i?rows[i-1]:0)-base,end=rows[i]-base;
    int dest=blockIdx.x?destrows[blockIdx.x-1]:0;
    for(int base=begin;base<end;base+=blockDim.x){
        int p=base+threadIdx.x,j=p<end?pairs[p]:0;
        int flag=p<end&&keep(sel[i],other[i],j,sel,other),offset,total;
        Scan(tmp).ExclusiveSum(flag,offset,total);
        if(flag)out[dest+offset]=j;
        dest+=total;__syncthreads();
    }
}
}
extern "C" int charmm_vdw_dense_cuda(int n,int capacity,const int* iblo,const int* inb14,
    int* rows,int* pairs,int* count,char* err,int len){
    if(err&&len>0)std::memset(err,0,len);
    auto fail=[&](int code,const char* msg){if(err&&len>0)std::snprintf(err,len,"%s",msg);return code;};
    if(n<1||capacity<0||!iblo||!inb14||!rows||!pairs||!count)return fail(2,"Invalid dense pair arguments");
    try{
        std::vector<int> offsets(n+1),exrows(n+1),ex,prefix(1,0);
        int prev=0;long long total=0;
        for(int i=1;i<=n;++i){
            if(iblo[i-1]<prev)return fail(2,"Nonmonotone exclusion offsets");
            int last=0,positive=0;
            for(int q=prev;q<iblo[i-1];++q){
                int v=inb14[q];if(v==INT_MIN)return fail(2,"Invalid exclusion");
                int j=abs(v);if(j<=last||j>n)return fail(2,"Unsorted exclusion row");last=j;
                if(j<=i)continue;
                ex.push_back(v);prefix.push_back(prefix.back()+(v>0));positive+=(v>0);
            }
            prev=iblo[i-1];exrows[i]=ex.size();total+=n-i-positive;
            if(total>INT_MAX)return fail(2,"Pair count exceeds 32-bit indexing");offsets[i]=int(total);
        }
        *count=int(total);if(total>capacity)return 3;
        if(!total){std::fill(rows,rows+n,0);return 0;}
        Device d;int *doff,*dexrow,*dex,*dpre,*dout;
        #define CHECK(expr) do{auto status=(expr);if(status!=cudaSuccess){if(err&&len>0)std::snprintf(err,len,"%s: %s",#expr,cudaGetErrorString(status));return 1;}}while(0)
        size_t free,totalmem;CHECK(cudaMemGetInfo(&free,&totalmem));
        size_t buffer=std::min<size_t>(total,chunk_pairs);
        size_t need=sizeof(int)*(buffer+offsets.size()+exrows.size()+ex.size()+prefix.size());
        if(free<need+(size_t(2)<<30))return fail(1,"Dense pairs require 2 GiB GPU reserve");
        CHECK(d.alloc(&doff,offsets.size()));CHECK(d.alloc(&dexrow,exrows.size()));
        CHECK(d.alloc(&dex,ex.size()));CHECK(d.alloc(&dpre,prefix.size()));CHECK(d.alloc(&dout,buffer));
        CHECK(cudaMemcpy(doff,offsets.data(),offsets.size()*sizeof(int),cudaMemcpyHostToDevice));
        CHECK(cudaMemcpy(dexrow,exrows.data(),exrows.size()*sizeof(int),cudaMemcpyHostToDevice));
        if(!ex.empty())CHECK(cudaMemcpy(dex,ex.data(),ex.size()*sizeof(int),cudaMemcpyHostToDevice));
        CHECK(cudaMemcpy(dpre,prefix.data(),prefix.size()*sizeof(int),cudaMemcpyHostToDevice));
        for(int first=0;first<n;){
            int last=first+1,base=offsets[first];
            while(last<n&&size_t(offsets[last+1]-base)<=buffer)++last;
            size_t amount=offsets[last]-base;
            if(amount>buffer)return fail(2,"Dense row exceeds bounded buffer");
            if(amount){
                dense_pairs<<<last-first,256>>>(n,doff,dexrow,dex,dpre,first,base,dout);CHECK(cudaGetLastError());
                // Transfer only this contiguous group of rows. Errors are fatal;
                // the native subprocess never accepts a partially generated list.
                CHECK(cudaMemcpy(pairs+base,dout,amount*sizeof(int),cudaMemcpyDeviceToHost));
            }
            first=last;
        }
        std::copy(offsets.begin()+1,offsets.end(),rows);return 0;
        #undef CHECK
    }catch(const std::bad_alloc&){return fail(1,"Dense-pair allocation failed");}
}
extern "C" int charmm_vdw_select_cuda(int n,int capacity,const int* rows,const int* pairs,
    const int* sel,const int* other,int* destrows,int* dest,char* err,int len){
    if(err&&len>0)std::memset(err,0,len);
    auto fail=[&](int code,const char* msg){if(err&&len>0)std::snprintf(err,len,"%s",msg);return code;};
    if(n<1||capacity<0||!rows||!pairs||!sel||!other||!destrows||!dest)return fail(2,"Invalid selection arguments");
    int prev=0;for(int i=0;i<n;++i){if(rows[i]<prev)return fail(2,"Nonmonotone pair rows");prev=rows[i];}
    int total=rows[n-1];if(!total){std::fill(destrows,destrows+n,0);return 0;}
    try{
        Device d;int *dr,*dp,*ds,*dt,*dc,*dout;std::vector<int> counts(n);
        #define CHECK(expr) do{auto status=(expr);if(status!=cudaSuccess){if(err&&len>0)std::snprintf(err,len,"%s: %s",#expr,cudaGetErrorString(status));return 1;}}while(0)
        size_t free,totalmem;CHECK(cudaMemGetInfo(&free,&totalmem));
        size_t buffer=std::min<size_t>(total,chunk_pairs);
        size_t need=sizeof(int)*(2*buffer+size_t(n)*4);
        if(free<need+(size_t(2)<<30))return fail(1,"Pair selection requires 2 GiB GPU reserve");
        CHECK(d.alloc(&dr,n));CHECK(d.alloc(&dp,buffer));CHECK(d.alloc(&dout,buffer));CHECK(d.alloc(&ds,n));CHECK(d.alloc(&dt,n));CHECK(d.alloc(&dc,n));
        CHECK(cudaMemcpy(dr,rows,size_t(n)*sizeof(int),cudaMemcpyHostToDevice));
        CHECK(cudaMemcpy(ds,sel,size_t(n)*sizeof(int),cudaMemcpyHostToDevice));
        CHECK(cudaMemcpy(dt,other,size_t(n)*sizeof(int),cudaMemcpyHostToDevice));
        int selected=0;
        std::vector<int> resultrows(n);
        for(int first=0;first<n;){
            int base=first?rows[first-1]:0,last=first+1;
            while(last<n&&size_t(rows[last]-base)<=buffer)++last;
            size_t amount=rows[last-1]-base;
            if(amount>buffer)return fail(2,"Selection row exceeds bounded buffer");
            if(amount){
                CHECK(cudaMemcpy(dp,pairs+base,amount*sizeof(int),cudaMemcpyHostToDevice));
                count_selected<<<last-first,256>>>(n,dr,dp,ds,dt,first,base,dc);CHECK(cudaGetLastError());
                CHECK(cudaMemcpy(counts.data(),dc,size_t(last-first)*sizeof(int),cudaMemcpyDeviceToHost));
                int local=0;for(int k=0;k<last-first;++k){
                    local+=counts[k];counts[k]=local;
                    if(size_t(selected)+local>size_t(capacity))return fail(2,"Selected pairs exceed capacity");
                    resultrows[first+k]=selected+local;
                }
                if(local){
                    CHECK(cudaMemcpy(dc,counts.data(),size_t(last-first)*sizeof(int),cudaMemcpyHostToDevice));
                    select_pairs<<<last-first,256>>>(n,dr,dp,ds,dt,dc,first,base,dout);CHECK(cudaGetLastError());
                    CHECK(cudaMemcpy(dest+selected,dout,size_t(local)*sizeof(int),cudaMemcpyDeviceToHost));
                }
                selected+=local;
            }else std::fill(resultrows.begin()+first,resultrows.begin()+last,selected);
            first=last;
        }
        std::copy(resultrows.begin(),resultrows.end(),destrows);return 0;
        #undef CHECK
    }catch(const std::bad_alloc&){return fail(1,"Pair-selection allocation failed");}
}
