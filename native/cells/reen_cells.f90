! Conservative spatial candidates for MAYER_REEN; no distance/SEGPAR arithmetic.
! Width is four times the largest positive expanded radius. Nonadjacent cells
! are at least that far apart; two intersecting spheres are <= twice that radius
! apart. The factor-two gap covers bin arithmetic rounding on the bridge's
! existing finite coordinate/radius domain. Selected indices remain in order.
module mb_reen_cell_index
  use iso_fortran_env, only: real64, int64
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
contains
  subroutine build_reen_cells(n,selected,x,y,z,radius,water,offsets,candidates,atom_cell,ok)
    integer,intent(in) :: n,selected(n)
    real(real64),intent(in) :: x(*),y(*),z(*),radius(*),water
    integer,allocatable,intent(out) :: offsets(:),candidates(:),atom_cell(:)
    logical,intent(out) :: ok
    integer,allocatable :: xyz(:,:),population(:),lengths(:)
    integer :: dims(3),i,j,a,c,k,ix,iy,iz,jx,jy,jz,neighbor,status,cursor,occupied
    integer(int64) :: ncells,entries
    real(real64) :: lower(3),upper(3),point(3),extent(3),largest,width
    ok=.false.
    if(n<256.or..not.ieee_is_finite(water).or.water<0) return
    lower=huge(0.0_real64);upper=-huge(0.0_real64);largest=0
    do j=1,n
      a=selected(j)
      if(a<1) return
      point=[x(a),y(a),z(a)]
      if(.not.all(ieee_is_finite(point)).or.any(abs(point)>1.0e6_real64)) return
      if(.not.ieee_is_finite(radius(a)+water)) return
      if(radius(a)+water<1.0e-4_real64.or.radius(a)+water>1.0e4_real64) return
      lower=min(lower,point);upper=max(upper,point)
      largest=max(largest,radius(a)+water)
    enddo
    width=4.0_real64*largest
    extent=(upper-lower)/width
    if(any(extent>999999.0_real64)) return
    dims=int(extent)+1
    ncells=int(dims(1),int64)*int(dims(2),int64)
    if(ncells>1000000_int64) return
    ncells=ncells*int(dims(3),int64)
    if(ncells>1000000_int64) return
    allocate(xyz(3,n),population(ncells),lengths(ncells),atom_cell(n),stat=status)
    if(status/=0) return
    population=0;lengths=0
    do j=1,n
      a=selected(j)
      xyz(:,j)=min(dims-1,max(0,int(([x(a),y(a),z(a)]-lower)/width)))
      c=1+xyz(1,j)+dims(1)*(xyz(2,j)+dims(2)*xyz(3,j))
      atom_cell(j)=c;population(c)=population(c)+1
    enddo
    occupied=count(population>0)
    ! Avoid paying for many sparsely occupied cells; the caller retains its
    ! original scan if this representation is not economical.
    if(occupied>n/8) return
    entries=0
    do iz=0,dims(3)-1
      do iy=0,dims(2)-1
        do ix=0,dims(1)-1
          c=1+ix+dims(1)*(iy+dims(2)*iz)
          if(population(c)==0) cycle
          k=0
          do jz=max(0,iz-1),min(dims(3)-1,iz+1)
            do jy=max(0,iy-1),min(dims(2)-1,iy+1)
              do jx=max(0,ix-1),min(dims(1)-1,ix+1)
                neighbor=1+jx+dims(1)*(jy+dims(2)*jz)
                k=k+population(neighbor)
              enddo
            enddo
          enddo
          lengths(c)=k;entries=entries+int(k,int64)
        enddo
      enddo
    enddo
    if(entries>16000000_int64.or.entries<1) return
    allocate(offsets(ncells+1),candidates(entries),stat=status)
    if(status/=0) return
    offsets(1)=1
    do c=1,int(ncells)
      offsets(c+1)=offsets(c)+lengths(c)
    enddo
    do i=1,n
      c=atom_cell(i)
      if(lengths(c)==0) cycle
      cursor=offsets(c)
      do j=1,n
        if(any(abs(xyz(:,j)-xyz(:,i))>1)) cycle
        candidates(cursor)=j;cursor=cursor+1
      enddo
      if(cursor/=offsets(c+1)) return
      lengths(c)=0
    enddo
    ok=.true.
  end subroutine build_reen_cells
end module mb_reen_cell_index
