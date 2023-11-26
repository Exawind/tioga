!
! Tioga is a library for overset grid assembly on parallel distributed systems
! Copyright (C) 2015 Jay Sitaraman
!
!
module gridtype
 type grid
    integer :: nv,n4,n5,n6,n8,nwbc,nobc,nvar,ncells,nghost,nmax,ndof
    integer :: mdim(3)
    real*8  :: dx(3),dx0(3)
    real*8,  pointer     :: x(:)=>null(),q(:)=>null(),s(:)=>null(),dq(:)=>null(),xcentroid(:)=>null(),q0(:)=>null()
    integer, allocatable :: bodytag(:)
    integer, allocatable :: iblank(:)
    integer, allocatable :: ghostData(:,:)
    integer, pointer     :: ndc4(:,:)=>null(),ndc5(:,:)=>null(),ndc6(:,:)=>null(),ndc8(:,:)=>null()
    integer, pointer     :: wbcnode(:)=>null(),obcnode(:)=>null()
    real*8,  allocatable :: scal(:)
    real*8,  allocatable :: nres(:),cres(:)
    integer :: nv4 = 4
    integer :: nv5 = 5
    integer :: nv6 = 6
    integer :: nv8 = 8
 end type grid
end module gridtype
