!--------------------------------------------------------------
! Solver for quasilinear DAEs of the following form
!
! DAEQ1 :  a(t,u)u' = h(t,u)
!          f(t,u)   = 0
!          u(t0) = u0, such that f(t0,u0) = 0
!
!          dim(rge a)= dim(rge h) = nd, dim(rge f) = na, 
!          nu = nd na, nw = 0, dim(u) = nu,
!
!          requires subroutine daefct with 
!            fname = 'amt ' for evaluating a(t,u)
!            fname = 'hf  ' for evaluating h(t,u)
!            fname = 'ff  ' for evaluating f(t,u)
!            fname = 'dff ' for evaluating df(t,u)
!
! DAEN1 :  u' = v 
!          f(t,u,v,w) = 0,
!          u(t0) = u0, v(t0) = up0, w(t0) = w0, 
!                      such that f(t0,u0,up0,w0) = 0
!
!          dim(u) = nu, dim(w) = nw,
!          nd = nu, na = nu + nw, dim(rge f)= na
!
!          requires subroutine daefct with 
!            fname = 'ff  ' for evaluating f(t,u)
!            fname = 'dff ' for evaluating df(t,u)
!
! DAEQ2    a(t,u)u' + b(t,u)w = h(t,u)
!          f(t,u)             = 0
!          u(t0) = u0, such that f(t0,u0) = 0
!
!          dim(rge a)= dim(rge h) = nd, dim(rge f)= na, 
!          dim(w) = nw, dim(u) = nu, nu + nw = nd + na
!
!          requires subroutine daefct with 
!            fname = 'amt ' for evaluating a(t,u)
!            fname = 'bmt ' for evaluating b(t,u)
!            fname = 'hf  ' for evaluating h(t,u)
!            fname = 'ff  ' for evaluating f(t,u)
!            fname = 'dff ' for evaluating df(t,u)
!
! NOHOL    a(u) u" + d_{u'}f(t,u,u')^T w = h(t,u,u')
!          f(t,u,u')                     = 0
!          u(t0) = u0, u'(t0) = up0
!                      such that f(t0,u0,up0) = 0
!
!          dim(rge a)= dim(rge h) = nd, dim(rge f)= na, 
!          dim(w) = na, dim(u) = nd,
!
!          requires subroutine daefct with 
!            fname = 'amt ' for evaluating a(u)
!            fname = 'hf  ' for evaluating h(t,u,u')
!            fname = 'ff  ' for evaluating f(t,u,u')
!            fname = 'dff ' for evaluating df(t,u,u')
!
! DAEQ3    a(t,u,u')u" + b(t,u,u')w = h(t,u,u')
!          f(t,u)                   = 0
!          u(t0)=u0, u'(t0)=up0, such that f(t0,u0) = 0
!
!          dim(rge f)= na, dim(w) = nw, 
!          dim(rge a)= dim(rge h) = nd, dim(u) = nd+na-nw
!
!          requires subroutine daefct with 
!            fname = 'amt ' for evaluating a(t,u,u')
!            fname = 'bmt ' for evaluating b(t,u,u')
!            fname = 'hf  ' for evaluating h(t,u,u')
!            fname = 'ff  ' for evaluating f(t,u)
!            fname = 'dff ' for evaluating df(t,u)
!            fname = 'd2ff' for evaluating d2f(t,u)(v,v)
!
! EULAG    a(u) u" + df(u)^T w = h(u,u')
!          f(u)                = 0
!          u(t0) = u0, such that f(u0) = 0
!
!          dim(rge a) = dim(rge h) = nd, dim(rge f) = na
!          dim(w) = na, dim(u) = nd
!
!          requires subroutine daefct with 
!            fname = 'amt ' for evaluating a(u)
!            fname = 'hf  ' for evaluating h(u,u')
!            fname = 'ff  ' for evaluating f(u)
!            fname = 'dff ' for evaluating df(u)
!            fname = 'd2ff' for evaluating d2f(u)(v,v)
!
! Version:  W. Rheinboldt, Univ. of Pittsburgh, Nov. 6, 2000
!
! The file contains the modules:
!
! module dae_dat      to set global data for the dae solver
! module dae_wrk      to define work arrays and to provide
!                     associated routines for setting and
!                     deleting these arrays
! module dae_man      to define local parametrizations on the
!                     manifolds and to provide the associated
!                     routines for establishing such bases
!                     and for the evauation of points on
!                     the manifolds
! module dae_loc      to provide subroutines for evaluating the
!                     local odes for the various dae types
! subroutine daeslv   global routine for solving the dae
!
!--------------------------------------------------------------
!
module dae_dat
  !
  !--------------------------------------------------------
  ! Module for the global data of the dae solver package
  ! and the following routine:
  !
  ! subroutine setup         for initializing the run
  !                          and setting dimensions
  ! subroutine close_run     for ending the run
  !--------------------------------------------------------
  !
  use def_kinds, ONLY : RP => RPX
  use file_sup
  implicit none
  !
  character(len=5) :: dnam  ! type of DAE
  !
  integer, save :: nd    ! number of differential equations
  integer, save :: na    ! number of algebraic equations
  integer, save :: nu    ! number of differential variables
  integer, save :: nw    ! number of algebraic variables
  !
  integer, save :: ka    ! ka+1 = na
  integer, save :: km    ! km+1 = manifold dimension
  integer, save :: ko    ! ko+1 = local ode dimension
  integer, save :: ks    ! ks+1 = ambient space dimension
  integer, save :: kx    ! kx+1 = global dimension of x
  !
  integer, save :: outlun                ! logical output unit
  !  
  real(kind=RP), dimension(:), pointer :: x  ! solution point
  !
  ! counters for statistics
  !
  integer, save :: nstep,naccpt,nrejct,nbas,nff,ndff
  !
  ! global error tolerances
  !
  real(kind=RP), parameter :: epmach = epsilon(1.0_RP)
  real(kind=RP), save      :: hitol  = 0.0_RP, &
                              abstol = 0.0_RP, &
                              reltol = 0.0_RP
  !
  ! global test constants
  !
  integer,       parameter :: maxit   = 50, &
                              kfailmx = 3
  real(kind=RP), parameter :: alphopt = 0.5_RP
  integer, parameter       :: itstpmx = 20, &
                              maxstp  = 30
  !
  ! global variables for step control
  !
  integer,       save :: nmax = 100000
  real(kind=RP), save :: atol,rtol
  real(kind=RP), save :: tout,posneg
  real(kind=RP), save :: h,hmin,hmax,hint
  ! 
contains
  !
  subroutine setup(dnm,pname,outunit,nvd,nva,nvu,nvw,info)
    !
    implicit none
    intrinsic log, ceiling, min
    !
    character(len=5), intent(in)  :: dnm    
    character(len=6), intent(in)  :: pname
    integer,          intent(in)  :: outunit
    integer,          intent(in)  :: nvd,nva,nvu,nvw
    integer,          intent(out) :: info
    !
    type(file)    :: output
    real(kind=RP) :: tmp
    !
    !-------------------------------------------------
    ! Opens the files for the log and the output,
    ! initializes the statistics counters, and sets
    ! the dimensions for the different DAEs
    ! Error returns:
    !  info =  0  no error
    !  info = -1  unknown type of DAE
    !  info = -2  dimension error
    !--------------------------------------------------
    !
    ! open the log file and the output file
    !
    outlun = outunit
    call open_log(pname)
    call open_file(output,outlun,'REPLACE',trim(pname)//'.out')
    !
    ! initialize data for the ode solver
    !
    hmin = reltol
    hint = 0.0_RP
    !
    ! initialize stat counters
    !
    nstep  = 0           ! total number of steps
    naccpt = 0           ! number of accepted steps
    nrejct = 0           ! number of rejected steps
    nbas   = 0           ! number of basis evaluations
    nff    = 0           ! number of calls of function f
    ndff   = 0           ! number of calls of derivative df 
    !
    ! set dimensions
    !
    info = 0
    dnam = trim(dnm)
    !
    select case(dnam)
      !
      case('daen1')
        !
        ! given nu and nw set na = nu + nw, nd = nu
        !
        nu = nvu
        nw = nvw
        na = nu + nw
        nd = nu
        if(na /= nva .or. nd /= nvu) then
          info = -2
          return
        endif
        ka = na-1
        km = nu
        ko = km
        ks = nu+nu+nw
        kx = ks
        !
      case('daeq1')
        !
        ! given nd and na set nw = 0, nu = nd+na
        !
        nd = nvd
        na = nva
        nw = 0
        nu = nd + na
        if(nw /= nvw .or. nu /= nvu) then
          info = -2
          return
        endif
        ka = na-1
        km = nd
        ko = km
        ks = nu
        kx = nu+nu 
        !
      case('daeq2')
        !
        ! given nd, na, and nu set nw = nd + na - nu
        !
        nd = nvd
        na = nva
        nu = nvu 
        nw = nd+na-nu
        if(nw /= nvw) then
          info = -2
          return
        endif
        ka = na-1
        km = nu-na
        ko = km
        ks = nu
        kx = nu+nu+nw
        !
      case('nohol')
        !
        ! given nd and na set nw = na, nu = nd
        !
        nd = nvd
        na = nva
        nw = na
        nu = nd
        if(nw /= nvw .or. nu /= nvu) then
          info = -2
          return
        endif
        ka = na-1
        km = nu+nu-na
        ko = km+km
        ks = nu+nu
        kx = nu+nu+nu+nw
        !
      case('daeq3')
        !
        ! given nd, na, and nu, set nw = nd + na - nu
        ! 
        nd = nvd
        na = nva
        nu = nvu
        nw = nd+na-nu
        if(nw /= nvw) then
          info = -2
          return
        endif
        ka = na-1
        km = nu-na
        ko = km+km
        ks = nu
        kx = nu+nu+nw 
        !
      case('eulag')
        !
        ! given nd and na set nw = na, nu = nd
        !
        nd = nvd
        na = nva
        nw = na
        nu = nd
        if(nw /= nvw .or. nu /= nvu) then
          info = -2
          return
        endif
        ka = na-1
        km = nu-na
        ko = km+km
        ks = nu
        kx = nu+nu+nu+nw 
        !
      case default
        !
        info = -1
        !
    end select
    !
    if(associated(x)) nullify(x)       ! allocate solution point
    allocate(x(0:kx))
    !
    ! write header
    !
    write(outlun,*)'problem name: '//pname
    write(outlun,*)'DAE type    : '//dnam
    write(outlun,10)nd,na,nu,nw
    10 format(' nv= ',I2,' na= ',I2,' nu= ',I2,' nw= ',I2)
    !
  end subroutine setup
  !
  subroutine close_run(info)
    !
    implicit none
    integer, intent(in) :: info
    !
    write(outlun,*)'Run completed with info= ',info
    call log_msg('main',info,'Run completed')
    !
    ! write statistics
    !  
    write(outlun,*)'Statistics'
    write(outlun,*)'total number of steps = ',nstep
    write(outlun,*)'       accepted steps = ',naccpt
    write(outlun,*)'       rejected steps = ',nrejct
    write(outlun,*)'local bases computed  = ',nbas
    write(outlun,*)'function evaluations  = ',nff
    write(outlun,*)'jacobian evaluations  = ',ndff
    !
    call close_file()
    nullify(x)
    !
  end subroutine close_run 
  !
end module dae_dat
!
module dae_wrk
  !
  !----------------------------------------------------------
  ! module for defining the work arrays for the dae solver
  ! and for the following routines:
  ! 
  ! subroutine loc_work        for allocating the arrays
  ! subroutine del_work        for disallocating the arrays
  !----------------------------------------------------------
  !
  use def_kinds, ONLY : RP => RPX
  implicit none
  !
  integer,       dimension(:),   allocatable, target :: jwrk
  real(kind=RP), dimension(:),   allocatable, target :: wkxc
  real(kind=RP), dimension(:),   allocatable, target :: wkxn
  real(kind=RP), dimension(:),   allocatable, target :: wkxi
  real(kind=RP), dimension(:),   allocatable, target :: wky
  real(kind=RP), dimension(:),   allocatable, target :: wkzp
  real(kind=RP), dimension(:),   allocatable, target :: wkzi
  real(kind=RP), dimension(:),   allocatable, target :: wrk1
  real(kind=RP), dimension(:),   allocatable, target :: wrk2
  real(kind=RP), dimension(:,:), allocatable, target :: wrkmt1
  real(kind=RP), dimension(:,:), allocatable, target :: wrkmt2
  real(kind=RP), dimension(:,:), allocatable, target :: wrkmt3
  real(kind=RP), dimension(:),   allocatable, target :: wk0
  real(kind=RP), dimension(:),   allocatable, target :: wk1
  real(kind=RP), dimension(:),   allocatable, target :: wk2
  real(kind=RP), dimension(:),   allocatable, target :: wk3
  real(kind=RP), dimension(:),   allocatable, target :: wk4
  real(kind=RP), dimension(:),   allocatable, target :: wk5
  real(kind=RP), dimension(:),   allocatable, target :: wk6
  real(kind=RP), dimension(:,:), allocatable, target :: wkint
  !
contains
  !
  subroutine loc_work(klge,ksml,info)
    !
    implicit none
    !
    integer, intent(in)  :: klge,ksml
    integer, intent(out) :: info
    !---------------------------------------------------
    ! Allocates the work arrays
    ! Error returns:
    !  info =  0  no error
    !  info = -3  allocation error
    !---------------------------------------------------
    ! 
    integer :: ialloc
    !
    info = 0
    if(allocated(jwrk))   deallocate(jwrk)
    if(allocated(wkxc))   deallocate(wkxc)
    if(allocated(wkxn))   deallocate(wkxn)
    if(allocated(wkxi))   deallocate(wkxi)
    if(allocated(wky))    deallocate(wky)
    if(allocated(wkzp))   deallocate(wkzp)
    if(allocated(wkzi))   deallocate(wkzi)
    if(allocated(wrk1))   deallocate(wrk1)
    if(allocated(wrk2))   deallocate(wrk2)
    if(allocated(wrkmt1)) deallocate(wrkmt1)
    if(allocated(wrkmt2)) deallocate(wrkmt2)
    if(allocated(wrkmt2)) deallocate(wrkmt3)
    if(allocated(wk0))    deallocate(wk0)
    if(allocated(wk1))    deallocate(wk1)
    if(allocated(wk2))    deallocate(wk2)
    if(allocated(wk3))    deallocate(wk3)
    if(allocated(wk4))    deallocate(wk4)
    if(allocated(wk5))    deallocate(wk5)
    if(allocated(wk6))    deallocate(wk6)
    if(allocated(wkint))  deallocate(wkint)
    ! 
    allocate(jwrk(0:klge),STAT=ialloc)
    allocate(wkxc(0:klge),STAT=ialloc)
    allocate(wkxn(0:klge),STAT=ialloc)
    allocate(wkxi(0:klge),STAT=ialloc)
    allocate(wky(0:klge), STAT=ialloc)
    allocate(wkzp(0:klge),STAT=ialloc)
    allocate(wkzi(0:klge),STAT=ialloc)
    allocate(wrk1(0:klge),STAT=ialloc)
    allocate(wrk2(0:klge),STAT=ialloc)
    !
    allocate(wrkmt1(0:klge,0:klge),STAT=ialloc)
    allocate(wrkmt2(0:klge,0:klge),STAT=ialloc)
    allocate(wrkmt3(0:klge,0:klge),STAT=ialloc)
    !
    allocate(wk0(0:ksml),STAT=ialloc)
    allocate(wk1(0:ksml),STAT=ialloc)
    allocate(wk2(0:ksml),STAT=ialloc)
    allocate(wk3(0:ksml),STAT=ialloc)
    allocate(wk4(0:ksml),STAT=ialloc)
    allocate(wk5(0:ksml),STAT=ialloc)
    allocate(wk6(0:ksml),STAT=ialloc)
    !
    allocate(wkint(0:4,0:ksml),STAT=ialloc)
    ! 
    if(ialloc /= 0) then
      info = -3
      return
    endif
    !
  end subroutine loc_work

  subroutine del_work
    !  
    implicit none
    !
    if(allocated(jwrk))   deallocate(jwrk)
    if(allocated(wkxc))   deallocate(wkxc)
    if(allocated(wkxn))   deallocate(wkxn)
    if(allocated(wkxi))   deallocate(wkxi)
    if(allocated(wky))    deallocate(wky)
    if(allocated(wkzp))   deallocate(wkzp)
    if(allocated(wkzi))   deallocate(wkzi)
    if(allocated(wrk1))   deallocate(wrk1)
    if(allocated(wrk2))   deallocate(wrk2)
    if(allocated(wrkmt1)) deallocate(wrkmt1)
    if(allocated(wrkmt2)) deallocate(wrkmt2)
    if(allocated(wrkmt3)) deallocate(wrkmt3)
    if(allocated(wk0))    deallocate(wk0)
    if(allocated(wk1))    deallocate(wk1)
    if(allocated(wk2))    deallocate(wk2)
    if(allocated(wk3))    deallocate(wk3)
    if(allocated(wk4))    deallocate(wk4)
    if(allocated(wk5))    deallocate(wk5)
    if(allocated(wk6))    deallocate(wk6)
    if(allocated(wkint))  deallocate(wkint)
    !
  end subroutine del_work
  !
end module dae_wrk

module dae_man
  !
  !----------------------------------------------------------
  ! Module for defining the basis of local parametrizations
  ! on the constraint manifolds and including the routines:
  !
  ! subroutine loc_basis      for allocating storage for 
  !                           a basis
  ! subroutine del_basis      for disallocating the basis 
  !                           storage
  ! subroutine new_basis      for computing a first basis
  !                           at a given global point
  ! subroutine reval_basis    for computing a new basis 
  !                           and orienting it in 
  !                           accordance with an earlier 
  !                           computed basis
  ! subroutine eval_point     for evaluating a global
  !                           point corresponding to a 
  !                           point with given local
  !                           coordinates in a given basis
  !----------------------------------------------------------
  !
  use def_kinds, ONLY : RP => RPX
  use dae_dat
  !
  interface
    subroutine daefct(fname,xc,vc,vec,mat,iret)
      !
      use def_kinds, ONLY : RP => RPX
      implicit none
      character(len=4),                intent(in)  :: fname
      real(kind=RP), dimension(0:),    intent(in)  :: xc
      real(kind=RP), dimension(0:),    intent(in)  :: vc
      real(kind=RP), dimension(0:),    intent(out) :: vec
      real(kind=RP), dimension(0:,0:), intent(out) :: mat
      integer,                         intent(out) :: iret
      optional vc,vec,mat
    end subroutine daefct
  end interface
  !
  type basis
    ! ambient space dimension
    integer                                :: namb
    !
    ! dimension of manifold
    integer                                :: mdim
    !
    ! base point
    real(kind=RP), dimension(:),   pointer :: baspt
    !
    ! matrix of basis vectors
    real(kind=RP), dimension(:,:), pointer :: basm
    !
    ! augmented matrix
    real(kind=RP), dimension(:,:), pointer :: augb
    !
    ! pivot array
    integer, dimension(:), pointer         :: jaugb
  end type basis
  !
contains
  !
  subroutine loc_basis(xbs,kbs,kbm)
    !
    integer,     intent(in)    :: kbs,kbm
    type(basis), intent(inout) :: xbs
    !
    xbs%namb = kbs
    xbs%mdim = kbm
    !
    if(associated(xbs%baspt)) nullify(xbs%baspt)
    if(associated(xbs%basm)) nullify(xbs%basm)
    if(associated(xbs%augb)) nullify(xbs%augb)
    if(associated(xbs%jaugb)) nullify(xbs%jaugb)
    allocate(xbs%baspt(0:kbs))
    allocate(xbs%basm(0:kbs,0:kbm))
    allocate(xbs%augb(0:kbs,0:kbs))
    allocate(xbs%jaugb(0:kbs))
    !
  end subroutine loc_basis
  !
  subroutine del_basis(xbs)
    !
    implicit none
    type(basis), intent(inout) :: xbs
    !
    xbs%namb = 0
    xbs%mdim = 0
    nullify(xbs%baspt,xbs%basm,xbs%augb,xbs%jaugb)
    !
  end subroutine del_basis
  !
  subroutine new_basis(xbs,x,info)
    !
    use dae_wrk
    use file_sup
    use lin_lib, ONLY : housl, lqf, luf
    implicit none
    !
    type(basis),                  intent(inout) :: xbs
    real(kind=RP), dimension(0:), intent(in)    :: x
    integer,                      intent(out)   :: info
    !
    ! pointers to work arrays
    !
    integer,       dimension(:),   pointer :: jpivt
    real(kind=RP), dimension(:),   pointer :: tau
    real(kind=RP), dimension(:),   pointer :: cnrm
    real(kind=RP), dimension(:,:), pointer :: dfm
    real(kind=RP), dimension(:,:), pointer :: xbasis
    !
    ! local variables and parameters
    !
    integer                  :: i,j,kbs,kbm,kba
    real(kind=RP), parameter :: zer=0.0_RP, one=1.0_RP
    !
    !--------------------------------------------------------
    ! The routine assumes that the pointer arrays for the
    ! basis_class and the targeted work arrays have been
    ! allocated by means of the subroutines loc_basis 
    ! and loc_work. It then evaluates a new basis at
    ! the global point x.
    ! Error returns:
    !   info =  0 successful return
    !   info = -5 zero pivot in basis computation
    !   info = -6 zero pivot in augmented matrix
    !--------------------------------------------------------
    !
    kbs = xbs%namb                       ! dimensions
    kbm = xbs%mdim
    kba = kbs - kbm - 1
    !
    dfm    => wrkmt1                      ! associate
    jpivt  => jwrk                        ! pointers
    tau    => wrk1                        ! with
    cnrm   => wrk2                        ! work
    xbasis => wrkmt2                      ! arrays
    !
    ! get jacobian and retain copy
    !
    call daefct(fname='dff ',xc=x(0:kbs), &
                mat=dfm(0:kba,0:kbs),iret=info)
    xbs%augb(0:kba,0:kbs) = dfm(0:kba,0:kbs)
    !
    ! compute lq-factorization of submatrix dfm(0:kba,1:kbs)
    !        
    call lqf(dfm(0:kba,1:kbs),jpivt(0:kba),tau(0:kba), &
                                       cnrm(0:kba),info)
    if(info /= 0) then
      info = -5
      call log_msg('new_basis',info,'zero pivot in &
                                &basis computation')
      nullify(dfm,jpivt,tau,cnrm,xbasis)
      return
    endif
    !
    xbasis(0:kbs,0:kbm) = zer             ! initialize the basis
    xbasis(0,0) = one
    if(kbm > 0) then
      do j = 1,kbm,1
        xbasis(kba+1+j,j) = one
      enddo
      !
      do j = kba,0,-1
        call housl(xbasis(j+1:kbs,1:kbm),dfm(j,j+1:kbs),tau(j),info)
      enddo
    endif
    !
    xbs%baspt(0:kbs) = x(0:kbs)             ! store basis in xbs
    xbs%basm(0:kbs,0:kbm) = xbasis(0:kbs,0:kbm)
    nullify(dfm,jpivt,tau,cnrm,xbasis)
    !
    ! set up and factor the augmented matrix
    !
    xbs%augb(kba+1:kba+1+kbm,0:kbs) = transpose(xbs%basm(0:kbs,0:kbm))
    call luf(xbs%augb(0:kbs,0:kbs),xbs%jaugb(0:kbs),info)
    if(info /= 0) then
      info = -6
      call log_msg('new_basis',info,'zero pivot in augmented matrix')
      return
    endif
    info = 0
    !
  end subroutine new_basis
  !
  subroutine reval_basis(xbs,x,info)
    !
    use dae_wrk
    use file_sup
    use lin_lib, ONLY : housl, lqf, luf
    implicit none
    intrinsic abs, sign
    !
    type(basis),                  intent(inout) :: xbs
    real(kind=RP), dimension(0:), intent(in)    :: x
    integer,                      intent(inout) :: info
    !
    ! pointers to work arrays
    !
    integer,       dimension(:),   pointer :: jpivt
    real(kind=RP), dimension(:),   pointer :: tau
    real(kind=RP), dimension(:),   pointer :: cnrm
    real(kind=RP), dimension(:,:), pointer :: dfm
    real(kind=RP), dimension(:,:), pointer :: prod
    real(kind=RP), dimension(:,:), pointer :: xbasis
    !
    ! local variables and parameters
    !
    integer                  :: kbs,kbm,kba
    integer                  :: i,j,jtmp,k,maxi
    real(kind=RP)            :: bsgn,esgn,tmp,emax
    real(kind=RP), parameter :: zer=0.0_RP, one=1.0_RP
    !
    !--------------------------------------------------------
    ! The routine assumes that the pointer arrays for the
    ! basis_class and the target work arrays have been
    ! allocated by means of the subroutines loc_basis 
    ! and loc_work. It then evaluates a basis and reorients
    ! it in accordance with the basis given in xbs. The
    ! new basis then replaces the old one in xbs.
    !
    ! The return indicator has the values
    !   info =  0 successful return
    !   info = -5 zero pivot in basis computation
    !   info = -6 zero pivot in augmented matrix
    !   info = -7 basis orientation failed
    !--------------------------------------------------------
    !
    kbs = xbs%namb                            ! dimensions
    kbm = xbs%mdim
    kba = kbs - kbm - 1
    !
    dfm    => wrkmt1                          ! associate
    jpivt  => jwrk                            ! pointers
    tau    => wrk1                            ! with 
    cnrm   => wrk2                            ! work
    xbasis => wrkmt2                          ! arrays
    !
    ! get jacobian and retain copy
    !
    call daefct(fname='dff ',xc=x(0:kbs), &
                mat=dfm(0:kba,0:kbs),iret=info)
    xbs%augb(0:kba,0:kbs) = dfm(0:kba,0:kbs)
    !
    ! compute lq-factorization of submatrix dfm(0:kba,1:kbs)
    !                                             
    call lqf(dfm(0:kba,1:kbs),jpivt(0:kba),tau(0:kba), &
                                       cnrm(0:kba),info)
    if(info /= 0) then
      info = -5
      call log_msg('reval_basis',info,'zero pivot &
                             &in basis computation')
      nullify(dfm,jpivt,tau,cnrm,xbasis)
      return
    endif
    !
    xbasis(0:kbs,0:kbm) = zer             ! initialize the basis
    xbasis(0,0) = one
    if(kbm > 0) then
      do j = 1,kbm,1
        xbasis(kba+1+j,j) = one
      enddo
      !
      do j = kba,0,-1
        call housl(xbasis(j+1:kbs,1:kbm),dfm(j,j+1:kbs),tau(j),info)
      enddo
    endif
    !
    ! reorient the basis 
    !
    prod => wrkmt1
    prod(0:kbm,0:kbm) = matmul(transpose(xbs%basm(0:kbs,0:kbm)),&
                                             xbasis(0:kbs,0:kbm))
    !
    if(kbm == 0) then
      bsgn = sign(one,prod(0,0))
    else      
      bsgn = one
      !
      ! in each row of the product matrix determine
      ! the largest, uncovered element
      !
      jpivt(0:kbm) = 0                   ! begin with zero vector              
      do i = 0,kbm,1
        emax = zer
        maxi = -1
        esgn = one
        do j = 0,kbm,1
          if(jpivt(j) == 0) then
            tmp = prod(i,j)
            if(abs(tmp) > emax) then
              emax = abs(tmp)
              maxi = j
              esgn = sign(one,tmp)
            endif
          endif
        enddo
        if(maxi >= 0) then
          jpivt(maxi) = i
        else
          info = -7
          call log_msg('reval_basis',info,'basis &
                             &orientation failed')
          exit                        ! reorientation failed
        endif
        if(esgn < zer) xbasis(0:kbs,maxi) = -xbasis(0:kbs,maxi)
      enddo
      !
      if(info == 0) then              ! if the reorientation
        !                             ! was successful then
        do j = 0,kbm,1                ! reorient the basis
          i = jwrk(j)
          if(j /= i) then
            do k = 0,kbs
              tmp = xbasis(k,i)
              xbasis(k,i) = xbasis(k,j)
              xbasis(k,j) = tmp
            enddo
            jtmp    = jpivt(i)
            jpivt(i) = jpivt(j)
            jpivt(j) = jtmp
          endif
        enddo
        !
      endif
      !
    endif
    !
    xbs%baspt(0:kbs) = x(0:kbs)                  ! store the basis
    xbs%basm(0:kbs,0:kbm) = bsgn*xbasis(0:kbs,0:kbm)
    !
    nullify(jpivt,xbasis,prod)
    !
    ! insert basis in augmented matrix and factor the matrix
    !
    xbs%augb(kba+1:kba+1+kbm,0:kbs) = transpose(xbs%basm(0:kbs,0:kbm))
    call luf(xbs%augb(0:kbs,0:kbs),xbs%jaugb(0:kbs),info)
    if(info /= 0) then
      info = -6
      call log_msg('reval_basis',info,'zero pivot in &
                                    &augmented matrix')
      return
    endif
    !
  end subroutine reval_basis
  !
  subroutine eval_point(xbs,y,xc,itstp,alpha,info)
    !
    use def_kinds, ONLY : RP => RPX
    use dae_dat
    use dae_wrk
    use file_sup
    use lin_lib, ONLY : lus
    implicit none
    intrinsic abs,ceiling,min,maxval,sqrt,log
    !
    type(basis),                  intent(in)  :: xbs    
    real(kind=RP), dimension(0:), intent(in)  :: y
    real(kind=RP), dimension(0:), intent(out) :: xc
    integer,                      intent(out) :: itstp,info
    real(kind=RP),                intent(out) :: alpha
    !
    real(kind=RP), dimension(:), pointer :: fxc   ! pointers
    real(kind=RP), dimension(:), pointer :: xc0   ! to work
    !
    ! parameters and local variables
    !
    real(kind=RP), parameter :: testhi   = 1.5_RP, &
                                testlo   = 0.95_RP, &
                                testerr  = 1.0e3_RP
    !
    integer       :: kbs,kbm,kba,kfail,kpred,i
    real(kind=RP) :: fnrm,stepl,fnrml,step,tolx
    real(kind=RP) :: ritst,xfact,cdist
    real(kind=RP) :: fnquot,stquot,quot,err
    !
    real(kind=RP), parameter :: zer = 0.0_RP, one = 1.0_RP
    !
    !-----------------------------------------------------------
    ! Routine applies chord newton iteration to evaluate the
    ! global point x on the manifold defined by the given local
    ! coordinates y in R^md.
    !
    ! There are three convergence indicators:
    !  itstp    counter of the number of steps taken 
    !  alpha    estimated contraction constant
    !
    ! The return indicator has the values:
    !  info =  0    successful iteration
    !  info = -18   the iteration matrix has a zero pivot
    !  info = -19   maximum number of iteration steps
    !  info = -20-k iteration diverged at step k
    !------------------------------------------------------------
    kbs = xbs%namb                      ! set
    kbm = xbs%mdim                      ! array
    kba = kbs - kbm - 1                 ! indices
    itstp = 0                           ! initialize
    alpha = zer                         ! data
    ritst = zer
    info = 0
    !
    if(maxval(abs(y(0:kbm))) < epmach) then   ! return
      xc(0:kbs) = xbs%baspt(0:kbs)            ! when y
      return                                  ! is small
    endif
    xfact = maxval(abs(xbs%baspt(0:kbs)))
    if(xfact > one) xfact = one
    !
    fxc => wrk1                               ! associate pointer
    xc0 => wrk2                               ! with work arrays
    !
    ! get starting point of iteration
    !
    xc(0:kbs) = xbs%baspt(0:kbs) + &
                 matmul(xbs%basm(0:kbs,0:kbm),y(0:kbm))
    !
    xc0(0:kbs) = xc(0:kbs)                    ! retain in xc0
    call daefct(fname='ff  ',xc=xc(0:kbs), &
                  vec=fxc(0:kba),iret=info)
    !
    fnrm = maxval(abs(fxc(0:kba)))
    if(fnrm <= hitol) then           ! quick exit for small fnrm
      nullify(fxc,xc0)
      return
    endif
    !
    step  = one
    kfail  = 0
    info = 0
    !
    do while (itstp <= maxstp)            ! iteration loop
      !
      if(itstp == maxstp) then
        kfail = 2
        exit
      endif
      !
      fxc(kba+1:kbs) = zer                ! augment function vector
      !
      call lus(xbs%augb(0:kbs,0:kbs),xbs%jaugb(0:kbs), &
               fxc(0:kbs),info)        ! solve chord newton equation
      if(info /= 0) then
        info = -18
        call log_msg('eval_point',info,'zero pivot in &
                                   &newton iteration')
        itstp = -itstp
        nullify(fxc,xc0)
        return
      endif
      !
      stepl = step
      step = maxval(abs(fxc(0:kbs)))              ! step norm
      !
      xc(0:kbs) = xc(0:kbs) - fxc(0:kbs)          ! next iterate
      !
      call daefct(fname='ff  ',xc=xc(0:kbs), &
                     vec=fxc(0:kba),iret=info)    ! next residual
      !
      fnrml = fnrm
      fnrm   = maxval(abs(fxc(0:kba)))            ! residual norm
      !
      cdist = maxval(abs(xc(0:kbs)-xc0(0:kbs)))
      quot = one
      if(cdist > zer)then
        quot = step/cdist
        if(ritst > zer .and. quot > zer) alpha = quot**(one/ritst)
      endif
      itstp = itstp + 1
      ritst = ritst + one
      err = max(step,reltol)
      kpred = itstp
      if(alpha > zer .and. alpha < one) then
        err = (alpha/(one-alpha))*err
        quot = xfact*reltol/err
        if(quot > zer .and. quot < one) quot = log(quot)/log(alpha)
        kpred = itstp + ceiling(quot)
      endif
      !
      ! convergence and divergence tests
      !
      select case(itstp)
        !
        case(1)
          if((fnrm <= hitol) .or. (step <= xfact*hitol)) exit
          if(fnrm > testhi*fnrml .or. step > testhi*stepl) then
            kfail = 1
            exit
          endif
        case(2:5)
          if(((fnrm <= abstol) .and. (step <= xfact*abstol)) .or. &
             ((err <= xfact*reltol) .and. (fnrm <= reltol))) exit
          if(fnrm > testlo*fnrml .or. step > testlo*stepl) then
            kfail = 1
            exit
          endif
        case(6:maxstp)
          if(((fnrm <= abstol) .and. (step <= xfact*abstol)) .or. &
             (err <= xfact*reltol)) exit
          if(fnrm > testlo*fnrml .or. step > testlo*stepl .or. &
             kpred > maxstp .or. alpha >= one) then
            kfail = 1
            exit
          endif
        case default
          kfail = 2
          exit
      end select
      !
    enddo                          ! end of iteration loop
    !
    select case(kfail)
      !
      case(1)
        info = -20-itstp
        call log_msg('eval_point',info,'divergence detected')
      case(2)
        info = -19
        call log_msg('eval_point',info,'maximum number of &
                                        &iteration steps')
    end select
    !
    nullify(fxc,xc0)               ! disconnect pointer assoc.    
    !
  end subroutine eval_point
  !
end module dae_man
!
module dae_loc
  !
  !------------------------------------------------------------
  ! Module for computing the local odes corresponding to a
  ! dae of the given type and for the following routines:
  !
  ! subroutine locode    generic interface routine to handle
  !                      the required polymorphism
  ! subroutine dyn1      for computing the local ode for
  !                      daes of the type daen1
  ! subroutine dyq1      for computing the local ode for
  !                      daes of the type daeq1
  ! subroutine dyq2      for computing the local ode for
  !                      daes of the type daeq2
  ! subroutine dynh      for computing the local ode for
  !                      daes of the type nohol
  ! subroutine dyq3      for computing the local ode for
  !                      daes of the type daeq3
  ! subroutine dyel      for computing the local ode for
  !                      daes of the type eulag
  !
  ! Error returns in the subroutines:
  !   info =  0  no error
  !   info = -1  unknown type of DAE
  !   info = -8  zero pvot in augmented matrix
  !   info = -9  zero pivot in local direction
  !------------------------------------------------------------
  !
  use def_kinds, ONLY : RP => RPX
  use file_sup
  use lin_lib, ONLY : luf, lus
  use dae_dat
  use dae_wrk
  use dae_man
  !   
  ! pointers to work arrays
  !
  integer,       dimension(:),   pointer :: jzm
  real(kind=RP), dimension(:),   pointer :: wq1
  real(kind=RP), dimension(:),   pointer :: wq2
  real(kind=RP), dimension(:,:), pointer :: zmat
  real(kind=RP), dimension(:,:), pointer :: emat
  real(kind=RP), dimension(:,:), pointer :: dphi
  !
  real(kind=RP), parameter :: zer=0.0_RP, one=1.0_RP
  !
contains
  ! 
  subroutine locode(dnam,xbs,y,xp,zp,info)
    !
    implicit none
    !
    character(len=5),             intent(in)    :: dnam
    type(basis),                  intent(in)    :: xbs
    real(kind=RP), dimension(0:), intent(inout) :: y
    real(kind=RP), dimension(0:), intent(inout) :: xp
    real(kind=RP), dimension(0:), intent(inout) :: zp
    integer,                      intent(inout) :: info
    !
    select case(dnam)
      !
      case ('daen1')
        !
        call dyn1(xbs,xp,zp,info)
        !
      case ('daeq1')
        !
        call dyq1(xbs,xp,zp,info)
        ! 
      case ('daeq2')
        !
        call dyq2(xbs,xp,zp,info)
        ! 
      case ('nohol')
        !
        call dynh(xbs,y,xp,zp,info)
        ! 
      case ('daeq3')
        !
        call dyq3(xbs,y,xp,zp,info)
        !
      case ('eulag')
        !
        call dyel(xbs,y,xp,zp,info)
        !
      case default
        !
        info = -1
        call log_msg('dae_loc',info,'wrong DAE type')
        return
        !
    end select
    !
  end subroutine locode
  !
  subroutine dyn1(xbs,xp,zp,info)
    !
    !-------------------------------------------------------
    ! Local ode for daen1
    !
    ! Evaluate y' as solution of the system
    !
    !   dphi_1 y' = phi_2
    !
    ! where phi denotes the local parametrization and 
    ! phi(y) = (t,u,u'w) = (phi_0,phi_1,phi_2,phi_3).
    !-------------------------------------------------------
    !
    implicit none
    !
    type(basis),                  intent(in)    :: xbs
    real(kind=RP), dimension(0:), intent(inout) :: xp
    real(kind=RP), dimension(0:), intent(inout) :: zp
    integer,                      intent(inout) :: info
    !
    integer :: j  
    !------------------------------------------------------
    !
    jzm  => jwrk                       ! associate 
    zmat => wrkmt1                     ! pointers
    dphi => wrkmt2                     ! with work arrays
    !
    ! evaluate and factor the current augmented matrix
    !
    call daefct(fname='dff ',xc=xp(0:ks), &
                mat=zmat(0:ka,0:ks),iret=info)
    zmat(na:na+km,0:ks) = transpose(xbs%basm(0:ks,0:km))
    !
    call luf(zmat(0:ks,0:ks),jzm(0:ks),info)
    if(info /= 0) then
      info = -8
      call log_msg('eval_field',info,'zero pivot in &
                                  &augmented matrix')
      nullify(jzm,zmat,dphi)
      return
    endif
    !
    ! initialize the columns of dphi and solve for dphi
    !
    dphi(0:ks,0:km) = zer
    do j = 0,km,1
      dphi(na+j,j) = one
    enddo
    call lus(zmat(0:ks,0:ks),jzm(0:ks),dphi(0:ks,0:km),info)
    !
    ! evaluate the local direction
    !
    zmat(0:km,0:km) = dphi(0:km,0:km)
    !
    call luf(zmat(0:km,0:km),jzm(0:km),info)
    if(info /= 0) then
      info = -9
      call log_msg('eval_field',info,'zero pivot in &
                                   &local direction')
      nullify(jzm,zmat,dphi)
      return
    endif
    !
    ! set up and solve for the local direction
    !
    zp(0)    = one
    zp(1:nu) = xp(km+1:km+nu)
    !
    call lus(zmat(0:km,0:km),jzm(0:km),zp(0:km),info)
    !
    nullify(jzm,zmat,dphi)
    info = 0
    !
  end subroutine dyn1
  !
  subroutine dyq1(xbs,xp,zp,info)
    !
    !-------------------------------------------------------
    ! Routine for evaluating the local direction for DAEs
    ! of the form daeq1.
    !-------------------------------------------------------
    !
    implicit none
    !
    type(basis),                  intent(in)    :: xbs
    real(kind=RP), dimension(0:), intent(inout) :: xp
    real(kind=RP), dimension(0:), intent(inout) :: zp
    integer,                      intent(inout) :: info
    !
    integer :: j 
    !------------------------------------------------------
    !
    jzm  => jwrk                       ! associate 
    zmat => wrkmt1                     ! pointers
    dphi => wrkmt2                     ! with work arrays
    !
    ! evaluate and factor current augm. matrix
    ! we use that na+km=ks
    !
    call daefct(fname='dff ',xc=xp(0:ks), &
                mat=zmat(0:ka,0:ks),iret=info)
    zmat(na:ks,0:ks) = transpose(xbs%basm(0:ks,0:km))
    !
    call luf(zmat(0:ks,0:ks),jzm(0:ks),info)
    if(info /= 0) then
      info = -8
      call log_msg('eval_field',info,'zero pivot in &
                                  &augmented matrix')
      nullify(jzm,zmat,dphi)
      return
    endif
    !
    ! initialize the columns of dphi and solve for dphi
    !
    dphi(0:ks,0:km) = zer
    do j = 0,km,1
      dphi(na+j,j) = one
    enddo
    call lus(zmat(0:ks,0:ks),jzm(0:ks),dphi(0:ks,0:km),info)
    !
    ! evaluate the local direction -- first get amt
    ! and multiply with dphi
    !
    zmat(0:ks,0:ks)     = zer
    call daefct(fname='amt ',xc=xp(0:ks), &
                mat=zmat(1:nd,1:nu),iret=info)
    zmat(0,0)       = one
    zmat(0:ks,0:km) = matmul(zmat(0:ks,0:nu),dphi(0:nu,0:km))
    !
    ! factor the matrix
    !
    call luf(zmat(0:km,0:km),jzm(0:km),info)
    if(info /= 0) then
      info = -9
      call log_msg('eval_field',info,'zero pivot in &
                                   &local direction')
      nullify(jzm,zmat,dphi)
      return
    endif
    !
    ! Evaluate the right side and solve for local direction
    !
    zp(0) = one
    call daefct(fname='hf  ',xc=xp(0:ks), &
                  vec=zp(1:nd),iret=info)
    call lus(zmat(0:km,0:km),jzm(0:km),zp(0:km),info)
    !
    ! Evaluate the global derivative and store in xp
    !
    xp(nu+1:nu+nu) = matmul(dphi(1:nu,0:km),zp(0:km))
    !
    nullify(jzm,zmat,dphi)
    info = 0
    !
  end subroutine dyq1
  !
  subroutine dyq2(xbs,xp,zp,info)   
    !
    !---------------------------------------------
    ! Routine for evaluating the local direction
    ! for DAEs of the form daeq2
    !---------------------------------------------
    !
    implicit none
    !
    type(basis),                   intent(in)    :: xbs
    real(kind=RP), dimension(0:),  intent(inout) :: xp
    real(kind=RP), dimension(0:),  intent(inout) :: zp
    integer,                       intent(inout) :: info
    !
    integer :: i,j
    !------------------------------------------------------
    !
    jzm  => jwrk                       ! associate 
    zmat => wrkmt1                     ! pointers
    dphi => wrkmt2                     ! with work arrays
    !
    ! evaluate and factor the current augmented matrix
    ! we use that ks = na + km
    !
    call daefct(fname='dff ',xc=xp(0:ks), &
                mat=zmat(0:ka,0:ks),iret=info)
    zmat(na:ks,0:ks) = transpose(xbs%basm(0:ks,0:km))
    !
    call luf(zmat(0:ks,0:ks),jzm(0:ks),info)
    if(info /= 0) then
      info = -8
      call log_msg('locode',info,'zero pivot in &
                               &augmented matrix')
      nullify(jzm,zmat,dphi)
      return
    endif
    !
    ! initialize the columns of dphi and solve for dphi
    !
    dphi(0:ks,0:km)  = zer
    do j = 0,km,1
      dphi(na+j,j) = one
    enddo
    call lus(zmat(0:ks,0:ks),jzm(0:ks),dphi(0:ks,0:km),info)
    !
    ! evaluate the local direction 
    ! first evaluate amt and multiply by dphi
    !
    zmat(0:ks,0:ks)     = zer
    call daefct(fname='amt ',xc=xp(0:ks), &
           mat=zmat(1:nd,1:nu),iret=info)          ! get amt
    zmat(0,0)           = one
    zmat(0:ks,0:km) = matmul(zmat(0:ks,0:nu),dphi(0:nu,0:km))
    !
    call daefct(fname='bmt ',xc=xp(0:ks), &
             mat=zmat(1:nd,km+1:km+nw),iret=info)  ! get bmt
    !
    ! factor the matrix
    !
    call luf(zmat(0:ks,0:ks),jzm(0:ks),info)
    if(info /= 0) then
      info = -9
      call log_msg('locode',info,'zero pivot in &
                                &local direction')
      nullify(jzm,zmat,dphi)
      return
    endif
    !
    ! Evaluate the right side and solve for local direction
    !
    zp(0)      = one
    call daefct(fname='hf  ',xc=xp(0:ks), &
                  vec=zp(1:nd),iret=info)         ! get hf
    call lus(zmat(0:ks,0:ks),jzm(0:ks),zp(0:ks),info)
    !
    ! Evaluate the global derivative and store in xp
    !
    xp(nu+1:nu+nu) = matmul(dphi(1:nu,0:km),zp(0:km))
    !
    ! Store algebraic variable
    !
    xp(nu+nu+1:nu+nu+nw) = zp(km+1:km+nw)
    !
    nullify(jzm,zmat,dphi)
    info = 0
    !
  end subroutine dyq2
  !  
  subroutine dynh(xbs,y,xp,zp,info)
    !
    !-------------------------------------------------------
    ! Routine for evaluating the local direction for DAEs 
    ! of the form eulag
    ! 
    ! For given y = (t,y1,y2) the routine evaluates 
    ! y' = (y1', y2') = (y2, z) and w by as solution of 
    ! the linear system
    !
    !   ( (1      0) 0    ) ( z ) = ( 1           )
    !   ( (0 A*dphi) Df^T ) ( w )   ( h - A*d2phi )
    !
    ! where phi denotes the local parametrization and the 
    ! arguments of the functions A, Df^T are (phi(y1)),
    ! while that of h is (phi(y1),dphi(y1)).
    !-------------------------------------------------------
    !
    implicit none
    !
    type(basis),                   intent(in)    :: xbs
    real(kind=RP), dimension(0:),  intent(inout) :: y
    real(kind=RP), dimension(0:),  intent(inout) :: xp
    real(kind=RP), dimension(0:),  intent(inout) :: zp
    integer,                       intent(inout) :: info
    !
    integer :: i,j,istart
    !------------------------------------------------------
    !
    jzm  => jwrk                       ! associate 
    zmat => wrkmt1                     ! pointers
    dphi => wrkmt2                     ! with work arrays
    emat => wrkmt3
    !
    ! evaluate and factor the current augmented matrix
    ! we use that na+km=ks
    !
    call daefct(fname='dff ',xc=xp(0:ks), &
                mat=zmat(0:ka,0:ks),iret=info)
    emat(0:nu,km+1:ks)    = zer
    emat(nu+1:ks,km+1:ks) = transpose(zmat(0:ka,nu+1:ks))
    zmat(na:ks,0:ks) = transpose(xbs%basm(0:ks,0:km))
    !
    call luf(zmat(0:ks,0:ks),jzm(0:ks),info)
    if(info /= 0) then
      info = -8
      call log_msg('locode',info,'zero pivot in &
                               &augmented matrix')
      nullify(jzm,zmat,dphi,emat)
      return
    endif
    !
    ! initialize the columns of dphi and solve for dphi
    !
    dphi(0:ks,0:km) = zer
    do j = 0,km,1
      dphi(na+j,j) = one
    enddo
    call lus(zmat(0:ks,0:ks),jzm(0:ks),dphi(0:ks,0:km),info)
    !
    ! evaluate the local direction 
    ! first evaluate amt and multiply by dphi
    !
    zmat(0:ks,0:ks) = zer
    do i=0,nu,1
      zmat(i,i) = one
    enddo
    call daefct(fname='amt ',xc=xp(0:ks), &
               mat=zmat(nu+1:ks,nu+1:ks),iret=info)    ! get amt
    emat(0:ks,0:km) = matmul(zmat(0:ks,0:ks),dphi(0:ks,0:km))
    !
    ! factor the matrix
    !
    call luf(emat(0:ks,0:ks),jzm(0:ks),info)
    if(info /= 0) then
      info = -9
      call log_msg('locode',info,'zero pivot in &
                                &local direction')
      nullify(jzm,zmat,dphi,emat)
      return
    endif
    !
    ! Evaluate right side and solve for local direction
    ! We use that nu+nu=ks and nu=nd
    !
    zp(0)    = one
    zp(1:nu) = xp(nu+1:ks)
    call daefct(fname='hf  ',xc=xp(0:ks), &
                  vec=zp(nu+1:ks),iret=info)         ! get hf
    call lus(emat(0:ks,0:ks),jzm(0:ks),zp(0:ks),info)
    !
    ! Evaluate the global derivative and store in xp
    !
    xp(ks+1:ks+nu) = matmul(dphi(nu+1:nu+nu,0:km),zp(0:km))
    !
    ! Store algebraic variable using that kx=ks+nw
    !
    xp(ks+nu+1:kx) = zp(km+1:km+nw)
    !
    nullify(jzm,zmat,dphi,emat)
    info = 0
    !
  end subroutine dynh
  !
  subroutine dyq3(xbs,y,xp,zp,info)
    !
    !-------------------------------------------------------
    ! Routine for evaluating the local direction for DAEs 
    ! of the form daeq2.
    ! 
    ! For given y = (t,y1,y2) the routine evaluates 
    ! y' = (y1', y2') = (y2, z) and w by as solution of 
    ! the linear system
    !
    !   ( (1      0) 0 ) ( z ) = ( 1           )
    !   ( (0 A*dphi) B ) ( w )   ( h - A*d2phi )
    !
    ! where phi denotes the local parametrization and the 
    ! arguments of the functions A, B, and h are 
    ! (phi(y1),dphi(y1)).
    !-------------------------------------------------------
    !
    implicit none
    !
    type(basis),                   intent(in)    :: xbs
    real(kind=RP), dimension(0:),  intent(inout) :: y
    real(kind=RP), dimension(0:),  intent(inout) :: xp
    real(kind=RP), dimension(0:),  intent(inout) :: zp
    integer,                       intent(inout) :: info
    !
    integer :: i,j,istart
    !------------------------------------------------------
    !
    jzm  => jwrk                       ! associate 
    wq1  => wrk1                       ! pointers
    wq2  => wrk2                       ! with
    zmat => wrkmt1                     ! work
    dphi => wrkmt2                     ! arrays
    !
    istart = 0
    if(info == 1) istart = 1
    !
    ! evaluate and factor the current augmented matrix
    !
    call daefct(fname='dff ',xc=xp(0:ks), &
                mat=zmat(0:ka,0:ks),iret=info)
    zmat(na:na+km,0:ks) = transpose(xbs%basm(0:ks,0:km))
    !
    call luf(zmat(0:ks,0:ks),jzm(0:ks),info)
    if(info /= 0) then
      info = -8
      call log_msg('locode',info,'zero pivot in &
                              &augmented matrix')
      nullify(jzm,zmat,dphi)
      return
    endif
    !
    ! initialize the columns of dphi and solve for dphi
    !
    dphi(0:ks,0:km)  = zer
    do j = 0,km,1
      dphi(na+j,j) = one
    enddo
    call lus(zmat(0:ks,0:ks),jzm(0:ks),dphi(0:ks,0:km),info)
    !
    ! evaluate first local derivative
    !
    if(istart == 1) then
      wq1(0)      = one 
      wq1(1:nu) = xp(nu+1:nu+nu)
      zp(0:km)    = matmul(transpose(xbs%basm(0:nu,0:km)), &
                                                wq1(0:nu))
      y(km+1:km+km) = zp(1:km)
    else
      zp(0)          = one
      zp(1:km)       = y(km+1:km+km)
      wq1(0:nu)      = matmul(dphi(0:nu,0:km),zp(0:km))
      xp(nu+1:nu+nu) = wq1(1:nu)
    endif
    ! 
    ! evaluate d2f and with it d2phi
    !
    call daefct(fname='d2ff',xc=xp(0:ks), &
           vc=wq1(0:ks),vec=wq2(0:ka),iret=info)
    wq2(0:ka)  = -wq2(0:ka)
    wq2(na:ks) = zer
    call lus(zmat(0:ks,0:ks),jzm(0:ks),wq2(0:ks),info)
    wq2(0) = zer
    !
    ! evaluate the local direction 
    ! first evaluate amt and multiply by dphi
    !
    zmat(0:nd,0:nd) = zer
    call daefct(fname='amt ',xc=xp(0:kx), &
           mat=zmat(1:nd,1:nu),iret=info)        ! get amt
    zmat(0,0) = one
    !
    ! evaluate A*d2dphi and A*dphi
    !
    wq1(1:nd)       = matmul(zmat(1:nd,1:ks),wq2(1:ks))
    zmat(0:nd,0:km) = matmul(zmat(0:nd,0:ks),dphi(0:ks,0:km))
    !
    call daefct(fname='bmt ',xc=xp(0:kx), &
           mat=zmat(1:nd,km+1:km+nw),iret=info)  ! get bmt
    !
    ! factor the matrix
    !
    call luf(zmat(0:nd,0:nd),jzm(0:nd),info)
    if(info /= 0) then
      info = -9
      call log_msg('locode',info,'zero pivot in &
                                &local direction')
      nullify(jzm,zmat,dphi)
      return
    endif
    !
    ! Evaluate the right side and solve for local direction
    !
    call daefct(fname='hf  ',xc=xp(0:kx), &
                vec=wq2(1:nd),iret=info)           ! get hf
    wq2(1:nd) = wq2(1:nd) - wq1(1:nd)
    wq2(0)    = zer
    call lus(zmat(0:nd,0:nd),jzm(0:nd),wq2(0:nd),info)
    zp(km+1:km+km) = wq2(1:km)
    !
    ! Store algebraic variable
    !
    xp(nu+nu+1:nu+nu+nw) = wq2(km+1:km+nw)
    !
    nullify(jzm,zmat,dphi)
    info = 0
    !
  end subroutine dyq3
  !
  subroutine dyel(xbs,y,xp,zp,info)
    !
    !-------------------------------------------------------
    ! Routine for evaluating the local direction for DAEs 
    ! of the form eulag
    ! 
    ! For given y = (t,y1,y2) the routine evaluates 
    ! y' = (y1', y2') = (y2, z) and w by as solution of 
    ! the linear system
    !
    !   ( (1      0) 0    ) ( z ) = ( 1           )
    !   ( (0 A*dphi) Df^T ) ( w )   ( h - A*d2phi )
    !
    ! where phi denotes the local parametrization and the 
    ! arguments of the functions A, Df^T are (phi(y1)),
    ! while that of h is (phi(y1),dphi(y1)).
    !-------------------------------------------------------
    !
    implicit none
    !
    type(basis),                   intent(in)    :: xbs
    real(kind=RP), dimension(0:),  intent(inout) :: y
    real(kind=RP), dimension(0:),  intent(inout) :: xp
    real(kind=RP), dimension(0:),  intent(inout) :: zp
    integer,                       intent(inout) :: info
    !
    integer :: i,j,istart
    !------------------------------------------------------
    !
    jzm  => jwrk                       ! associate 
    wq1  => wrk1                       ! pointers
    wq2  => wrk2                       ! with
    zmat => wrkmt1                     ! work
    dphi => wrkmt2                     ! arrays
    emat => wrkmt3
    !
    istart = 0
    if(info == 1) istart = 1
    !
    ! evaluate and factor the current augmented matrix
    !
    call daefct(fname='dff ',xc=xp(0:ks), &
                mat=zmat(0:ka,0:ks),iret=info)
    emat(0:nd,0:nd)   = zer
    emat(1:nd,km+1:nd)= transpose(zmat(0:ka,1:nd))
    zmat(na:ks,0:ks)  = transpose(xbs%basm(0:ks,0:km))
    !
    call luf(zmat(0:ks,0:ks),jzm(0:ks),info)
    if(info /= 0) then
      info = -8
      call log_msg('locode',info,'zero pivot in &
                              &augmented matrix')
      nullify(jzm,zmat,emat,dphi)
      return
    endif
    !
    ! initialize the columns of dphi and solve for dphi
    !
    dphi(0:ks,0:km)  = zer
    do j = 0,km,1
      dphi(na+j,j) = one
    enddo
    call lus(zmat(0:ks,0:ks),jzm(0:ks),dphi(0:ks,0:km),info)
    !
    ! evaluate first local derivative
    !
    if(istart == 1) then
      wq1(0)    = one 
      wq1(1:nu) = xp(nu+1:nu+nu)
      zp(0:km)  = matmul(transpose(xbs%basm(0:nu,0:km)), &
                                              wq1(0:nu))
      y(km+1:km+km) = zp(1:km)
    else
      zp(0)          = one
      zp(1:km)       = y(km+1:km+km)
      wq1(0:nu)      = matmul(dphi(0:nu,0:km),zp(0:km))
      xp(nu+1:nu+nu) = wq1(1:nu)
    endif
    ! 
    ! evaluate d2f and with it d2phi
    !
    call daefct(fname='d2ff',xc=xp(0:ks), &
           vc=wq1(0:ks),vec=wq2(0:ka),iret=info)
    wq2(0:ka)  = -wq2(0:ka)
    wq2(na:ks) = zer
    call lus(zmat(0:ks,0:ks),jzm(0:ks),wq2(0:ks),info)
    wq2(0) = zer
    !
    ! evaluate the local direction 
    ! first evaluate amt and multiply by dphi
    !
    call daefct(fname='amt ',xc=xp(0:ks), &
           mat=zmat(1:nd,1:nu),iret=info)        ! get amt
    !
    ! evaluate A*d2dphi and A*dphi
    !
    wq1(1:nd)       = matmul(zmat(1:nd,1:ks),wq2(1:ks))
    emat(1:nd,1:km) = matmul(zmat(1:nd,1:ks),dphi(1:ks,1:km))
    emat(0,0) = one
    !
    ! factor the matrix
    !
    call luf(emat(0:nd,0:nd),jzm(0:nd),info)
    if(info /= 0) then
      info = -9
      call log_msg('locode',info,'zero pivot in &
                                &local direction')
      nullify(jzm,zmat,emat,dphi)
      return
    endif
    !
    ! Evaluate the right side and solve for local direction
    !
    call daefct(fname='hf  ',xc=xp(0:kx), &
                vec=wq2(1:nd),iret=info)           ! get hf
    wq2(1:nd) = wq2(1:nd) - wq1(1:nd)
    wq2(0)      = zer
    call lus(emat(0:nd,0:nd),jzm(0:nd),wq2(0:nd),info)
    zp(km+1:km+km) = wq2(1:km)
    !
    ! Store accelerations and algebraic variables
    !
    xp(nu+nu+1:nu+nu+nu) = wq1(1:nu) &
                + matmul(dphi(1:nu,0:km),wq2(0:km))
    xp(nu+nu+nu+1:nu+nu+nu+nw) = wq2(km+1:nd)
    !
    nullify(jzm,zmat,emat,dphi)
    info = 0
    !
  end subroutine dyel
  !
end module dae_loc

subroutine daeslv(info)
  !
  use def_kinds, ONLY : RP => RPX
  use file_sup
  use dae_dat
  use dae_loc
  use dae_wrk
  use dae_man
  implicit none
  !
  integer, intent(out) :: info
  !
  !--------------------------------------------------------------
  ! solver for daes of type daen1, daeq1, daeq2, nohol,
  ! daeq3, eulag and for the following local subroutines:
  !
  ! subroutine loc_storage  for allocating storage for the solver
  ! subroutine del_storage  for disallocating this storage
  ! subroutine outp         for handling all output at computed 
  !                         points and performing interpolatioj
  !                         on the ode output when required
  ! subroutine dopstn       Runge-Kutta ode solver
  !
  ! Error returns
  !   info =   0  no error
  !   info = -10  ode step-count exceeded
  !   info = -11  ode step fell below given minimum
  !   info = -12  ode routine failed
  !--------------------------------------------------------------
  !
  intrinsic sqrt
  !
  real(kind=RP), dimension(:),   pointer :: xc    ! pointers
  real(kind=RP), dimension(:),   pointer :: xn    !
  real(kind=RP), dimension(:),   pointer :: xi    ! to
  real(kind=RP), dimension(:),   pointer :: zp    !
  real(kind=RP), dimension(:),   pointer :: zi    ! work
  real(kind=RP), dimension(:),   pointer :: y     !
  real(kind=RP), dimension(:,:), pointer :: yint  ! arrays
  !
  ! local variables
  !
  integer                :: i,j
  character(LEN=6), save :: task,ptyp
  type(basis),      save :: xbs
  integer,          save :: itstp,kbas,kfail
  logical,          save :: presup, presw,interp,last
  real(kind=RP),    save :: hold,tc,tprev,tloc,tnext,tprnt
  real(kind=RP),    save :: alpha
  !
  ! parameters
  !
  real(kind=RP), parameter :: hfact  = 0.1_RP, &
                              tfact  = 1.01_RP
  !  
  !----------------------------------------------------------
  !
  call loc_storage()
  !
  ! check data
  !
  if(reltol <= zer) reltol = sqrt(epmach)
  if(abstol <= zer) abstol = reltol*1.0e-2_RP
  if(hitol  <= zer) hitol  = abstol*1.0e-1_RP
  rtol = reltol
  atol = rtol*0.1_RP
  !
  interp = .false.
  if(hint /= zer) interp = .true.
  !
  posneg = one
  if(tout < x(0)) posneg = -posneg
  if(hmin <= zer) hmin = reltol
  if(hmax == zer) hmax = abs(tout - x(0))
  if(hmax < zer) hmax = -hmax
  if(posneg*h < zer) h = -h
  if(h == zer) h = posneg*sqrt(real(ko+1,kind=RP))
  !
  ! set up xc
  !
  xc(0:kx) = x(0:kx)
  !
  ! establish the first basis 
  !
  nstep = 0
  call new_basis(xbs,xc(0:ks),info)
  nbas = nbas + 1
  ndff = ndff + 1
  if(info /= 0) then
    call log_msg('daeslv',info,'basis evaluation failed')
    call del_storage()
    return
  endif
  !
  tloc    = zer                      ! set local variables
  y(0:ko) = zer                      ! to zero
  !
  info = 1
  call locode(dnam,xbs,y(0:ko),xc(0:kx),zp(0:ko+nw),info)  
                                     ! evaluate local ode
  ndff = ndff + 1
  if(info /= 0) then
    call log_msg('daeslv',info,'local ode evaluation failed')
    call del_storage()
    return
  endif
  !
  tc    = xc(0)
  tnext = tc
  ptyp  = 'start'
  call outp(xc,ptyp,info)         ! print first point
  if(info == 1) then
    call del_storage()            ! return requested by user
    return
  endif
  !
  last = .false.                  ! last point indicator
  nstep = 0
  hold = h                        ! save the current step 
  kfail = 0
  !
  !----------------------------------------------------
  ! do loop for steps
  !----------------------------------------------------
  !
  do while (.not. last)
    !
    tprev = tc
    !
    ! check if we are close to the terminal value of t
    !
    if((tc + tfact*h - tout)*posneg > zer) then
      h = tout - tc 
      last = .true.
    endif
    !
    task = 'step'
    !
    ! loop point for calls to the ode solver
    !
    do while (task == 'step  ' .or. task == 'eval  ' &
                               .or. task == 'reduce')
      !
      if(interp) then
        call dopstn(task,tloc,y(0:ko),zp(0:ko),yint)
      else
        call dopstn(task,tloc,y(0:ko),zp(0:ko))
      endif
      !
      if(task /= 'eval') exit
      ! 
      ! compute global point
      !
      call eval_point(xbs,y(0:km),xn(0:ks), &
                        itstp,alpha,info)
      nff = nff + itstp
      if(info /= 0) then
        if(abs(h) < hfact*abs(hold) .or. kfail >= kfailmx) then
          kbas = 1
          task = 'newbas'
        else
          call log_msg('daeslv',info,'global point evaluation &
                                         &failed - reduce step')
          task = 'reduce'
          kfail = kfail + 1
        endif
        cycle
      endif
      !
      info = 0
      call locode(dnam,xbs,y(0:ko),xn(0:kx),zp(0:ko+nw),info)        
                                            ! compute local ode
      ndff = ndff + 1
      if(info /= 0) then
        call log_msg('daeslv',info,'local ode evaluation failed')
        call del_storage()
        return
      endif
      !
    enddo
    !
    if(task == 'done') then
      !
      ! output of new global point
      !
      tnext = xn(0)
      ptyp  = 'reg'
      call outp(xn(0:kx),ptyp,info) 
      if(info == 1) then
        call del_storage()
        return
      endif
      !
      xc(0:kx) = xn(0:kx)
      tc       = xc(0)
      if(last) then
        if(tc < tout) then
          last = .false.
        else
          exit
        endif
      endif
      !
      ! Check for new basis computation
      !
      kbas = 0
      if(abs(h) < hfact*abs(hold)) then    
        kbas = 1
        write(loglun,'(A,e16.8,A,e16.8)')'Basis check: &
                                     &h=',h,' hold=',hold
      endif     
      if(itstp > itstpmx) then
        kbas = 1
        write(loglun,'(A,I4)')'Basis check: itstp=',itstp
      endif      
      if(alpha > alphopt) then
        kbas = 1
        write(loglun,'(A,e16.8)')'Basis check: alpha=',alpha
      endif
    else
      if(task /= 'newbas') then
        !
        select case(task) 
          case ('stpcnt')
            info = -10
            call log_msg('daeslv',info,'step count exceeded')
          case ('minstp')
            info = -11
            call log_msg('daeslv',info,'step fell below hmin')
          case ('fail')
            info = -12
            call log_msg('daeslv',info,'ode routine failed')
          case default
            info = -12
            call log_msg('daeslv',info,'ode routine: unknown &
                                           &task-return')
        end select
        call del_storage()
        return
        !
      endif
    endif
    if(kbas /= 0) then
      write(loglun,'(A,I6)')'  compute new basis at nstep= ',nstep
      call reval_basis(xbs,xc(0:ks),info)
      nbas = nbas + 1
      ndff = ndff + 1
      if(info /= 0) then
        call log_msg('daeslv',info,'basis evaluation failed')
        call del_storage()
        return
      endif
      !
      tloc    = zer
      y(0:ko) = zer
      !h = (h+hold)*0.5_RP
      h = hold
      kfail = 0
      info = 1
      call locode(dnam,xbs,y(0:ko),xc(0:kx),zp(0:ko+nw),info)
      ndff = ndff + 1
      if(info /= 0) then
        call log_msg('daeslv',info,'local ode evaluation failed')
        call del_storage()
        return
      endif
    endif
    hold = h
    !
  enddo
  !
  if(abs(tprnt-tnext) > epmach) then
    ptyp = 'last'
    call outp(xc(0:kx),ptyp,info)
  endif
  !
  x(0:kx) = xc(0:kx)
  !
  info = 0
  call del_storage()
  return
  !
contains
  !  
  subroutine loc_storage()
    !
    ! allocate storage
    !
    call loc_basis(xbs,ks,km)
    call loc_work(kx,ko,info)
    !
    ! associate pointers with work arrays
    !
    xc   => wkxc
    xn   => wkxn
    xi   => wkxi
    y    => wky
    zp   => wkzp
    zi   => wkzi
    yint => wkint
    !
  end subroutine loc_storage
  !
  subroutine del_storage()
    !
    nullify(xc,xn,xi,y,yint)
    call del_basis(xbs)
    call del_work
    !
  end subroutine del_storage
  !
  subroutine outp(x,ptyp,info)
    ! 
    real(kind=RP), dimension(0:), intent(in)    :: x
    character(LEN=*),             intent(inout) :: ptyp
    integer,                      intent(out)   :: info
    !
    interface
      subroutine solout(nstp,xc,hinch,info)
        use def_kinds, ONLY : RP => RPX
        integer,                      intent(in) :: nstp
        real(kind=RP), dimension(0:), intent(in) :: xc
        real(kind=RP), intent(out)               :: hinch
        integer, intent(out)                     :: info
      end subroutine solout
    end interface
    !
    integer       :: itstp
    real(kind=RP) :: alpha,hstp,hinch,s,tinxt
    !
    !-----------------------------------------------------
    !
    hstp = tnext - tprev
    !
    ! check when only the current point is to be printed
    !
    if(ptyp == 'start' .or. ptyp == 'last' .or. &
           hstp == zer .or. (.not. interp)) then
      tprnt = tnext
      hinch = hint
      call solout(naccpt,x,hinch,info)
      if(info == 1) return
      if(hinch /= hint) then
        hint = hinch
        if(hint /= zer) interp = .true. 
      endif
      return
    endif 
    !
    ! Loop for interpolation in current interval
    !
    tinxt = tprnt + hint
    s = (tinxt - tprev)/hstp
    !
    do while (s > zer .and. s < one)
      !
      zi(0:ko) = yint(0,0:ko) &
            + hstp*s*(yint(1,0:ko) + s*(yint(2,0:ko) &
            + s*(yint(3,0:ko) + s*yint(4,0:ko))))
      call eval_point(xbs,zi,xi,itstp,alpha,info)
      nff = nff + itstp
      if(info /= 0) then
        call log_msg('outp',info,'global point &
           &interpolation failed -- skip output')
        return
      endif
      !
      info = 0
      if(ks < kx) call locode(dnam,xbs,zi(0:km),&
                           xi(0:kx),zi(0:ko),info)
      ptyp = 'int'
      tprnt = xi(0)
      hinch = hint
      call solout(naccpt,xi,hinch,info)
      if(info == 1) return
      if(hinch /= hint) then
        hint = hinch
        if(hint /= zer) interp = .true. 
      endif
      tinxt = tprnt + hint
      s = (tinxt - tprev)/hstp
      !
    enddo
    if(s < zer .or. s > one .or. tprnt == tnext) return  ! no printout
    if(s == one) then
      tprnt = tnext           ! print only the right endpoint
      hinch = hint
      call solout(naccpt,x,hinch,info)
      if(info == 1) return
      if(hinch /= hint) then
        hint = hinch
        if(hint /= zer) interp = .true. 
      endif
      return
    endif
    !
  end subroutine outp
  !
  subroutine dopstn(task,t,y,yp,yint)
    !
    implicit none
    !
    intrinsic abs, sqrt, max, min
    !
    character(LEN=6),                intent(inout) :: task
    real(kind=RP),                   intent(inout) :: t
    real(kind=RP), dimension(0:),    intent(inout) :: y
    real(kind=RP), dimension(0:),    intent(in)    :: yp
    real(kind=RP), dimension(0:,0:), intent(inout) :: yint
    optional :: yint
    !                   
    real(kind=RP), dimension(:), pointer :: w0,w1,w2,w3,w4,w5,w6
    !
    ! Dormand-Prince RK-5 coefficients for non-autonomous ODEs
    ! with interpolation
    !
    real(kind=RP), parameter :: a21 = 0.2_RP, a31 = 3.0_RP/40.0_RP, &
                       a32 = 9.0_RP/40.0_RP, a41 = 44.0_RP/45.0_RP, &
                      a42 = -56.0_RP/15.0_RP, a43 = 32.0_RP/9.0_RP, &
           a51 = 19372.0_RP/6561.0_RP, a52 = -25360.0_RP/2187.0_RP, &
              a53 = 64448.0_RP/6561.0_RP, a54 = -212.0_RP/729.0_RP, &
                a61 = 9017.0_RP/3168.0_RP, a62 = -355.0_RP/33.0_RP, &
                a63 = 46732.0_RP/5247.0_RP, a64 = 49.0_RP/176.0_RP, &
               a65 = -5103.0_RP/18656.0_RP, a71 = 35.0_RP/384.0_RP, &
                 a73 = 500.0_RP/1113.0_RP, a74 = 125.0_RP/192.0_RP, &
                 a75 = -2187.0_RP/6784.0_RP, a76 = 11.0_RP/84.0_RP, & 
         c2 = 0.2_RP, c3 = 0.3_RP, c4 = 0.8_RP, c5 = 8.0_RP/9.0_RP, &
               e1  = 71.0_RP/57600.0_RP, e3  = -71.0_RP/16695.0_RP, &
            e4  = 71.0_RP/1920.0_RP, e5  = -17253.0_RP/339200.0_RP, &
                     e6  = 22.0_RP/525.0_RP, e7  = -1.0_RP/40.0_RP, &
                  d21=-1337.0_RP/480.0_RP, d23=4216.0_RP/1113.0_RP, &
                   d24=-135.0_RP/80.0_RP, d25=-2187.0_RP/8480.0_RP, &
                       d26=66.0_RP/70.0_RP, d31=1039.0_RP/360.0_RP, &
                    d33=-468200.0_RP/83475.0_RP, d34=9.0_RP/2.0_RP, &
               d35=400950.0_RP/318000.0_RP, d36=-638.0_RP/210.0_RP, &
               d41=-1163.0_RP/1152.0_RP, d43=37900.0_RP/16695.0_RP, &
              d44=-415.0_RP/192.0_RP, d45=-674325.0_RP/508800.0_RP, &
                                             d46=374.0_RP/168.0_RP
    ! 
    ! Step parameters
    !
    real(kind=RP), parameter :: beta=0.04_RP, fac1=0.2_RP, &
                  fac2=10.0_RP, facmin=1.0E-4_RP, safe=0.9_RP
    real(kind=RP), parameter :: expo1=0.2_RP - beta*0.75_RP, &
                  facin1 = 1.0_RP/fac1, facin2 = 1.0_RP/fac2
    !
    real(kind=RP), parameter :: zer=0.0_RP,one=1.0_RP
    real(kind=RP), save :: t0
    real(kind=RP), save :: facold=1.0E-4_RP
    logical,       save :: reject
    !
    ! other local variables
    !
    integer             :: i
    integer,       save :: icall,my
    logical,       save :: interp
    real(kind=RP), save :: err,errfac,fac,hnew,sk,tmp
    !
    !-------------------------------------------------------------
    ! Subroutine for taking a step with the Dormand-Prince
    ! Runge Kutta method of order 5 in solving the
    ! NONAUTONOMOUS initial value problem
    !
    !        y' = f(t,y),  y(t0) = y0,
    !
    !        t    time
    !        y    vector of dimension N,
    !        f    mapping from R^N into R^N
    !
    ! The routine uses reverse communication and returns for all
    ! evaluations of the user function f. This is done under the
    ! control of the character variable task.
    !
    ! To take an RK step from a point (t, y), call the routine 
    ! with task = 'step'
    !
    ! The routine returns with task = 'eval' signifying a request
    ! for the evaluation of yp = f(t,y).
    ! 
    ! Return to the routine with the computed vector yp and without  
    ! any changes of task or any other variable.
    !
    ! Alternately, without computing yp, a forced step reduction
    ! may be requested by returning to the routine with 
    ! task = 'reduce' which prompts a reduction of the stepsize
    ! by a factor of 1/2. This feature is included to allow for the 
    ! solution of systems defined in terms of local coordinates on 
    ! a manifold.
    !
    ! Upon final return task may have either one of the following 
    ! values:
    !    task = 'minstp'  step fell below minimum steplength hmin
    !    task = 'stpcnt'  maximal number of steps nmax exceeded
    !    task = 'fail'   error condition.
    !------------------------------------------------------------------
    !
    my = size(y)-1                      ! index defining the dimension
    !
    interp = .false.                    ! set the interpolation
    if(present(yint)) interp = .true.   ! indicator
    !
    w0 => wk0                           ! associate
    w1 => wk1                           ! pointers
    w2 => wk2                           ! with
    w3 => wk3                           ! work
    w4 => wk4                           ! arrays
    w5 => wk5
    w6 => wk6
    !
    ! case selections on task
    !
    if(task /= 'step' .and. task /= 'reduce' &
                      .and. task /= 'eval') then
      !
      write(6,*)'wrong task for ode'
      task = 'fail'                     ! input error
      return
    endif
    !
    select case(task)
      !
      case('step')
        !
        reject = .false.
        t0 = t                        ! store initial t in t0
        w0(0:my) = y(0:my)            ! store initial y in w0
        w1(0:my) = yp(0:my)           ! store initial yp in w1
        icall = 0
      !
      case('reduce')
        !
        h = h*0.5_RP                  ! reduce step
        reject = .true.
        if(naccpt .ge. 1) nrejct = nrejct + 1
        icall = 0
      !
    end select
    !
    ! Stage evaluations
    !
    do                         ! loop until step is accepted
      !
      select case(icall)
        !
        case(0)
          !
          if(nstep >= nmax) then
            task = 'stpcnt'         !  maximal number of steps
            return
          endif
          !
          if(abs(h) < hmin) then
            task = 'minstp'         !  minimal stepsize
            return
          endif
          !
          ! First stage evaluation
          !
          y(0:my)  = w0(0:my) + h*a21*w1(0:my)
          t = t0 + c2*h 
          nstep = nstep + 1
          icall = 1
          task = 'eval'
          return
          !
        case(1)                     ! evaluate second stage
          !
          w2(0:my) = yp(0:my) 
          y(0:my)  = w0(0:my) + h*(a31*w1(0:my)+a32*w2(0:my))
          t = t0 + c3*h
          icall = 2
          return
          !
        case(2)                     ! evaluate third stage
          !    
          w3(0:my) = yp(0:my) 
          y(0:my)  = w0(0:my) + h*(a41*w1(0:my)+a42*w2(0:my) &
                                               +a43*w3(0:my))
          t = t0 + c4*h
          icall = 3
          return
          !
        case(3)                     ! evaluate fourth stage
          !
          w4(0:my) = yp(0:my) 
          y(0:my) = w0(0:my) + h*(a51*w1(0:my)+a52*w2(0:my) &
                             + a53*w3(0:my)+a54*w4(0:my))
          t = t0 + c5*h
          icall = 4
          return
          !
        case(4)                     ! evaluate fifth stage
          !
          w5(0:my) = yp(0:my) 
          y(0:my) = w0(0:my) + h*(a61*w1(0:my)+a62*w2(0:my) &
                    + a63*w3(0:my)+a64*w4(0:my)+a65*w5(0:my))
          t = t0 + h
          icall = 5
          return
          !
        case(5)                     ! evaluate last stage
          !
          w6(0:my) = yp(0:my) 
          y(0:my) = w0(0:my) + h*(a71*w1(0:my)+a73*w3(0:my) &
                     + a74*w4(0:my)+a75*w5(0:my)+a76*w6(0:my))
          !
          ! if continuous output is requested calculate uint 
          !
          if(interp) then 
            yint(0,0:my) = w0(0:my)
            yint(1,0:my) = w1(0:my)
            yint(2,0:my) = d21*w1(0:my) + d23*w3(0:my) + d24*w4(0:my) &
                                        + d25*w5(0:my) + d26*w6(0:my)
            yint(3,0:my) = d31*w1(0:my) + d33*w3(0:my) + d34*w4(0:my) &
                                        + d35*w5(0:my) + d36*w6(0:my)
            yint(4,0:my) = d41*w1(0:my) + d43*w3(0:my) + d44*w4(0:my) &
                                        + d45*w5(0:my) + d46*w6(0:my)
          endif
          icall = 6
          return
          !
        case(6)
          ! 
          w4(0:my) = h*(e1*w1(0:my)+e3*w3(0:my)+e4*w4(0:my) &
                        +e5*w5(0:my)+e6*w6(0:my)+e7*yp(0:my))
          ! 
          ! error estimation
          !
          err = zer
          do i=0,my 
            sk  = atol + rtol*max(abs(w0(i)),abs(y(i)))
            tmp = w4(i)/sk
            err = err + tmp*tmp
          enddo
          err = sqrt(err/real(my+1,kind=RP))
          !
          ! estimation of the next step including lund
          ! stabilization to ensure that fac1 <= hnew/h <= fac2 
          !
          errfac = err**expo1
          fac    = errfac/(facold**beta)
          fac    = max(facin2, min(facin1,fac/safe))
          hnew   = h/fac
          !
          ! error test
          !
          if(err .le. one) then
            !
            ! accepted step
            !
            facold = max(err,facmin)
            naccpt = naccpt+1
            !
            ! set the next stepsize
            !
            if(abs(hnew) .gt. abs(hmax) &
                    .and. abs(hmax) .gt. zer)hnew = posneg*hmax  
            if(reject)hnew = posneg*min(abs(hnew),abs(h))
            h = hnew
            !
            ! successful return
            !
            task = 'done'
            return
            !
          else
            !
            ! rejected step. try smaller h
            !
            hnew   = h/min(facin1,errfac/safe)
            reject = .true.
            if(naccpt .ge. 1) nrejct = nrejct+1 
            h = hnew
            icall = 0
            cycle
          endif
          !  
        case default
          !
          ! icall has an improper value
          !
          write(6,*)'wrong icall in ode'
          task = 'fail'
          return
          !  
      end select
      !
    enddo
    !
  end subroutine dopstn
  !
end subroutine daeslv