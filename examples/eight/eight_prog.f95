!-----------------------------------------------------------------
! Driver for applying the DAE-solver with daen1 for solving
! the "figure 8" problem
!
!   u^2 + u'^2 - 1  = 0
!   2*u*u' - w      = 0 
!   u(0) = 0, u'(0) = 1, w(0) = 0
!
! with the exact solution
!
!   u(t) = sin t, w(t) = sin 2t
!
! Note that here the (nu + nw) x (nu + nw) matrix (D_pF, D_wF)
! has the form
!                 (2*p    0)
!                 (2*u   -1)
! and hence is nonsingular only as long as  p .ne. 0. This does
! not hold for t = k*pi/2. These points are removable 
! singularities of the solution curve which the code passes 
! easily but where dassl fails.
!-----------------------------------------------------------------
!
module prob_dat
  !
  ! no problem data for the "figure 8" problem
  !
end module prob_dat

program eight
  !
  use def_kinds, ONLY : RP => RPX
  use dae_dat
  use prob_dat
  implicit none
  !
  interface
    subroutine daeslv(info)
      use def_kinds, ONLY : RP => RPX
      implicit none
      integer, intent(out) :: info
    end subroutine daeslv
  end interface
  !
  integer          :: info               ! return indicator
  character(len=5) :: dnm     = 'daen1'  ! dae type
  character(len=6) :: pname   = 'eight'  ! problem name
  integer          :: outunit = 10       ! output unit
  integer          :: nvd     = 1, &     ! no. of diff. equ.               
                      nva     = 2, &     ! no. of alg. equ.
                      nvu     = 1, &     ! dim. of u
                      nvw     = 1        ! dim of w
  !
  !-----------------------------------------------------------
  !
  call setup(dnm,pname,outunit,nvd,nva,nvu,nvw,info)
  if(info /= 0) then
    if(info == -1) write(6,*)'unknown type of DAE'
    if(info == -2) write(6,*)'dimension error'
    write(6,*)'Stop due to input error'
    stop
  endif
  !
  x(0:kx) = 0.0_RP                     ! set initial point
  x(2)    = 1.0_RP
  !
  ! set global data
  !
  reltol = 1.0e-8_RP                   ! relative tolerance
  abstol = reltol*1.0e-1_RP            ! absolute tolerance
  h    = 0.001_RP                      ! first step
  hint = 0.5_RP                        ! interpolation step
  tout = 6.0_RP                        ! set terminal time
  hmax = 1.0_RP                        ! maximal step
  !
  ! print any desired information about the run, such as
  !
  write(outunit,20)reltol,abstol
  20 format(' tolerances: reltol=',e14.6,' abstol= ',e14.6)
  write(outunit,30)h,tout
  30 format(' h= ',e14.6,' tout= ',e14.6)
  !
  call daeslv(info)                     ! call the solver
  !
  call close_run(info)                  ! terminate run
  !
  write(6,*)'Stop with info= ',info     ! write stop message
  !
end program eight
!
subroutine daefct(fname,xc,vc,vec,mat,iret)
  !
  use def_kinds, ONLY : RP => RPX
  use dae_dat
  use prob_dat
  implicit none
  !
  character(len=4),                intent(in)  :: fname
  real(kind=RP), dimension(0:),    intent(in)  :: xc
  real(kind=RP), dimension(0:),    intent(in)  :: vc
  real(kind=RP), dimension(0:),    intent(out) :: vec
  real(kind=RP), dimension(0:,0:), intent(out) :: mat
  integer,                         intent(out) :: iret
  optional vc,vec,mat
  !
  !------------------------------------------------------
  ! Routine for evaluating the parts of the dae under
  ! control of the identifier fname.
  ! Return indicator:
  ! iret = 0   no error
  ! iret = 1   not available part as requested by fname
  ! iret = 2   wrong dimension for output
  !------------------------------------------------------
  !
  iret = 0 
  select case(fname)
    !
    case ('ff  ')
      if(.not.present(vec)) then
        iret = 2
        return
      endif
      vec(0) = xc(1)*xc(1) + xc(2)*xc(2) - 1.0_RP
      vec(1) = 2.0_RP*xc(1)*xc(2) - xc(3)
      !
    case ('dff ')
      if(.not.present(mat)) then
        iret = 2
        return
      endif
      mat(0,0) = 0.0_RP
      mat(0,1) = 2.0_RP*xc(1)
      mat(0,2) = 2.0_RP*xc(2)
      mat(0,3) = 0.0_RP
      mat(1,0) = 0.0_RP
      mat(1,1) = 2.0_RP*xc(2)
      mat(1,2) = 2.0_RP*xc(1)
      mat(1,3) = -1.0_RP
      !
    case default
      iret = 1
      !
  end select
  !
end subroutine daefct
!
subroutine solout(nstp,xc,hinch,info)
  !
  use def_kinds, ONLY : RP => RPX
  use dae_dat
  use prob_dat
  implicit none
  intrinsic sin, cos
  integer,                      intent(in)  :: nstp
  real(kind=RP), dimension(0:), intent(in)  :: xc
  real(kind=RP),                intent(out) :: hinch
  integer,                      intent(out) :: info
  !
  real(kind=RP) :: t, ue, upe, we 
  !-----------------------------------------------------------
  ! Routine for intermediate printout of the solution vector
  ! x = (t,u(0:ku),up(0:ku)).
  !
  ! nstp   current count of accepted steps
  ! hinch  a new value of the interpolation step or
  !        zero if no change is desired
  ! info   nonzero if integration is to stop or
  !        zero otherwise
  !-----------------------------------------------------------
  t = xc(0)
  ue  = sin(t)
  upe = cos(t)
  we  = sin(t+t)
  !
  write(outlun,10) nstp,t
  10 format(1X/' nstp = ',I6,' t= ',e18.11)
  !
  write(outlun,20)xc(1),xc(2),xc(3)
  20 format('  x = ',3e18.11)
  write(outlun,30)ue,upe,we
  30 format(' xe = ',3e18.11)
  !
end subroutine solout
!
!-------------------------------------------------------------
!
! Output from this program:
!
!
! log file eight.log
!
! open_log :: info=   0 opened  eight.log
! open_file :: info=   0 opened  eight.out
! Basis check: itstp=  22
!   compute new basis at nstep=     12
! Basis check: itstp=  22
!   compute new basis at nstep=     21
! Basis check: itstp=  22
!   compute new basis at nstep=     30
! Basis check: itstp=  23
!   compute new basis at nstep=     37
! Basis check: itstp=  23
!   compute new basis at nstep=     44
! Basis check: itstp=  23
!   compute new basis at nstep=     54
! Basis check: itstp=  22
!   compute new basis at nstep=     62
! Basis check: itstp=  23
!   compute new basis at nstep=     69
! Basis check: itstp=  30
! Basis check: alpha=  0.52826816E+00
!   compute new basis at nstep=     77
! Basis check: itstp=  23
!   compute new basis at nstep=     86
! Basis check: itstp=  22
!   compute new basis at nstep=     95
! Basis check: itstp=  24
!   compute new basis at nstep=    106
! Basis check: itstp=  21
!   compute new basis at nstep=    114
! Basis check: itstp=  22
!   compute new basis at nstep=    123
! main :: info=   0 Run completed
! close_file :: info=   0  closed  eight.out
!
!
! output file eight.out
!
! problem name: eight 
! DAE type    : daen1
! nv=  1 na=  2 nu=  1 nw=  1
! tolerances: reltol=  0.100000E-07 abstol=   0.100000E-08
! h=   0.100000E-02 tout=   0.600000E+01
!
! nstp =      0 t=  0.00000000000E+00
!  x =  0.00000000000E+00 0.10000000000E+01 0.00000000000E+00
! xe =  0.00000000000E+00 0.10000000000E+01 0.00000000000E+00
!
! nstp =     12 t=  0.50000000000E+00
!  x =  0.47942553502E+00 0.87758256789E+00 0.84147099215E+00
! xe =  0.47942553860E+00 0.87758256189E+00 0.84147098481E+00
!
! nstp =     25 t=  0.10000000000E+01
!  x =  0.84147098725E+00 0.54030230530E+00 0.90929743161E+00
! xe =  0.84147098481E+00 0.54030230587E+00 0.90929742683E+00
!
! nstp =     36 t=  0.15000000000E+01
!  x =  0.99749498777E+00 0.70737201590E-01 0.14112001066E+00
! xe =  0.99749498660E+00 0.70737201668E-01 0.14112000806E+00
!
! nstp =     43 t=  0.20000000000E+01
!  x =  0.90929743074E+00-0.41614683210E+00-0.75680249425E+00
! xe =  0.90929742683E+00-0.41614683655E+00-0.75680249531E+00
!
! nstp =     57 t=  0.25000000000E+01
!  x =  0.59847214705E+00-0.80114361421E+00-0.95892427941E+00
! xe =  0.59847214410E+00-0.80114361555E+00-0.95892427466E+00
!
! nstp =     65 t=  0.30000000000E+01
!  x =  0.14112000989E+00-0.98999249772E+00-0.27941550522E+00
! xe =  0.14112000806E+00-0.98999249660E+00-0.27941549820E+00
!
! nstp =     72 t=  0.35000000000E+01
!  x = -0.35078321565E+00-0.93645669457E+00 0.65698658697E+00
! xe = -0.35078322769E+00-0.93645668729E+00 0.65698659872E+00
!
! nstp =     85 t=  0.40000000000E+01
!  x = -0.75680249355E+00-0.65364362459E+00 0.98935825168E+00
! xe = -0.75680249531E+00-0.65364362086E+00 0.98935824662E+00
!
! nstp =     89 t=  0.45000000000E+01
!  x = -0.97753012206E+00-0.21079579200E+00 0.41211846880E+00
! xe = -0.97753011767E+00-0.21079579943E+00 0.41211848524E+00
!
! nstp =     96 t=  0.50000000000E+01
!  x = -0.95892427817E+00 0.28366217804E+00-0.54402110143E+00
! xe = -0.95892427466E+00 0.28366218546E+00-0.54402111089E+00
!
! nstp =    109 t=  0.55000000000E+01
!  x = -0.70554032450E+00 0.70866977879E+00-0.99999021653E+00
! xe = -0.70554032557E+00 0.70866977429E+00-0.99999020655E+00
!
! nstp =    114 t=  0.60000000000E+01
!  x = -0.27941551091E+00 0.96017029359E+00-0.53657292500E+00
! xe = -0.27941549820E+00 0.96017028665E+00-0.53657291800E+00
! Run completed with info=  0
! Statistics
! total number of steps =  128
!        accepted steps =  114
!        rejected steps =  14
! local bases computed  =  15
! function evaluations  =  8514
! jacobian evaluations  =  798
!----------------------------------------------------------------
