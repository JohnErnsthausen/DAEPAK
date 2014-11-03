!-------------------------------------------------------------------
!  Driver for applying the DAE-solver with daeq1 to solve the
!  daeq1 model example 
!
!       u1' = u2
!        0  = u1 - t*t + 2*u2
!
!  subject to the initial conditions
!
!       u1(1) = 11, u2(1) = -5
!
!  The exact solution is
!
!       u1(t) = 6 exp[(1-t)/2] + t^2 - 4*t + 8
!       u2(t) = -3 exp[(1-t)/2] + 2*t - 4
!------------------------------------------------------------------
!
module prob_dat
  !
  ! problem data for the model problem
  !
  use def_kinds, ONLY : RP => RPX
  implicit none
  !   
  real(kind=RP), parameter :: zer   = 0.0_RP
  real(kind=RP), parameter :: half  = 0.5_RP
  real(kind=RP), parameter :: one   = 1.0_RP
  real(kind=RP), parameter :: two   = 2.0_RP
  real(kind=RP), parameter :: three = 3.0_RP
  real(kind=RP), parameter :: four  = 4.0_RP
  real(kind=RP), parameter :: six   = 6.0_RP
  real(kind=RP), parameter :: eight = 8.0_RP
  !
end module prob_dat

program model_q1
  !
  use def_kinds, ONLY : RP => RPX
  use dae_dat
  use prob_dat
  implicit none
  intrinsic abs
  !
  interface
    subroutine daeslv(info)
      use def_kinds, ONLY : RP => RPX
      implicit none
      integer, intent(out) :: info
    end subroutine daeslv
  end interface
  !
  integer          :: info                  ! return indicator
  character(len=5) :: dnm     = 'daeq1'     ! dae type
  character(len=6) :: pname   = 'model'     ! problem name
  integer          :: outunit = 10          ! output unit
  integer          :: nvd     = 1, &        ! no. of diff. equ.               
                      nva     = 1, &        ! no. of alg. equ.
                      nvu     = 2, &        ! dim. of u
                      nvw     = 0           ! dim of w
  real(kind=RP)    :: t
  !
  !-----------------------------------------------------------
  !
  call setup(dnm,pname,outunit,nvd,nva,nvu,nvw,info)
  if(info /= 0) then
    if(info == -2) write(6,*)'unknown type of DAE'
    if(info == -1) write(6,*)'dimension error'
    write(6,*)'Stop due to input error'
    stop
  endif
  !
  ! set initial point, kx is set in setup
  !
  x(0:kx) = 0.0_RP
  t       = one
  x(0)    = t
  x(1)    = 11.0_RP                      ! u(1)
  x(2)    =-5.0_RP                       ! u(2)
  !
  ! set global data
  !
  reltol = 1.0e-8_RP                     ! relative tolerance
  abstol = reltol*1.0e-2_RP              ! absolute tolerance
  h    = 5.0e-3_RP                       ! first step
  hint = 0.1_RP                          ! interpolation step
  tout = 3.0_RP                          ! terminal time
  hmax = abs(tout-1.0_RP)                ! maximum step
  !
  ! print any desired information about the run, such as
  !
  write(outunit,20)reltol,abstol
  20 format(' tolerances: reltol=',e14.6,' abstol= ',e14.6)
  write(outunit,30)h,tout
  30 format(' h= ',e14.6,' tout= ',e14.6)
  !
  call daeslv(info)                      ! call the solver
  !
  call close_run(info)                   ! terminate run
  !
  write(6,*)'Stop with info= ',info      ! write stop message
  !
end program model_q1
!
subroutine daefct(fname,xc,vc,vec,mat,iret)
  !
  use def_kinds, ONLY : RP => RPX
  use dae_dat
  use prob_dat
  implicit none
  intrinsic exp
  !
  character(len=4),                intent(in)  :: fname
  real(kind=RP), dimension(0:),    intent(in)  :: xc
  real(kind=RP), dimension(0:),    intent(in)  :: vc
  real(kind=RP), dimension(0:),    intent(out) :: vec
  real(kind=RP), dimension(0:,0:), intent(out) :: mat
  integer,                         intent(out) :: iret
  optional vc,vec,mat
  !------------------------------------------------------
  ! Routine for evaluating the parts of the dae under
  ! control of the identifier fname.
  ! Return indicator:
  ! iret = 0   no error
  ! iret = 1   not available part as requested by fname
  ! iret = 2   wrong dimension for output
  !------------------------------------------------------
  !
  real(kind=RP) :: t
  !
  iret = 0
  !
  ! The vector vc is not used in DAEQ1.
  !
  if(present(vc)) then
    iret = 2
    return
  endif
  select case(fname)
    !
    case ('amt ')
      if(.not.present(mat)) then
        iret = 2
        return
      endif
      mat      = zer
      mat(0,0) = one
      !
    case ('hf  ')
      if(.not.present(vec)) then
        iret = 2
        return
      endif
      vec(0) = xc(2)
      !
    case ('ff  ')
      if(.not.present(vec)) then
        iret = 2
        return
      endif
      t = xc(0)
      vec(0) = xc(1) - t*t + two*xc(2)
      !
    case ('dff ')
      if(.not.present(mat)) then
        iret = 2
        return
      endif
      t = xc(0)
      mat(0,0)  =-two*t
      mat(0,1)  = one
      mat(0,2)  = two
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
  use lin_lib
  implicit none
  intrinsic exp, abs
  integer,                      intent(in) :: nstp
  real(kind=RP), dimension(0:), intent(in) :: xc
  real(kind=RP), intent(out)               :: hinch
  integer, intent(out)                     :: info
  !-----------------------------------------------------------
  ! Routine for intermediate printout of the solution vector
  ! x = (t,u(0:ku),up(0:ku),w(0:kw)).
  !
  ! nstp   current count of accepted steps
  ! hinch  a new value of the interpolation step or
  !        zero if no change is desired
  ! info   nonzero if integration is to stop or
  !        zero otherwise
  !-----------------------------------------------------------
  ! 
  integer :: i
  real(kind=RP) :: t,tmp
  real(kind=RP), dimension(0:2) :: tvec
  !
  t   = xc(0)
  tmp =  half*(one-t)
  tvec(0) = t
  tvec(1) =  six*exp(tmp) + ((t-four)*t + eight)
  tvec(2) = -three*exp(tmp) + two*t - four
  !
  write(outlun,10) nstp,t
  10 format(1X/' nstp = ',I6,' t= ',e18.11)
  !
  write(outlun,50)(xc(i), i=0,2,1)
  50 format(' u = ',3e16.8)
  write(outlun,60)(tvec(i),i=0,2,1)
  60 format(' ue= ',3e16.8)
  !
  write(outlun,70) ddist(xc(0:2),tvec(0:2))
  70 format(' ||u-ue||_2 = ',e24.16)
  !
  hinch = hint
  info = 0
  !
end subroutine solout

