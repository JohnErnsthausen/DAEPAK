!-------------------------------------------------------------------
!  Driver for applying the DAE-solver with daeq2 to the 
!  Ascher-Petzold problem of index two
!
!  u1' - a*(2 - t) w = [a - 1/(2-t)]*u1 + [(3-t)/(2-t)]*exp(t)
!  u2' -   (a - 1) w = [(a - 1)*u1]/(2-t) - u2 + 2*exp(t)
!  0 = (2 + t)*u1 + (t^2 - 4)*u2 - (t^2 + t - 2)*exp(t)
!
!  with a = 50 and the initial conditions
!
!  u1(0) = 1, u2(0) = 1, w(t0) = -1/2
!
!  Reference: U. Ascher and L. Petzold
!             Lawrence Livermore Nat'l Lab.
!             Num. Math. Group, Tech. Rep. UCRL-JC-107441, 1991
!
!  The exact solution is
!
!  u1(t) = exp[t], u2(t) = exp[t], w = -exp[t]/(2-t)  
!------------------------------------------------------------------
!
module prob_dat
  !
  ! problem data for the index-two Ascher-Petzold problem
  !
  use def_kinds, ONLY : RP => RPX
  implicit none
  !   
  real(kind=RP), parameter :: ap = 50.0_RP
  !
end module prob_dat

program ap_q2
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
  integer          :: info               ! return indicator
  character(len=5) :: dnm     = 'daeq2'  ! dae type
  character(len=6) :: pname   = 'ap_q2'  ! problem name
  integer          :: outunit = 10       ! output unit
  integer          :: nvd     = 2, &     ! no. of diff. equ.               
                      nva     = 1, &     ! no. of alg. equ.
                      nvu     = 2, &     ! dim. of u
                      nvw     = 1        ! dim of w
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
  ! set initial point
  !
  x(0:kx) = 0.0_RP
  x(1)    = 1.0_RP                       ! u(1)
  x(2)    = 1.0_RP                       ! u(2)
  !
  ! set global data
  !
  reltol = 1.0e-8_RP                     ! relative tolerance
  abstol = reltol*1.0e-2_RP              ! absolute tolerance
  h    = 5.0e-3_RP                       ! first step
  hint = 0.1_RP                          ! interpolation step
  tout = 1.0_RP                          ! terminal time
  hmax = abs(tout-x(0))                  ! maximum step
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
end program ap_q2
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
  !
  real(kind=RP) :: t, et, twmt
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
    case ('amt ')
      if(.not.present(mat)) then
        iret = 2
        return
      endif
      mat      = 0.0_RP
      mat(0,0) = 1.0_RP
      mat(1,1) = 1.0_RP
      !
    case ('bmt ')
      if(.not.present(mat)) then
        iret = 2
        return
      endif
      mat(0,0) = -ap*(2.0_RP - xc(0))
      mat(1,0) = -ap + 1.0_RP
      !
    case ('hf  ')
      if(.not.present(vec)) then
        iret = 2
        return
      endif
      t      = xc(0)
      twmt   = 2.0_RP - t
      et     = exp(t)
      vec(0) = (ap - (1.0_RP/twmt))*xc(1) &
              + (et*(3.0_RP - t))/twmt
      vec(1) = ((ap - 1.0_RP)*xc(1))/twmt - xc(2) + 2.0_RP*et
      !
    case ('ff  ')
      if(.not.present(vec)) then
        iret = 2
        return
      endif
      t = xc(0)
      vec(0) = (2.0_RP + t)*xc(1) + (t*t - 4.0_RP)*xc(2) &
                        - (t*(t + 1.0_RP) - 2.0_RP)*exp(t)
      !
    case ('dff ')
      if(.not.present(mat)) then
        iret = 2
        return
      endif
      t = xc(0)
      mat(0,0)  = xc(1) + 2.0_RP*t*xc(2) &
                   - (t*(t + 3.0_RP) - 1.0_RP)*exp(t)
      mat(0,1)  = 2.0_RP + t
      mat(0,2)  = t*t - 4.0_RP
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
  intrinsic exp
  integer,                      intent(in) :: nstp
  real(kind=RP), dimension(0:), intent(in) :: xc
  real(kind=RP), intent(out)               :: hinch
  integer, intent(out)                     :: info
  ! 
  real(kind=RP) :: t,ue1,ue2,we
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
  t   = xc(0)
  ue1 = exp(t)
  ue2 = ue1
  we  = ue1/(t - 2.0_RP)
  !
  write(outlun,10) nstp,t
  10 format(1X/' nstp = ',I6,' t= ',e18.11)
  !
  write(outlun,50)xc(1),xc(2),xc(5)
  50 format(' u = ',2e16.8,' w = ',e16.8)
  write(outlun,60)xc(3),xc(4)
  60 format(' up= ',2e16.8)
  write(outlun,70)ue1,ue2,we
  70 format(' ue= ',2e16.8,' we= ',e16.8)
  !
end subroutine solout
!
!-------------------------------------------------------------
!
! Output from this program:
!
!
! log file  ap_q2.log
!
! open_log :: info=   0 opened  ap_q2.log
! open_file :: info=   0 opened  ap_q2.out
! main :: info=   0 Run completed
! close_file :: info=   0  closed  ap_q2.out
!
!
! output file  ap_q2.out
!
! problem name: ap_q2 
! DAE type    : daeq2
! nv=  2 na=  1 nu=  2 nw=  1
! tolerances: reltol=  0.100000E-07 abstol=   0.100000E-09
! h=   0.500000E-02 tout=   0.100000E+01
!
! nstp =      0 t=  0.00000000000E+00
! u =   0.10000000E+01  0.10000000E+01 w =  -0.50000000E+00
! up=   0.10000000E+01  0.10000000E+01
! ue=   0.10000000E+01  0.10000000E+01 we=  -0.50000000E+00
!
! nstp =     18 t=  0.10000000000E+00
! u =   0.11051709E+01  0.11051709E+01 w =  -0.58166890E+00
! up=   0.11051710E+01  0.11051710E+01
! ue=   0.11051709E+01  0.11051709E+01 we=  -0.58166890E+00
!
! nstp =     34 t=  0.20000000000E+00
! u =   0.12214028E+01  0.12214028E+01 w =  -0.67855709E+00
! up=   0.12214028E+01  0.12214028E+01
! ue=   0.12214028E+01  0.12214028E+01 we=  -0.67855709E+00
!
! nstp =     48 t=  0.30000000000E+00
! u =   0.13498588E+01  0.13498588E+01 w =  -0.79403460E+00
! up=   0.13498586E+01  0.13498587E+01
! ue=   0.13498588E+01  0.13498588E+01 we=  -0.79403459E+00
!
! nstp =     62 t=  0.40000000000E+00
! u =   0.14918247E+01  0.14918247E+01 w =  -0.93239044E+00
! up=   0.14918244E+01  0.14918245E+01
! ue=   0.14918247E+01  0.14918247E+01 we=  -0.93239044E+00
!
! nstp =     75 t=  0.50000000000E+00
! u =   0.16487213E+01  0.16487213E+01 w =  -0.10991475E+01
! up=   0.16487211E+01  0.16487212E+01
! ue=   0.16487213E+01  0.16487213E+01 we=  -0.10991475E+01
!
! nstp =     88 t=  0.60000000000E+00
! u =   0.18221188E+01  0.18221188E+01 w =  -0.13015134E+01
! up=   0.18221186E+01  0.18221187E+01
! ue=   0.18221188E+01  0.18221188E+01 we=  -0.13015134E+01
!
! nstp =    101 t=  0.70000000000E+00
! u =   0.20137527E+01  0.20137527E+01 w =  -0.15490406E+01
! up=   0.20137521E+01  0.20137522E+01
! ue=   0.20137527E+01  0.20137527E+01 we=  -0.15490405E+01
!
! nstp =    113 t=  0.80000000000E+00
! u =   0.22255409E+01  0.22255409E+01 w =  -0.18546174E+01
! up=   0.22255422E+01  0.22255420E+01
! ue=   0.22255409E+01  0.22255409E+01 we=  -0.18546174E+01
!
! nstp =    126 t=  0.90000000000E+00
! u =   0.24596031E+01  0.24596031E+01 w =  -0.22360029E+01
! up=   0.24596027E+01  0.24596027E+01
! ue=   0.24596031E+01  0.24596031E+01 we=  -0.22360028E+01
!
! nstp =    138 t=  0.10000000000E+01
! u =   0.27182818E+01  0.27182818E+01 w =  -0.27182818E+01
! up=   0.27182818E+01  0.27182818E+01
! ue=   0.27182818E+01  0.27182818E+01 we=  -0.27182818E+01
! Run completed with info=  0
! Statistics
! total number of steps =  138
!        accepted steps =  138
!        rejected steps =  0
! local bases computed  =  1
! function evaluations  =  4166
! jacobian evaluations  =  830
!-----------------------------------------------------------------