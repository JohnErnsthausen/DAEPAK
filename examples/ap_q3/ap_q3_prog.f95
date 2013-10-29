!-------------------------------------------------------------
! Driver for applying the DAE solver with daeq3 for solving
! the second order problem 
!
! u1" + u1*w = exp[t]*(1 + sin[t])
! u2" + u2*w = 1/(1+t)*(2/(1+t)^2 + sin[t])
! 0  = exp[t]/(1+t)^3 - 3*exp[t]/(1+t)^2 + exp[1/(1+t)]/(1+t)
!                     -(u1*u2^2 - 3*u1*u2 + exp[u2])*u2 
!
! for which
!
!   u1 = exp[t], u2 = (1+t)^(-1), w = sin[t]
!
! is seen to be a solution.
!
! The DAE is a modified version of a problem given in
!
! U.M. Ascher, and L.R. Petzold, 
! Projected Implicit Runge-Kutta Methods for Differential 
! Algebraic Equations, 
! Lawrence Livermore National Lab., Num. Math. Group, 
! Technical Report UCRL-JC-104037
!
! The modification was constructed so that the existence 
! condition
!
!             ( A(t,u,p)     B(t,u,p) )
!        rank (                       )  = ku + kw 
!             ( D_uF(t,u)     0       )
!
! is violated at u1 = exp[s], u2 = (1+s)^(-1) with
! s .= 0.0375627.. Thus, when started at t = 0 the process
! fails near t = s. When started beyond the singular point
! there are no difficulties. Here only a run for the
! initial conditions
!
!  u1(0.5) =exp[0.5], u2(0.5) =2/3, 
!  u1'(0.5)=exp[0.5], u2'(0.5)=4/9
!
! is shown. For specified values of t, the output files 
! exhibits computed vectors x = (u1,u2,u1',u2',w) and
! the corresponding exact solution values
!    
!--------------------------------------------------------------
!
module prob_dat
  !
  ! problem data for the index-three Ascher-Petzold problem
  !
  use def_kinds, ONLY : RP => RPX
  implicit none
  !   
end module prob_dat

program ap_q3
  !
  use def_kinds, ONLY : RP => RPX
  use dae_dat
  use prob_dat
  implicit none
  intrinsic exp
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
  character(len=5) :: dnm     = 'daeq3'  ! dae type
  character(len=6) :: pname   = 'ap_q3'  ! problem name
  integer          :: outunit = 10       ! output unit
  integer          :: nvd     = 2, &     ! no. of diff. equ.               
                      nva     = 1, &     ! no. of alg. equ.
                      nvu     = 2, &     ! dim. of u
                      nvw     = 1        ! dim of w
  !
  real(kind=RP) :: t0
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
  ! initial point
  !
  x(0:kx) = 0.0_RP
  !
  t0   = 0.5_RP                        ! initial
  x(0) = t0                            ! time
  x(1) = exp(t0)                       ! u(1)
  x(2) = 1.0_RP/(1.0_RP + t0)          ! u(2)
  x(3) = x(1)                          ! u'(1)
  x(4) = -(x(2)*x(2))                  ! u'(2)
  !
  ! global data
  !
  reltol = 5.0e-9_RP                     ! relative tolerance
  abstol = reltol*1.0e-1_RP              ! absolute tolerance
  h    = 0.01_RP                         ! first step
  hint = 0.25_RP                         ! interpolation step
  tout = 2.0_RP                          ! terminal time
  hmax = 0.05                            ! maximum step
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
end program ap_q3
!
subroutine daefct(fname,xc,vc,vec,mat,iret)
  !
  use def_kinds, ONLY : RP => RPX
  use dae_dat
  use prob_dat
  implicit none
  intrinsic exp, sin
  !
  character(len=4),                intent(in)  :: fname
  real(kind=RP), dimension(0:),    intent(in)  :: xc
  real(kind=RP), dimension(0:),    intent(in)  :: vc
  real(kind=RP), dimension(0:),    intent(out) :: vec
  real(kind=RP), dimension(0:,0:), intent(out) :: mat
  integer,                         intent(out) :: iret
  optional vc,vec,mat
  !
  real(kind=RP) :: t,ot,sit,u2,d2f00,d2f12,d2f22
  real(kind=RP), parameter :: zer=0.0_RP, one=1.0_RP, &
                              two=2.0_RP, three=3.0_RP
  !------------------------------------------------------
  ! Routine for evaluating the parts of the dae under
  ! control of the identifier fname.
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
      mat(0,0) = one
      mat(0,1) = zer
      mat(1,0) = zer
      mat(1,1) = one
      !
    case ('bmt ')
      if(.not.present(mat)) then
        iret = 2
        return
      endif
      mat(0,0) = xc(1)
      mat(1,0) = xc(2)
      !
    case ('hf  ')
      if(.not.present(vec)) then
        iret = 2
        return
      endif
      t   = xc(0)
      ot  = one/(one + t)
      sit = sin(t)
      vec(0) = exp(t)*(one + sit)
      vec(1) = (two*ot*ot + sit)*ot
      !
    case ('ff  ')
      if(.not.present(vec)) then
        iret = 2
        return
      endif
      t  = xc(0)
      ot  = one/(one + t)
      u2 = xc(2)
      vec(0) = ot*(exp(t)*ot*(ot - three) + exp(ot)) &
                 - (xc(1)*u2*(u2 - three) + exp(u2))*u2
      !
    case ('dff ')
      if(.not.present(mat)) then
        iret = 2
        return
      endif
      t = xc(0)
      ot  = one/(one + t)
      u2 = xc(2)
      mat(0,0) = ot*ot*(exp(t)*((7.0_RP-three*ot)*ot-three) &
                                           -exp(ot)*(ot+one))
      mat(0,1) = u2*u2*(three - u2)
      mat(0,2) = three*xc(1)*u2*(two-u2) - (one+u2)*exp(u2)
      !
    case ('d2ff')
      if(.not.present(vc) .or. .not.present(vec)) then
        iret = 2
        return
      endif
      t = xc(0)
      ot = one/(one + t)
      u2 = xc(2)
      d2f00 = ot*ot*(exp(t)*(((ot-two)*12.0_RP*ot &
                               + 13.0_RP)*ot-3.0_RP) &
                 + exp(ot)*ot*((ot + 4.0_RP)*ot + two))
      d2f12 = 3.0_RP*u2*(two-u2)
      d2f22 = 6.0_RP*xc(1)*(one-u2) - exp(u2)*(two+u2)
      !
      vec(0) = d2f00*vc(0)*vc(0) + two*d2f12*vc(1)*vc(2) &
                                 + d2f22*vc(2)*vc(2)
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
  intrinsic exp,sin
  !
  integer,                      intent(in)  :: nstp
  real(kind=RP), dimension(0:), intent(in)  :: xc
  real(kind=RP),                intent(out) :: hinch
  integer,                      intent(out) :: info
  !
  real(kind=RP) :: t,ue1,ue2,upe1,upe2,we
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
  t    = xc(0)
  ue1  = exp(t)
  ue2  = 1.0_RP/(1.0_RP + t)
  upe1 = ue1
  upe2 = -(ue2*ue2)
  we   = sin(t)
  !
  write(outlun,10) nstp,t
  10 format(1X/' nstp = ',I6,' t= ',e18.11)
  !
  write(outlun,50)xc(1:5)
  50 format(' x = '/(3e16.8))
  write(outlun,60)ue1,ue2,upe1,upe2,we
  60 format(' xe= '/(3e16.8))
  !
end subroutine solout
!
!-------------------------------------------------------------
!
! Output from this program:
!
!
! log file ap_q3.log
!
! open_log :: info=   0 opened  ap_q3.log
! open_file :: info=   0 opened  ap_q3.out
! Basis check: itstp=  25
!   compute new basis at nstep=      8
! Basis check: itstp=  24
!   compute new basis at nstep=     16
! Basis check: itstp=  25
!   compute new basis at nstep=     25
! Basis check: itstp=  21
!   compute new basis at nstep=     34
! main :: info=   0 Run completed
! close_file :: info=   0  closed  ap_q3.out
!
!
! output file ap_q3.out
!
! problem name: ap_q3 
! DAE type    : daeq3
! nv=  2 na=  1 nu=  2 nw=  1
! tolerances: reltol=  0.500000E-08 abstol=   0.500000E-09
! h=   0.100000E-01 tout=   0.200000E+01
!
! nstp =      0 t=  0.50000000000E+00
! x = 
!  0.16487213E+01  0.66666667E+00  0.16487213E+01 
! -0.44444444E+00  0.47942554E+00
! xe= 
!  0.16487213E+01  0.66666667E+00  0.16487213E+01 
! -0.44444444E+00  0.47942554E+00
!
! nstp =      9 t=  0.75000000000E+00
! x = 
!  0.21170000E+01  0.57142857E+00  0.21170000E+01 
! -0.32653061E+00  0.68163876E+00
! xe= 
!  0.21170000E+01  0.57142857E+00  0.21170000E+01 
! -0.32653061E+00  0.68163876E+00
!
! nstp =     15 t=  0.10000000000E+01
! x = 
!  0.27182818E+01  0.50000000E+00  0.27182818E+01 
! -0.25000000E+00  0.84147099E+00
! xe= 
!  0.27182818E+01  0.50000000E+00  0.27182818E+01 
! -0.25000000E+00  0.84147098E+00
!
! nstp =     21 t=  0.12500000000E+01
! x = 
!  0.34903429E+01  0.44444444E+00  0.34903429E+01 
! -0.19753086E+00  0.94898463E+00
! xe= 
!  0.34903430E+01  0.44444444E+00  0.34903430E+01 
! -0.19753086E+00  0.94898462E+00
!
! nstp =     27 t=  0.15000000000E+01
! x = 
!  0.44816891E+01  0.40000000E+00  0.44816891E+01 
! -0.16000000E+00  0.99749498E+00
! xe= 
!  0.44816891E+01  0.40000000E+00  0.44816891E+01 
! -0.16000000E+00  0.99749499E+00
!
! nstp =     32 t=  0.17500000000E+01
! x = 
!  0.57546027E+01  0.36363636E+00  0.57546027E+01 
! -0.13223141E+00  0.98398595E+00
! xe= 
!  0.57546027E+01  0.36363636E+00  0.57546027E+01 
! -0.13223140E+00  0.98398595E+00
!
! nstp =     38 t=  0.20000000000E+01
! x = 
!  0.73890561E+01  0.33333333E+00  0.73890561E+01 
! -0.11111111E+00  0.90929743E+00
! xe= 
!  0.73890561E+01  0.33333333E+00  0.73890561E+01 
! -0.11111111E+00  0.90929743E+00
! Run completed with info=  0
! Statistics
! total number of steps =  42
!        accepted steps =  38
!        rejected steps =  4
! local bases computed  =  5
! function evaluations  =  2803
! jacobian evaluations  =  262
!-------------------------------------------------------------
