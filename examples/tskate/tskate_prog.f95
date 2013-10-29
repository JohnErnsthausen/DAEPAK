!----------------------------------------------------------
! Driver for applying the DAE-solver with nohol to the 
! model of a scating body moving on a tilted plane in 
! u, phi coordinates.
!
!  u1'' - sin alph - w sin t = 0
!  u2''            + w cos t = 0
!  -u1'*sin t + u2'*cos t    = 0
!
!  u1(0) = 0.0, u2(0) = 0.0, u3(0) = 0.0, u4(0) = 1.0  
!
! The exact solution is
!
!  u1(t) = (1/2)*(sin alph)*(sin t)^2   
!  u2(t) = (1/2)*(sin alph)*(t - (1/2)sin 2t)
!  w(t)  = -2*(sin alph)*sin t
!
!  The problem was given in
!
!     Ju. I. Neimark and N. A. Fufaev
!     Dynamics of nonholonomic systems
!     Trans;. of Mathem. Monographs Vol 33
!     Amer. Math. Soc., Providence, RI 1972
!     p 108
!----------------------------------------------------------------
!
module prob_dat
  !
  use def_kinds, ONLY : RP => RPX
  implicit none
  real(kind=RP), parameter :: alph = 0.8_RP
  real(kind=RP), save      :: sina  
  !   
end module prob_dat

program tskate
  !
  use def_kinds, ONLY : RP => RPX
  use dae_dat
  use prob_dat
  implicit none
  intrinsic sin
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
  character(len=5) :: dnm     = 'nohol'  ! dae type
  character(len=6) :: pname   = 'tskate' ! problem name
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
    if(info == -1) write(6,*)'unknown type of DAE'
    if(info == -2) write(6,*)'dimension error'
    write(6,*)'Stop due to input error'
    stop
  endif
  !
  ! set initial point
  !
  x(0:4) = 0.0_RP
  x(4) = 1.0_RP
  !
  ! set global data
  !
  sina = sin(alph)
  !
  reltol = 1.0e-8_RP                     ! relative tolerance
  abstol = reltol*1.0e-1_RP              ! absolute tolerance
  h    = 0.05_RP                         ! first step
  hint = 0.5_RP                          ! interpolation step
  tout = 10.0_RP                         ! terminal time
  hmax = tout-x(0)                       ! maximum step
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
end program tskate
!
subroutine daefct(fname,xc,vc,vec,mat,iret)
  !
  use def_kinds, ONLY : RP => RPX
  use dae_dat
  use prob_dat
  implicit none
  intrinsic sin, cos
  !
  character(len=4),                intent(in)  :: fname
  real(kind=RP), dimension(0:),    intent(in)  :: xc
  real(kind=RP), dimension(0:),    intent(in)  :: vc
  real(kind=RP), dimension(0:),    intent(out) :: vec
  real(kind=RP), dimension(0:,0:), intent(out) :: mat
  integer,                         intent(out) :: iret
  optional vc,vec,mat
  !
  real(kind=RP)            :: si3, co3
  real(kind=RP), parameter :: zer=0.0_RP, one=1.0_RP
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
      mat(0,0) = one
      mat(1,1) = one
      !
    case ('hf  ')
      if(.not.present(vec)) then
        iret = 2
        return
      endif
      vec(0) = sin(alph)
      vec(1) = zer
      !
    case ('ff  ')
      if(.not.present(vec)) then
        iret = 2
        return
      endif
      vec(0) = -xc(3)*sin(xc(0)) + xc(4)*cos(xc(0))
      !
    case ('dff ')
      if(.not.present(mat)) then
        iret = 2
        return
      endif
      si3 = sin(xc(0))
      co3 = cos(xc(0))
      !
      mat(0,0)  = -xc(3)*co3 - xc(4)*si3
      mat(0,1)  = zer
      mat(0,2)  = zer
      mat(0,3)  = -si3
      mat(0,4)  = co3
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
  integer, intent(out)                      :: info
  !
  real(kind=RP) :: ue1,ue2,upe1,upe2,we,t,tt,sit
  !----------------------------------------------------------
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
  t    = xc(0)
  tt   = t+t
  sit  = sin(t)
  ue1  = 0.5_RP*sina*sit*sit
  upe1 = sina*sit*cos(t)
  ue2 = 0.5_RP*sina*(t - 0.5_RP*sin(tt))
  upe2 = 0.5_RP*sina*(1.0_RP - cos(tt))
  we  = -2.0_RP*sina*sit
  !
  write(outlun,10) nstp,t
  10 format(1X/' nstp = ',I6,' t= ',e18.11)
  !
  write(outlun,20)xc(1:nu)
  20 format(' u  = ',2e16.8)
  write(outlun,25)ue1,ue2
  25 format(' ue = ',2e16.8)
  write(outlun,30)xc(nu+1:ks),xc(ks+nu+1:kx)
  30 format(' up = ',2e16.8,' w  = ',e16.8)
  write(outlun,35)upe1,upe2,we
  35 format(' upe= ',2e16.8,' we = ',e16.8)
  !
end subroutine solout
!
!-------------------------------------------------------------
!
! Output from this program:
!
!
! log file tskate.log
!
! open_log :: info=   0 opened  tskate.log
! open_file :: info=   0 opened  tskate.out
! Basis check: itstp=  23
!   compute new basis at nstep=      8
! Basis check: itstp=  24
!   compute new basis at nstep=     20
! Basis check: itstp=  22
!   compute new basis at nstep=     33
! Basis check: itstp=  22
!   compute new basis at nstep=     48
! Basis check: itstp=  22
!   compute new basis at nstep=     61
! Basis check: itstp=  24
!   compute new basis at nstep=     74
! Basis check: itstp=  23
!   compute new basis at nstep=     88
! Basis check: itstp=  23
!   compute new basis at nstep=    101
! Basis check: itstp=  22
!   compute new basis at nstep=    114
! Basis check: itstp=  21
!   compute new basis at nstep=    129
! main :: info=   0 Run completed
! close_file :: info=   0  closed  tskate.out
!
!
! output file tskate.out
!
! problem name: tskate
! DAE type    : nohol
! nv=  2 na=  1 nu=  2 nw=  1
! tolerances: reltol=  0.100000E-07 abstol=   0.100000E-08
! h=   0.500000E-01 tout=   0.100000E+02
!
! nstp =      0 t=  0.00000000000E+00
! u  =   0.00000000E+00  0.00000000E+00
! ue =   0.00000000E+00  0.00000000E+00
! up =   0.00000000E+00  0.00000000E+00 w  =   0.00000000E+00
! upe=   0.00000000E+00  0.00000000E+00 we =   0.00000000E+00
!
! nstp =      7 t=  0.50000000000E+00
! u  =   0.82441736E-01  0.28430437E-01
! ue =   0.82441735E-01  0.28430439E-01
! up =   0.30181717E+00  0.16488347E+00 w  =  -0.68783766E+00
! upe=   0.30181717E+00  0.16488347E+00 we =  -0.68783766E+00
!
! nstp =     13 t=  0.10000000000E+01
! u  =   0.25397039E+00  0.19560553E+00
! ue =   0.25397039E+00  0.19560553E+00
! up =   0.32614503E+00  0.50794078E+00 w  =  -0.12072687E+01
! upe=   0.32614502E+00  0.50794078E+00 we =  -0.12072687E+01
!
! nstp =     20 t=  0.15000000000E+01
! u  =   0.35688331E+00  0.51270875E+00
! ue =   0.35688331E+00  0.51270874E+00
! up =   0.50616651E-01  0.71376661E+00 w  =  -0.14311182E+01
! upe=   0.50616649E-01  0.71376662E+00 we =  -0.14311182E+01
!
! nstp =     26 t=  0.20000000000E+01
! u  =   0.29656284E+00  0.85308031E+00
! ue =   0.29656283E+00  0.85308031E+00
! up =  -0.27144844E+00  0.59312567E+00 w  =  -0.13045801E+01
! upe=  -0.27144844E+00  0.59312566E+00 we =  -0.13045801E+01
!
! nstp =     32 t=  0.25000000000E+01
! u  =   0.12846733E+00  0.10686677E+01
! ue =   0.12846732E+00  0.10686677E+01
! up =  -0.34394508E+00  0.25693465E+00 w  =  -0.85863527E+00
! upe=  -0.34394508E+00  0.25693465E+00 we =  -0.85863528E+00
!
! nstp =     37 t=  0.30000000000E+01
! u  =   0.71430359E-02  0.11261442E+01
! ue =   0.71430219E-02  0.11261442E+01
! up =  -0.10022021E+00  0.14286050E-01 w  =  -0.20246660E+00
! upe=  -0.10022020E+00  0.14286044E-01 we =  -0.20246659E+00
!
! nstp =     44 t=  0.35000000000E+01
! u  =   0.44134944E-01  0.11375498E+01
! ue =   0.44134929E-01  0.11375498E+01
! up =   0.23564666E+00  0.88269850E-01 w  =   0.50327296E+00
! upe=   0.23564667E+00  0.88269858E-01 we =   0.50327297E+00
!
! nstp =     51 t=  0.40000000000E+01
! u  =   0.20543287E+00  0.12572816E+01
! ue =   0.20543286E+00  0.12572816E+01
! up =   0.35486108E+00  0.41086571E+00 w  =   0.10857938E+01
! upe=   0.35486108E+00  0.41086571E+00 we =   0.10857938E+01
!
! nstp =     57 t=  0.45000000000E+01
! u  =   0.34274025E+00  0.15401423E+01
! ue =   0.34274023E+00  0.15401423E+01
! up =   0.14781786E+00  0.68548048E+00 w  =   0.14024744E+01
! upe=   0.14781785E+00  0.68548047E+00 we =   0.14024744E+01
!
! nstp =     64 t=  0.50000000000E+01
! u  =   0.32981730E+00  0.18909544E+01
! ue =   0.32981729E+00  0.18909544E+01
! up =  -0.19512842E+00  0.65963457E+00 w  =   0.13757803E+01
! upe=  -0.19512843E+00  0.65963458E+00 we =   0.13757803E+01
!
! nstp =     70 t=  0.55000000000E+01
! u  =   0.17854534E+00  0.21520665E+01
! ue =   0.17854532E+00  0.21520665E+01
! up =  -0.35867452E+00  0.35709064E+00 w  =   0.10122473E+01
! upe=  -0.35867453E+00  0.35709064E+00 we =   0.10122473E+01
!
! nstp =     76 t=  0.60000000000E+01
! u  =   0.28003097E-01  0.22482967E+01
! ue =   0.28003078E-01  0.22482967E+01
! up =  -0.19245690E+00  0.56006151E-01 w  =   0.40088079E+00
! upe=  -0.19245693E+00  0.56006157E-01 we =   0.40088082E+00
!
! nstp =     82 t=  0.65000000000E+01
! u  =   0.16598437E-01  0.22560549E+01
! ue =   0.16598404E-01  0.22560549E+01
! up =   0.15070469E+00  0.33196803E-01 w  =  -0.30863527E+00
! upe=   0.15070469E+00  0.33196808E-01 we =  -0.30863527E+00
!
! nstp =     88 t=  0.70000000000E+01
! u  =   0.15481674E+00  0.23330917E+01
! ue =   0.15481670E+00  0.23330918E+01
! up =   0.35530911E+00  0.30963341E+00 w  =  -0.94258668E+00
! upe=   0.35530911E+00  0.30963341E+00 we =  -0.94258668E+00
!
! nstp =     95 t=  0.75000000000E+01
! u  =   0.31558075E+00  0.25734633E+01
! ue =   0.31558071E+00  0.25734634E+01
! up =   0.23324397E+00  0.63116141E+00 w  =  -0.13457600E+01
! upe=   0.23324397E+00  0.63116142E+00 we =  -0.13457600E+01
!
! nstp =    102 t=  0.80000000000E+01
! u  =   0.35108477E+00  0.29210566E+01
! ue =   0.35108474E+00  0.29210567E+01
! up =  -0.10326459E+00  0.70216947E+00 w  =  -0.14194443E+01
! upe=  -0.10326460E+00  0.70216948E+00 we =  -0.14194443E+01
!
! nstp =    108 t=  0.85000000000E+01
! u  =   0.22868659E+00  0.32211795E+01
! ue =   0.22868655E+00  0.32211795E+01
! up =  -0.34483217E+00  0.45737310E+00 w  =  -0.11455992E+01
! upe=  -0.34483217E+00  0.45737309E+00 we =  -0.11455992E+01
!
! nstp =    114 t=  0.90000000000E+01
! u  =   0.60918516E-01  0.33627837E+01
! ue =   0.60918470E-01  0.33627837E+01
! up =  -0.26936264E+00  0.12183694E+00 w  =  -0.59127141E+00
! upe=  -0.26936264E+00  0.12183694E+00 we =  -0.59127141E+00
!
! nstp =    120 t=  0.95000000000E+01
! u  =   0.20257594E-02  0.33805626E+01
! ue =   0.20257027E-02  0.33805626E+01
! up =   0.53757645E-01  0.40513995E-02 w  =   0.10782021E+00
! upe=   0.53757665E-01  0.40514055E-02 we =   0.10782023E+00
!
! nstp =    126 t=  0.10000000000E+02
! u  =   0.10615403E+00  0.34230537E+01
! ue =   0.10615398E+00  0.34230537E+01
! up =   0.32745341E+00  0.21230796E+00 w  =   0.78051371E+00
! upe=   0.32745342E+00  0.21230797E+00 we =   0.78051371E+00
! Run completed with info=  0
! Statistics
! total number of steps =  135
!        accepted steps =  126
!        rejected steps =  9
! local bases computed  =  11
! function evaluations  =  8496
! jacobian evaluations  =  832
!-----------------------------------------------------------------
