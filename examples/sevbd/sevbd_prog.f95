!--------------------------------------------------------------------
! Driver for applying the DAE-solver with eulag to the 
! seven body mechanism described in 
!
! Hairer, E., Wanner, G.: Solving Ordinary Differential Equations II.
! Berlin, Heidelberg: Springer 1991.
!
! In the notation of this reference the variables are
!
!  x(1) = beta,  x(2) = Theta, x(3) = gamma, x(4) = Phi,
!  x(5) = delta, x(6) = Omega, x(7) = epsilon
!--------------------------------------------------------------------
!
module prob_dat
  !
  ! problem data for the seven body mechanism
  !
  use def_kinds, ONLY : RP => RPX
  implicit none
  !
  ! location of fixed points
  !   
  real(kind=RP), parameter :: xpa=-0.06934_RP,  ypa=-0.00227_RP, &
                              xpb=-0.03635_RP,  ypb= 0.03273_RP, &
                              xpc= 0.014_RP,    ypc= 0.072_RP
  !
  ! geometrical parameters
  !
  real(kind=RP), parameter :: d  = 0.028_RP,    da = 0.0115_RP,  &
                              e  = 0.02_RP,     ea = 0.01421_RP, &
                              zf = 0.02_RP,     fa = 0.01421_RP, &
                              rr = 0.007_RP,    ra = 0.00092_RP, &
                              ss = 0.035_RP,    sa = 0.01874_RP, &
                              sb = 0.01043_RP,  sc = 0.018_RP,   &
                              sd = 0.02_RP,     zt = 0.04_RP,    &    
                              ta = 0.02308_RP,  tb = 0.00916_RP, &
                              u  = 0.04_RP,     ua = 0.01228_RP, &
                              ub = 0.00449_RP,  c0 = 4530.0_RP,  &
                              l0 = 0.07785_RP,  mom= 0.033_RP
  !
  ! masses of the bodies
  !
  real(kind=RP), parameter :: m1 = 0.04325_RP,  m2 = 0.00365_RP, &
                              m3 = 0.02373_RP,  m4 = 0.00706_RP, &
                              m5 = 0.07050_RP,  m6 = 0.00706_RP, & 
                              m7 = 0.05498_RP
  !
  ! inertias of the bodies
  !
  real(kind=RP), parameter :: i1 = 2.194e-6_RP, i2 = 4.41e-7_RP, &
                              i3 = 5.255e-6_RP, i4 = 5.667e-7_RP,&
                              i5 = 1.169e-5_RP, i6 = 5.667e-7_RP,&
                              i7 = 1.912e-5_RP
  !   
end module prob_dat
!
program sevbd
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
  character(len=5) :: dnm     = 'eulag'  ! dae type
  character(len=6) :: pname   = 'sevbd'  ! problem name
  integer          :: outunit = 10       ! output unit
  integer          :: nvd     = 7, &     ! no. of diff. equ.               
                      nva     = 6, &     ! no. of alg. equ.
                      nvu     = 7, &     ! dim. of u
                      nvw     = 6        ! dim of w
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
  !
  x(1) = -0.617138900142764496358948458001e-1_RP
  x(2) =  0.0_RP
  x(3) =  0.455279819163070380255912382449_RP
  x(4) =  0.222668390165885884674473185609_RP
  x(5) =  0.487364979543842550225598953530_RP
  x(6) = -0.222668390165885884674473185609_RP
  x(7) =  0.123054744454982119249735015568e+1_RP
  !
  ! set global data
  !
  reltol = 1.0e-8_RP                     ! relative tolerance
  abstol = reltol*1.0e-2_RP              ! absolute tolerance
  h    = 1.0e-3_RP                       ! first step
  hint = 0.0025_RP                       ! interpolation step
  tout = 0.03_RP                         ! terminal time
  hmax = tout                            ! maximum step
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
end program sevbd
!
subroutine daefct(fname,xc,vc,vec,mat,iret)
  !
  use def_kinds, ONLY : RP => RPX
  use dae_dat
  use prob_dat
  implicit none
  intrinsic sin, cos, sqrt
  !
  character(len=4),                intent(in)  :: fname
  real(kind=RP), dimension(0:),    intent(in)  :: xc
  real(kind=RP), dimension(0:),    intent(in)  :: vc
  real(kind=RP), dimension(0:),    intent(out) :: vec
  real(kind=RP), dimension(0:,0:), intent(out) :: mat
  integer,                         intent(out) :: iret
  optional vc,vec,mat
  !
  real(kind=RP) :: xp1,xp2,xp3,xp4,xp5,xp6,xp7, &
                   xpd,ypd,lang,force,fx,fy, &
                   si1,si2,si3,si4,si5,si6,si7, &
                   co1,co2,co3,co4,co5,co6,co7, &
                   co12,si12,co45,si45,co67,si67, &
                   v1sq,v3sq,v5sq,v7sq,v12sq,v45sq,v67sq
  !
  real(kind=RP), parameter :: zer=0.0_RP
  !--------------------------------------------------------
  ! Routine for evaluating the parts of the dae under
  ! control of the identifier fname.
  ! iret = 0   no error
  ! iret = 1   not available part as requested by fname
  ! iret = 2   wrong dimension for output
  !--------------------------------------------------------
  !
  iret = 0 
  select case(fname)
    !
    case ('amt ')
      if(.not.present(mat)) then
        iret = 2
        return
      endif
      !
      co2 = cos(xc(2))
      si4 = sin(xc(4))
      si6 = sin(xc(6))
      !   
      mat(0:6,0:6) = zer
      !
      mat(0,0) = m1*ra*ra + m2*(rr*(rr-da*co2) + da*(da-rr*co2)) &
                 + i1 + i2
      mat(0,1) = m2*da*(da - rr*co2) + i2
      mat(1,1) = m2*da*da + i2
      mat(2,2) = m3*(sa*sa + sb*sb) + i3
      mat(3,3) = m4*(e-ea)**2 + i4
      mat(3,4) = m4*(e-ea)*((e-ea) + zt*si4) + i4
      mat(4,4) = m4*(zt*(zt+(e-ea)*si4) + (e-ea)*((e-ea)+zt*si4)) &
                 + m5*(ta*ta + tb*tb) + i4 + i5
      mat(5,5) = m6*(zf-fa)**2 + i6
      mat(5,6) = m6*(zf-fa)*((zf-fa)-u*si6) + i6
      mat(6,6) = m6*((zf-fa)*((zf-fa)-u*si6) + u*(u-(zf-fa)*si6)) &
                 + m7*(ua*ua + ub*ub) + i6 + i7
      !
    case ('hf  ')
      if(.not.present(vec)) then
        iret = 2
        return
      endif
      ! 
      si2 = sin(xc(2))
      si3 = sin(xc(3))
      co3 = cos(xc(3))
      co4 = cos(xc(4))
      co6 = cos(xc(6))   
      ! 
      xp1 = xc(8)
      xp2 = xc(9)
      xp3 = xc(10)
      xp4 = xc(11)
      xp5 = xc(12)
      xp6 = xc(13)
      xp7 = xc(14)
      !
      xpd   = sd*co3 + sc*si3 + xpb
      ypd   = sd*si3 - sc*co3 + ypb
      lang  = sqrt((xpd-xpc)**2 + (ypd-ypc)**2)
      force = - c0*(lang - l0)/lang
      fx    = force*(xpd-xpc)
      fy    = force*(ypd-ypc)
      !
      vec(0) = mom - m2*da*rr*xp2*(xp2 + 2.0_RP*xp1)*si2
      vec(1) = m2*da*rr*xp1*xp1*si2
      vec(2) = fx*(sc*co3 - sd*si3) + fy*(sd*co3 + sc*si3)
      vec(3) = m4*zt*(e-ea)*xp5*xp5*co4
      vec(4) = - m4*zt*(e-ea)*xp4*(xp4 + 2.0_RP*xp5)*co4
      vec(5) = - m6*u*(zf-fa)*xp7*xp7*co6
      vec(6) = m6*u*(zf-fa)*xp6*(xp6 + 2.0_RP*xp7)*co6
      !
    case ('ff  ')
      if(.not.present(vec)) then
        iret = 2
        return
      endif
      !
      si1 = sin(xc(1))
      co1 = cos(xc(1))
      si3 = sin(xc(3))
      co3 = cos(xc(3))
      si5 = sin(xc(5))
      co5 = cos(xc(5))
      si7 = sin(xc(7))
      co7 = cos(xc(7))
      !
      co12 = cos(xc(1)+xc(2))
      si12 = sin(xc(1)+xc(2))
      co45 = cos(xc(4)+xc(5))
      si45 = sin(xc(4)+xc(5))
      co67 = cos(xc(6)+xc(7))
      si67 = sin(xc(6)+xc(7))
      !
      vec(0) = rr*co1 - d*co12 - ss*si3 - xpb
      vec(1) = rr*si1 - d*si12 + ss*co3 - ypb
      vec(2) = rr*co1 - d*co12 - e*si45 - zt*co5 - xpa
      vec(3) = rr*si1 - d*si12 + e*co45 - zt*si5 - ypa
      vec(4) = rr*co1 - d*co12 - zf*co67 - u*si7 - xpa
      vec(5) = rr*si1 - d*si12 - zf*si67 + u*co7 - ypa
      !
    case ('dff ')
      if(.not.present(mat)) then
        iret = 2
        return
      endif
      !
      si1 = sin(xc(1))
      co1 = cos(xc(1))
      si3 = sin(xc(3))
      co3 = cos(xc(3))
      si5 = sin(xc(5))
      co5 = cos(xc(5))
      si7 = sin(xc(7))
      co7 = cos(xc(7))
      !
      co12 = cos(xc(1)+xc(2))
      si12 = sin(xc(1)+xc(2))
      co45 = cos(xc(4)+xc(5))
      si45 = sin(xc(4)+xc(5))
      co67 = cos(xc(6)+xc(7))
      si67 = sin(xc(6)+xc(7))
!
      mat(0:5,0:7) = zer
      !
      mat(0,1) = - rr*si1 + d*si12
      mat(0,2) = d*si12
      mat(0,3) = - ss*co3
      mat(1,1) = rr*co1 - d*co12
      mat(1,2) = - d*co12
      mat(1,3) = - ss*si3
      mat(2,1) = - rr*si1 + d*si12
      mat(2,2) = d*si12
      mat(2,4) = - e*co45
      mat(2,5) = - e*co45 + zt*si5
      mat(3,1) = rr*co1 - d*co12
      mat(3,2) = - d*co12
      mat(3,4) = - e*si45
      mat(3,5) = - e*si45 - zt*co5
      mat(4,1) = - rr*si1 + d*si12
      mat(4,2) = d*si12
      mat(4,6) = zf*si67
      mat(4,7) = zf*si67 - u*co7
      mat(5,1) = rr*co1 - d*co12
      mat(5,2) = - d*co12
      mat(5,6) = - zf*co67
      mat(5,7) = - zf*co67 - u*si7
      !
    case ('d2ff')
      if(.not.present(vc) .or. .not.present(vec)) then
        iret = 2
        return
      endif
      !
      co1 = cos(xc(1))
      si1 = sin(xc(1))
      co3 = cos(xc(3))
      si3 = sin(xc(3))
      co5 = cos(xc(5))
      si5 = sin(xc(5))
      co7 = cos(xc(7))
      si7 = sin(xc(7))
      !
      co12 = cos(xc(1)+xc(2))
      si12 = sin(xc(1)+xc(2))
      co45 = cos(xc(4)+xc(5))
      si45 = sin(xc(4)+xc(5))
      co67 = cos(xc(6)+xc(7))
      si67 = sin(xc(6)+xc(7))
      !
      v1sq  = vc(1)**2
      v3sq  = vc(3)**2
      v5sq  = vc(5)**2
      v7sq  = vc(7)**2
      v12sq = (vc(1)+vc(2))**2
      v45sq = (vc(4)+vc(5))**2   
      v67sq = (vc(4)+vc(5))**2
      !
      vec(0) = - rr*co1*v1sq   + d*co12*v12sq  + ss*si3*v3sq
      vec(1) = - rr*si1*v1sq   + d*si12*v12sq  - ss*co3*v3sq
      vec(2) = - rr*co1*v1sq   + d*co12*v12sq  &
                               + e*si45*v45sq  + zt*co5*v5sq
      vec(3) = - rr*si1*v1sq   + d*si12*v12sq &
                               - e*co45*v45sq  + zt*si5*v5sq
      vec(4) = - rr*co1*v1sq   + d*co12*v12sq &
                               + zf*co67*v67sq + u*si7*v7sq
      vec(5) = - rr*si1*v1sq   + d*si12*v12sq &
                               + zf*si67*v67sq - u*co7*v7sq
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
  intrinsic atan
  !
  integer,                      intent(in)  :: nstp
  real(kind=RP), dimension(0:), intent(in)  :: xc
  real(kind=RP),                intent(out) :: hinch
  integer,                      intent(out) :: info
  !
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
  write(outlun,10) nstp,xc(0)
  10 format(1X/' nstp = ',I6,' t= ',e18.11)
  !
  write(outlun,20)xc(1:7)
  20 format(' xc  = '/(4e16.8))
  write(outlun,30)xc(8:14)
  30 format(' xp  = '/(4e16.8))
  write(outlun,40)xc(15:21)
  40 format(' xpp = '/(4e16.8))
  write(outlun,50)xc(22:27)
  50 format(' lambda = '/(4e16.8))
  !
end subroutine solout
!
!-------------------------------------------------------------
!
! Output from this program:
!
!
! log file sevbd.log
!
! open_log :: info=   0 opened  sevbd.log
! open_file :: info=   0 opened  sevbd.out
! Basis check: h=  0.96482170E-04 hold=  0.10000000E-02
!   compute new basis at nstep=      4
! Basis check: itstp=  23
!   compute new basis at nstep=     48
! eval_point :: info= -21 divergence detected
! daeslv :: info= -21 global point evaluation failed - reduce step
! eval_point :: info= -21 divergence detected
! daeslv :: info= -21 global point evaluation failed - reduce step
! eval_point :: info= -21 divergence detected
! daeslv :: info= -21 global point evaluation failed - reduce step
! eval_point :: info= -21 divergence detected
!   compute new basis at nstep=     96
! Basis check: itstp=  21
!   compute new basis at nstep=    113
! eval_point :: info= -21 divergence detected
! daeslv :: info= -21 global point evaluation failed - reduce step
! eval_point :: info= -21 divergence detected
! daeslv :: info= -21 global point evaluation failed - reduce step
! eval_point :: info= -21 divergence detected
! daeslv :: info= -21 global point evaluation failed - reduce step
! eval_point :: info= -21 divergence detected
!   compute new basis at nstep=    158
! main :: info=   0 Run completed
! close_file :: info=   0  closed  sevbd.out
!
!
! output file sevbd.out
!
! problem name: sevbd 
! DAE type    : eulag
! nv=  7 na=  6 nu=  7 nw=  6
! tolerances: reltol=  0.100000E-07 abstol=   0.100000E-09
! h=   0.100000E-02 tout=   0.300000E-01
!
! nstp =      0 t=  0.00000000000E+00
! xc  = 
! -0.61713890E-01  0.00000000E+00  0.45527982E+00  0.22266839E+00
!  0.48736498E+00 -0.22266839E+00  0.12305474E+01
! xp  = 
!  0.00000000E+00  0.00000000E+00  0.00000000E+00  0.00000000E+00
!  0.00000000E+00  0.00000000E+00  0.00000000E+00
! xpp = 
!  0.11816556E+05 -0.88624168E+04  0.00000000E+00  0.00000000E+00
!  0.00000000E+00  0.00000000E+00  0.00000000E+00
! lambda = 
!  0.98698596E+02 -0.63917453E+01  0.00000000E+00  0.00000000E+00
!  0.00000000E+00  0.00000000E+00
!
! nstp =     14 t=  0.25000000000E-02
! xc  = 
! -0.21086646E-01 -0.30385589E-01  0.45513825E+00  0.22238250E+00
!  0.48739698E+00 -0.22238250E+00  0.12304393E+01
! xp  = 
!  0.35559810E+02 -0.26524086E+02 -0.24710583E+00 -0.49896213E+00
!  0.55839904E-01  0.49896213E+00 -0.18864563E+00
! xpp = 
!  0.19213487E+05 -0.14331353E+05 -0.13351838E+03 -0.26959685E+03
!  0.30174415E+02  0.26959672E+03 -0.10193477E+03
! lambda = 
!  0.97999126E+02 -0.53669515E+01 -0.11866125E+00  0.10179702E+00
! -0.16607141E+00 -0.23933107E+00
!
! nstp =     25 t=  0.50000000000E-02
! xc  = 
!  0.15070512E+00 -0.15736183E+00  0.45150726E+00  0.21506175E+00
!  0.48821229E+00 -0.21506175E+00  0.12276767E+01
! xp  = 
!  0.11349901E+03 -0.83281900E+02 -0.39664317E+01 -0.79851186E+01
!  0.88497387E+00  0.79851186E+01 -0.30077522E+01
! xpp = 
!  0.46936067E+05 -0.34440165E+05 -0.16402999E+04 -0.33021465E+04
!  0.36599950E+03  0.33021453E+04 -0.12438795E+04
! lambda = 
!  0.93351756E+02 -0.27817502E+00 -0.12146646E+01  0.10475203E+01
! -0.17064245E+01 -0.24737384E+01
!
! nstp =     34 t=  0.75000000000E-02
! xc  = 
!  0.62617129E+00 -0.50617022E+00  0.41848874E+00  0.14944218E+00
!  0.49516498E+00 -0.14944218E+00  0.12034023E+01
! xp  = 
!  0.28402688E+03 -0.21133158E+03 -0.28815161E+02 -0.56561513E+02
!  0.57164386E+01  0.56561513E+02 -0.20520146E+02
! xpp = 
!  0.90620261E+05 -0.67426459E+05 -0.91937866E+04 -0.18046262E+05
!  0.18239504E+04  0.18046257E+05 -0.65472715E+04
! lambda = 
!  0.76042755E+02  0.12094783E+02 -0.48467887E+01  0.43116977E+01
! -0.70106676E+01 -0.10644679E+02
!
! nstp =     42 t=  0.10000000000E-01
! xc  = 
!  0.16339143E+01 -0.13421300E+01  0.26278455E+00 -0.14277476E+00
!  0.51810214E+00  0.14277476E+00  0.11097833E+01
! xp  = 
!  0.53345760E+03 -0.50699464E+03 -0.10395801E+03 -0.18815921E+03
!  0.10320679E+02  0.18815921E+03 -0.50441080E+02
! xpp = 
!  0.15457059E+06 -0.14690308E+06 -0.30122197E+05 -0.54519635E+05
!  0.29900489E+04  0.54519640E+05 -0.14614597E+05
! lambda = 
!  0.52600115E+02  0.95934063E+01 -0.12208566E+01 -0.34658643E+01
! -0.88961577E+00  0.53397801E+01
!
! nstp =     58 t=  0.12500000000E-01
! xc  = 
!  0.33875600E+01 -0.33686125E+01  0.44749516E-01 -0.52786319E+00
!  0.52460748E+00  0.52786319E+00  0.10482168E+01
! xp  = 
!  0.81978140E+03 -0.10206586E+04  0.36914725E+02  0.64539583E+02
!  0.17942102E+01 -0.64539583E+02  0.14752331E+01
! xpp = 
! -0.17796776E+06  0.22157684E+06 -0.80111660E+04 -0.14010798E+05
! -0.38904699E+03  0.14010798E+05 -0.31984558E+03
! lambda = 
!  0.12376127E+03 -0.19854167E+02  0.11324782E+02  0.12818871E+02
!  0.11239431E+02  0.11523478E+02
!
! nstp =     71 t=  0.15000000000E-01
! xc  = 
!  0.50213913E+01 -0.52458612E+01  0.32516307E+00 -0.28466306E-01
!  0.51073962E+00  0.28466306E-01  0.11430517E+01
! xp  = 
!  0.48996989E+03 -0.49229128E+03  0.98951289E+02  0.18381308E+03
! -0.13552577E+02 -0.18381308E+03  0.57402890E+02
! xpp = 
! -0.11176233E+06  0.11229196E+06 -0.22571151E+05 -0.41927892E+05
!  0.30911193E+04  0.41927892E+05 -0.13093232E+05
! lambda = 
!  0.61489850E+02 -0.62106275E+01 -0.33353801E+01  0.39512229E+00
! -0.43956629E+01 -0.36166239E+01
!
! nstp =     80 t=  0.17500000000E-01
! xc  = 
!  0.60495979E+01 -0.61524727E+01  0.45267141E+00  0.21740651E+00
!  0.48795202E+00 -0.21740651E+00  0.12285604E+01
! xp  = 
!  0.36680666E+03 -0.28330692E+03  0.11254505E+02  0.22678976E+02
! -0.25213404E+01 -0.22678976E+02  0.85527667E+01
! xpp = 
!  0.50411552E+04 -0.38935719E+04  0.15426873E+03  0.31165206E+03
! -0.34275870E+02 -0.31166754E+03  0.11676654E+03
! lambda = 
!  0.80390562E+02  0.41546900E+01 -0.83085823E+01  0.71412832E+01
! -0.11657933E+02 -0.16853179E+02
!
! nstp =     92 t=  0.20000000000E-01
! xc  = 
!  0.70459034E+01 -0.68918829E+01  0.40349359E+00  0.12016192E+00
!  0.49805998E+00 -0.12016192E+00  0.11928804E+01
! xp  = 
!  0.46518586E+03 -0.35283453E+03 -0.54898101E+02 -0.10664700E+03
!  0.10309914E+02  0.10664700E+03 -0.37948277E+02
! xpp = 
!  0.94079869E+05 -0.71357831E+05 -0.11103055E+05 -0.21568490E+05
!  0.20852472E+04  0.21568481E+05 -0.76751125E+04
! lambda = 
!  0.65244564E+02  0.15571947E+02 -0.72413278E+01  0.61910365E+01
! -0.10527300E+02 -0.15877815E+02
!
! nstp =    102 t=  0.22500000000E-01
! xc  = 
!  0.85063202E+01 -0.82367334E+01  0.14658167E+00 -0.34964058E+00
!  0.52552735E+00  0.34964058E+00  0.10649699E+01
! xp  = 
!  0.74091552E+03 -0.83133649E+03 -0.13850714E+03 -0.24347361E+03
!  0.37797794E+01  0.24347361E+03 -0.38462632E+02
! xpp = 
!  0.21115230E+06 -0.23692159E+06 -0.39472102E+05 -0.69387149E+05
!  0.10762567E+04  0.69387159E+05 -0.10958966E+05
! lambda = 
!  0.67950229E+02 -0.14340718E+02  0.46847287E+01 -0.20378356E+02
!  0.86723403E+01  0.51368232E+02
!
! nstp =    121 t=  0.25000000000E-01
! xc  = 
!  0.10589184E+02 -0.10756890E+02  0.17157649E+00 -0.30561377E+00
!  0.52464109E+00  0.30561377E+00  0.10725313E+01
! xp  = 
!  0.72434598E+03 -0.82925823E+03  0.14927154E+03  0.26351830E+03
! -0.64963770E+01 -0.26351830E+03  0.48788586E+02
! xpp = 
! -0.25204300E+06  0.28854851E+06 -0.51939523E+05 -0.91693687E+05
!  0.22593637E+04  0.91693701E+05 -0.16973589E+05
! lambda = 
!  0.73355921E+02 -0.44171574E+02  0.55683109E+01 -0.22035611E+02
!  0.10716224E+02  0.54972370E+02
!
! nstp =    132 t=  0.27500000000E-01
! xc  = 
!  0.12045244E+02 -0.12206308E+02  0.43606653E+00  0.18416869E+00
!  0.49156533E+00 -0.18416869E+00  0.12161353E+01
! xp  = 
!  0.51546954E+03 -0.42573761E+03  0.43924402E+02  0.87349957E+02
! -0.92803555E+01 -0.87349957E+02  0.32356723E+02
! xpp = 
! -0.11887263E+04  0.98187289E+03 -0.10212913E+03 -0.20151183E+03
!  0.22042033E+02  0.20148338E+03 -0.76007690E+02
! lambda = 
!  0.59832458E+02  0.15113424E+02 -0.14597172E+02  0.12630914E+02
! -0.20796406E+02 -0.30613842E+02
!
! nstp =    146 t=  0.30000000000E-01
!
! xc  = 
!  0.13374173E+02 -0.13209393E+02  0.39805436E+00  0.10961459E+00
!  0.49907125E+00 -0.10961459E+00  0.11891411E+01
! xp  = 
!  0.57565777E+03 -0.43998153E+03 -0.70956720E+02 -0.13734656E+03
!  0.13059397E+02  0.13734656E+03 -0.48510878E+02
! xpp = 
!  0.96535432E+05 -0.73783190E+05 -0.11899711E+05 -0.23032524E+05
!  0.21901976E+04  0.23032510E+05 -0.81355891E+04
! lambda = 
!  0.58425040E+02  0.16903271E+02 -0.87861430E+01  0.72338846E+01
! -0.12754839E+02 -0.18960466E+02
! Run completed with info=  0
! Statistics
! total number of steps =  163
!        accepted steps =  146
!        rejected steps =  12
! local bases computed  =  6
! function evaluations  =  9042
! jacobian evaluations  =  955
!--------------------------------------------------------------------