!-----------------------------------------------------------------
!  Driver for applying the DAE-solver with daeq1 for solving
!  the catalytic converter problem
!
!  u1' = -jd*(u1 - u6)
!  u2' = -jh*(u2 - u5)
!  u3' = u4 
!  u4' = cc1*(2.0*u3 - u5 - theta)
!   0  = jh*(u2 -u5) + jk*(u3 -u5) 
!                     - cc3*u6*exp(gc*(1.0 - 1.0/u5))
!   0  = u6*(1.0 + cc2*exp(gc*(1.0 - 1.0/u5))) - u1 
!
!-----------------------------------------------------------------
!
module prob_dat
  !
  ! problem data for the catalytic converter problem
  !
  use def_kinds, ONLY : RP => RPX
  implicit none
  !
  real(kind=RP), parameter :: jk    = 143.8_RP, &
                              jd    = 2.63_RP,  &
                              jh    = 5.26_RP,  &
                              gc    = 10.0_RP,  &
                              theta =  1.5_RP,  &
                              gamma =  0.478_RP,&
                              dc    =  0.249_RP,&
                              pe    =  0.2_RP
  real(kind=RP), parameter :: cc1 = pe*jk, &
                              cc2 = dc/jd, &
                              cc3 = gamma*dc
  !
end module prob_dat
!
program catal
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
  character(len=5) :: dnm     = 'daeq1'  ! dae type
  character(len=6) :: pname   = 'catal'  ! problem name
  integer          :: outunit = 10       ! output unit
  integer          :: nvd     = 4, &     ! no. of diff. equ.               
                      nva     = 2, &     ! no. of alg. equ.
                      nvu     = 6, &     ! dim. of u
                      nvw     = 0        ! dim of w
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
  x(0:kx) = initpt()
  !
  ! set global data
  !
  reltol = 1.0e-8_RP                     ! relative tolerance
  abstol = reltol*1.0e-2_RP              ! absolute tolerance
  h    = 1.0e-2_RP                       ! first step
  hint = 0.05_RP                         ! interpolation step
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
contains
  !
  function initpt() result(x0)
    !    
    use def_kinds, ONLY : RP => RPX
    implicit none
    real(kind=RP), dimension(0:kx) :: x0
    intrinsic abs, exp
    !
    integer                  :: it, itmax=20
    real(kind=RP)            :: step,w1,w2,ew1,dew1
    real(kind=RP), parameter :: zer=0.0_RP, one=1.0_RP
    real(kind=RP), parameter :: epmach=epsilon(one)
    real(kind=RP), parameter :: tol = epmach*1.0e-3_RP
    !
    x0(0:kx) = zer
    x0(1) = one
    x0(2) = one
    x0(3) = 1.486788276286347_RP
    !
    ! find x0(5) via newton's method
    !
    w1   = one
    step = one
    it   = 0
    !
    do while (abs(step) > tol .and. it < itmax)
      !
      ew1  = exp(gc*(one-one/w1))
      dew1 = ew1*gc/(w1*w1)
      w2   = x0(1)/(one + cc2*ew1)
      step = (jh*(x0(2)-w1)+ &
              jk*(x0(3)-w1)-cc3*w2*ew1)/(jh+jk+cc3*w2*dew1)
      w1   = w1 + step
      it   = it + 1
      !
    enddo
    !
    x0(5) = w1
    x0(6) = w2
    !
  end function initpt
  !
end program catal
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
  real(kind=RP) :: ew1,dew1
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
      mat(0:3,0:5) = 0.0_RP
      mat(0,0) = 1.0_RP
      mat(1,1) = 1.0_RP
      mat(2,2) = 1.0_RP
      mat(3,3) = 1.0_RP
      !
    case ('hf  ')
      if(.not.present(vec)) then
        iret = 2
        return
      endif
      vec(0) = jd*(xc(6) - xc(1))
      vec(1) = jh*(xc(5) - xc(2))
      vec(2) = xc(4)
      vec(3) = cc1*(2.0_RP*xc(3) - xc(5) - theta)
      !
    case ('ff  ')
      if(.not.present(vec)) then
        iret = 2
        return
      endif
      ew1  = exp(gc*(1.0_RP - 1.0_RP/xc(5)))
      vec(0) = jh*(xc(2)-xc(5)) + jk*(xc(3)-xc(5)) &
                                - cc3*xc(6)*ew1
      vec(1) = xc(6)*(1.0_RP + cc2*ew1) - xc(1)
      !
    case ('dff ')
      if(.not.present(mat)) then
        iret = 2
        return
      endif
      mat = 0.0_RP
      ew1 = exp(gc*(1.0_RP - 1.0_RP/xc(5)))
      dew1 = ew1*gc/(xc(5)*xc(5))
      !
      mat(0,2) = jh
      mat(0,3) = jk
      mat(0,5) =-(jh + jk) - cc3*xc(6)*dew1
      mat(0,6) = -cc3*ew1
      mat(1,1) = -1.0_RP
      mat(1,5) = xc(6)*cc2*dew1
      mat(1,6) = 1.0_RP + cc2*ew1
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
  integer,                      intent(in)  :: nstp
  real(kind=RP), dimension(0:), intent(in)  :: xc
  real(kind=RP),                intent(out) :: hinch
  integer,                      intent(out) :: info
  ! 
  !-----------------------------------------------------------
  ! Routine for intermediate printout of the solution vector
  ! x = (t,u(0:5),up(0:5)).
  !
  ! nstp   current count of accepted steps
  ! hinch  a new value of the interpolation step or
  !        zero if no change is desired
  ! info   nonzero if integration is to stop or
  !        zero otherwise
  !-----------------------------------------------------------
  write(outlun,10) nstp,xc(0)
  10 format(1X/' nstp = ',I6,' t= ',e18.11)
  !
  write(outlun,20)xc(1:6)
  20 format(' u = '/(4e16.8))
  write(outlun,30)xc(7:12)
  30 format(' up= '/(4e16.8))
  !
end subroutine solout
!
!-------------------------------------------------------------
!
! Output from this program:
!
!
! log file catal.log
!
! open_log :: info=   0 opened  catal.log
! open_file :: info=   0 opened  catal.out
! main :: info=   0 Run completed
! close_file :: info=   0  closed  catal.out
!
!
! output file catal.out
!
! problem name: catal 
! DAE type    : daeq1
! nv=  4 na=  2 nu=  6 nw=  0
! tolerances: reltol=  0.100000E-07 abstol=   0.100000E-09
! h=   0.100000E-01 tout=   0.100000E+01
!
! nstp =      0 t=  0.00000000000E+00
! u = 
!  0.10000000E+01  0.10000000E+01  0.14867883E+01  0.00000000E+00
!  0.14637716E+01  0.30766746E+00
! up= 
! -0.18208346E+01  0.24394386E+01  0.11102230E-15  0.28199034E+00
!  0.95910106E-01 -0.65556009E+00
!
! nstp =      5 t=  0.50000000000E-01
! u = 
!  0.91270378E+00  0.11078112E+01  0.14870885E+01  0.11059107E-01
!  0.14683406E+01  0.27669337E+00
! up= 
! -0.16727074E+01  0.18963848E+01  0.11059107E-01  0.16785242E+00
!  0.86764241E-01 -0.58468707E+00
!
! nstp =      8 t=  0.10000000000E+00
! u = 
!  0.83256257E+00  0.11916909E+01  0.14878162E+01  0.17415071E-01
!  0.14724455E+01  0.24907117E+00
! up= 
! -0.15345824E+01  0.14767691E+01  0.17415071E-01  0.91652899E-01
!  0.77456997E-01 -0.52145135E+00
!
! nstp =     11 t=  0.15000000000E+00
! u = 
!  0.75908403E+00  0.12570665E+01  0.14887785E+01  0.20666933E-01
!  0.14760927E+01  0.22442753E+00
! up= 
! -0.14061466E+01  0.11520779E+01  0.20666933E-01  0.42113222E-01
!  0.68525574E-01 -0.46545017E+00
!
! nstp =     13 t=  0.20000000000E+00
! u = 
!  0.69179209E+00  0.13081143E+01  0.14898500E+01  0.21935974E-01
!  0.14793097E+01  0.20241635E+00
! up= 
! -0.12870582E+01  0.90048766E+00  0.21935974E-01  0.11226307E-01
!  0.60284156E-01 -0.41603498E+00
!
! nstp =     16 t=  0.25000000000E+00
! u = 
!  0.63022883E+00  0.13480528E+01  0.14909522E+01  0.22005569E-01
!  0.14821355E+01  0.18272636E+00
! up= 
! -0.11769315E+01  0.70527524E+00  0.22005569E-01 -0.66436347E-02
!  0.52901217E-01 -0.37248099E+00
!
! nstp =     18 t=  0.30000000000E+00
! u = 
!  0.57395668E+00  0.13793660E+01  0.14920398E+01  0.21422623E-01
!  0.14846154E+01  0.16508245E+00
! up= 
! -0.10753392E+01  0.55361196E+00  0.21422623E-01 -0.15411391E-01
!  0.46454062E-01 -0.33407705E+00
!
! nstp =     20 t=  0.35000000000E+00
! u = 
!  0.52256030E+00  0.14039743E+01  0.14930901E+01  0.20571288E-01
!  0.14867969E+01  0.14924376E+00
! up= 
! -0.98182251E+00  0.43564684E+00  0.20571288E-01 -0.17737284E-01
!  0.40967313E-01 -0.30016882E+00
!
! nstp =     22 t=  0.40000000000E+00
! u = 
!  0.47564781E+00  0.14233652E+01  0.14940971E+01  0.19726925E-01
!  0.14887281E+01  0.13500036E+00
! up= 
! -0.89590280E+00  0.34380925E+00  0.19726925E-01 -0.15358596E-01
!  0.36440327E-01 -0.27017557E+00
!
! nstp =     24 t=  0.45000000000E+00
! u = 
!  0.43285151E+00  0.14386930E+01  0.14950664E+01  0.19096106E-01
!  0.14904568E+01  0.12216940E+00
! up= 
! -0.81709394E+00  0.27227811E+00  0.19096106E-01 -0.93214349E-02
!  0.32867128E-01 -0.24359290E+00
!
! nstp =     26 t=  0.50000000000E+00
! u = 
!  0.39382798E+00  0.14508562E+01  0.14960130E+01  0.18847013E-01
!  0.14920308E+01  0.11059147E+00
! up= 
! -0.74491202E+00  0.21657810E+00  0.18847013E-01 -0.13600449E-03
!  0.30251603E-01 -0.21998825E+00
!
! nstp =     27 t=  0.55000000000E+00
! u = 
!  0.35825776E+00  0.14605569E+01  0.14969600E+01  0.19133557E-01
!  0.14934984E+01  0.10012711E+00
! up= 
! -0.67888362E+00  0.17327226E+00  0.19133557E-01  0.12125058E-01
!  0.28620063E-01 -0.19899345E+00
!
! nstp =     29 t=  0.60000000000E+00
! u = 
!  0.32584468E+00  0.14683458E+01  0.14979380E+01  0.20115800E-01
!  0.14949102E+01  0.90653859E-01
! up= 
! -0.61855187E+00  0.13972849E+00  0.20115800E-01  0.27775270E-01
!  0.28032822E-01 -0.18029639E+00
!
! nstp =     30 t=  0.65000000000E+00
! u = 
!  0.29631490E+00  0.14746588E+01  0.14989862E+01  0.21978796E-01
!  0.14963208E+01  0.82063603E-01
! up= 
! -0.56348092E+00  0.11394240E+00  0.21978796E-01  0.47501208E-01
!  0.28596222E-01 -0.16363295E+00
!
! nstp =     32 t=  0.70000000000E+00
! u = 
!  0.26941579E+00  0.14798440E+01  0.15001543E+01  0.24951659E-01
!  0.14977917E+01  0.74260411E-01
! up= 
! -0.51325866E+00  0.94404948E-01  0.24951659E-01  0.72386472E-01
!  0.30476450E-01 -0.14877967E+00
!
! nstp =     34 t=  0.75000000000E+00
! u = 
!  0.24491469E+00  0.14841847E+01  0.15015047E+01  0.29328775E-01
!  0.14993944E+01  0.67158634E-01
! up= 
! -0.46749842E+00  0.80003203E-01  0.29328775E-01  0.10397003E+00
!  0.33916568E-01 -0.13554732E+00
!
! nstp =     35 t=  0.80000000000E+00
! u = 
!  0.22259761E+00  0.14879166E+01  0.15031169E+01  0.35494780E-01
!  0.15012150E+01  0.60681350E-01
! up= 
! -0.42583976E+00  0.69949329E-01  0.35494780E-01  0.14434320E+00
!  0.39258152E-01 -0.12377530E+00
!
! nstp =     37 t=  0.85000000000E+00
! u = 
!  0.20226794E+00  0.14912434E+01  0.15050924E+01  0.43955516E-01
!  0.15033596E+01  0.54759056E-01
! up= 
! -0.38794836E+00  0.63731339E-01  0.43955516E-01  0.19629102E+00
!  0.46969461E-01 -0.11332668E+00
!
! nstp =     39 t=  0.90000000000E+00
! u = 
!  0.18374514E+00  0.14943491E+01  0.15075617E+01  0.55377389E-01
!  0.15059618E+01  0.49328597E-01
! up= 
! -0.35351552E+00  0.61082852E-01  0.55377389E-01  0.26348773E+00
!  0.57682377E-01 -0.10408358E+00
!
! nstp =     41 t=  0.95000000000E+00
! u = 
!  0.16686350E+00  0.14974105E+01  0.15106939E+01  0.70638047E-01
!  0.15091918E+01  0.44332335E-01
! up= 
! -0.32225695E+00  0.61969760E-01  0.70638047E-01  0.35075962E+00
!  0.72240849E-01 -0.95942558E-01
!
! nstp =     43 t=  0.10000000000E+01
! u = 
!  0.15147091E+00  0.15006082E+01  0.15147086E+01  0.90892323E-01
!  0.15132685E+01  0.39717543E-01
! up= 
! -0.29391135E+00  0.66592894E-01  0.90892323E-01  0.46443450E+00
!  0.91764634E-01 -0.88809610E-01
! Run completed with info=  0
! Statistics
! total number of steps =  43
!        accepted steps =  43
!        rejected steps =  0
! local bases computed  =  1
! function evaluations  =  1763
! jacobian evaluations  =  260
!-------------------------------------------------------------------
