!-----------------------------------------------------------------
! File containing the following three library modules for use with
! the dae package:
!
! module def_kinds    to define kind-parameters for
!                     integer and real arithmetic
!
! module file_sup     programs for handlingfile support and, 
!                     in particular, the log file
! module lin_lib      a small linear algebra library
!------------------------------------------------------------------
!
module def_kinds
  !
  implicit none
  public
  integer, parameter ::   I1B = selected_int_kind(2), &
                          I2B = selected_int_kind(4), &
                          I4B = selected_int_kind(9), &
                          I8B = selected_int_kind(18)
  !
  integer, parameter :: R4B = selected_real_kind(P =  6, R =  37), &
                        R8B = selected_real_kind(P = 15, R = 307), &
                        R16B = selected_real_kind(P = 33, R = 4931)
  !
  integer, parameter :: RPX = R8B
  !
end module def_kinds
!
module file_sup
!
integer, parameter :: name_len = 12      ! max length of file names
integer, parameter :: loglun = 4         ! the unit of the logfile
integer, save      :: lunit              ! the unit of the outfile
!
!
contains
  !
  subroutine open_log(name)
    !
    implicit none
    !
    character(len=*), intent(in) :: name
    !
    open(unit=loglun, file=trim(name)//'.log', status='REPLACE',action= 'WRITE')
    !
  end subroutine open_log
  !
  subroutine open_file(lunit,name)
    !
    implicit none
    !
    integer,            intent(in)    :: lunit
    character(len=*),   intent(in)    :: name
    !
    if(lunit == loglun) then
      write(*,'(a,i3)')'logical unit must differ from',loglun
      return
    endif
    !
    if(lunit == 6) return
    !
    open(unit=lunit, file=trim(name)//'.txt', status='REPLACE')
    !
  end subroutine open_file
  !
  subroutine close_file()
    !
    implicit none
    !
    close(unit=loglun,status='KEEP')
    !
    if(lunit == 6) return
    !
    close(unit=lunit,status='KEEP')
    return
    !
  end subroutine close_file
  !
  subroutine log_msg(srname,infon,msg)
    !
    ! routine to write messages and error indicators onto log file
    !
    implicit none
    !
    character(len=*), intent(in) :: srname      ! subroutine name
    integer,          intent(in) :: infon       ! info number
    character(len=*), intent(in) :: msg         ! the message
    !
    write(unit= loglun, fmt='(a,i3,a)')trim(srname)// ' :: &
                                    &info= ',infon,' '//trim(msg)
    !
  end subroutine log_msg
  !
end module file_sup
!
!-------------------------------------------------------------------
!
module lin_lib
  !
  use def_kinds, ONLY : RP => RPX
  !
  interface lus
    module procedure lus1, lusk
  end interface
  !
contains
  !
  function ddist(x,y)
    !
    use def_kinds, ONLY : RP => RPX
    implicit none
    intrinsic abs, sqrt, size
    !
    real(kind=RP)                            :: ddist
    real(kind=RP), dimension(0:), intent(in) :: x, y
    optional                                 :: y
    !
    ! local variables and parameters
    !
    integer       :: j, kvec, level, kxy
    real(kind=RP) :: diff, dx, dy, hitest, sum, tmp, trm, xmax
    !
    real(kind=RP), parameter :: zer=0.0_RP, one=1.0_RP
    real(kind=RP), parameter :: cutlo = 8.232e-11_RP, &
                                cuthi = 1.304e19_RP
    !----------------------------------------------------------
    ! Computes either the euclidean distance between two 
    ! n-dimensional vectors x and y or the euclidean norm 
    ! of one such vector x. 
    !
    ! For the euclidean norm of a vector x use call ddist(x)
    ! for the euclidean distance between x and y use call(x,y)
    !
    ! If size(x) /= size (y) then zero is returned.
    !
    ! The algorithm follows the four-phase method of
    ! C. L. Lawson in the lapack routine dnrm2.f. As in
    ! dnrm2.f two built-in constants are used that are 
    ! hopefully applicable to all machines:
    !   cutlo = maximum of sqrt(u/eps) over known machines
    !   cuthi = minimum of sqrt(v)     over known machines
    ! where
    !   eps = smallest number such that 1.0 + eps .gt. 1.0
    !   u   = smallest positive number  (underflow limit)
    !   v   = largest  number           (overflow  limit)
    !
    ! Values for cutlo and cuthi listed in dnrm2.f are:
    !
    !   cutlo, s.p. u/eps = 2**(-102) for honeywell. close 
    !              seconds are univac and dec at 2**(-103)
    !              thus cutlo = 2**(-51) = 4.44089e-16
    !   cuthi, s.p. v = 2**127 for univac, honeywell, and dec
    !              thus cuthi = 2**(63.5) = 1.30438e19
    !   cutlo, d.p.  u/eps = 2**(-67) for honeywell and dec.
    !              thus cutlo = 2**(-33.5) = 8.23181d-11
    !   cuthi, d.p.  same as s.p.  cuthi = 1.30438d19
    !
    ! Thus two choices are
    !   cutlo := 8.232d-11,  cuthi := 1.304d19
    !   cutlo := 4.441e-16,  cuthi := 1.304e19
    !
    ! (NOTE: For the Pentium II the numbers are
    !   rad =  2
    !   u =  4.450147717014403E-308
    !   eps =  2.220446049250313E-16
    !   v =  1.797693134862316E+308
    !   cutlo =  1.415686533102923E-146
    !   cuthi =  1.340780792994260E+154)
    !
    ! In line with the four phases of dnrm2.f, the algorithm 
    ! uses four states identified by level = 0,1,2,3, 
    ! respectively, which correspond to the following cases:
    !
    !   level = 0  only zero terms have been found so far
    !   level = 1  all nonzero terms encountered so far do not
    !              exceed cutlo in modulus
    !   level = 2  there are some terms that are larger than
    !              cutlo in modulus but none exceeds 
    !              hitest = cutlo/real(n) 
    !   level = 3  there are terms that exceed hitest in modulus.
    !
    !  all state transitions can only increase the level.
    !------------------------------------------------------------
    !
    kxy = size(x)-1                  ! dimension index
    if(kxy < 0) then                 ! dimension error
      ddist = zer
      return
    endif
    !
    kvec = 1
    if(present(y)) then
      kvec = 2
      if(kxy+1 /= size(y)) then       ! vectors of different size
        ddist = zer
        return
      endif
    endif 
    !
    hitest = cuthi/real(kxy+1,kind=RP)
    if(hitest < cutlo) hitest = cutlo
    !
    sum = zer
    xmax = zer
    level = 0
    !
    do j = 0,kxy                      ! main loop over all terms
      dx = x(j)
      if(kvec == 1) then
        diff = dx                     ! one vector
      else
        dy = y(j)                     ! two vectors
        if(sign(one, dx) /= sign(one, dy)) then
          diff = dx - dy
        else
          if(abs(dx) < abs(dy)) then
            tmp = dy
            dy  = dx
            dx  = tmp
          endif
          if(dx == zer) then
            diff = zer
          else
            diff = dx*(one - dy/dx)
          endif
        endif
      endif
      !
      ! summation of the squares of the nonzero terms
      !
      if(diff /= zer) then
        trm = abs(diff) 
        if(trm <= cutlo) then
          !
          ! small nonzero terms -- transition to level 1
          !
          select case(level)
            !
            case(0)
              !
              level = 1
              xmax = trm
              tmp = trm/xmax
              sum = sum + tmp*tmp
            case(1)
              !
              if(trm > xmax ) then
                tmp = xmax/trm
                sum = one + sum*tmp*tmp
                xmax = trm
              else
                tmp = trm/xmax
                sum = sum + tmp*tmp
              endif
            case default
              ! 
              sum = sum + trm*trm
          end select
          !
        else
          !
          ! mid-sized terms -- transition to level 2
          !
          select case(level)
            !
            case(0)
              !
              level = 2
            case(1)
              !
              level = 2
              sum = (sum * xmax) * xmax
          end select        
          !
          if(trm <= hitest) then
            sum = sum + trm*trm
          else
            !
            ! large terms
            !
            if(level == 2) then
              !
              ! transition to level = 3
              !
              level = 3
              sum = (sum / trm) / trm
              xmax = trm
            endif
            if(trm > xmax) then
              tmp = xmax/trm
              sum = one + sum*tmp*tmp
              xmax = trm
            else
              tmp = trm/xmax
              sum = sum + tmp*tmp
            endif
          endif
        endif
      endif
    enddo                                ! end of main loop
    !
    if((xmax == zer) .or. (level == 2)) then
      ddist = sqrt(sum)
    else
      ddist = xmax*sqrt(sum)
    endif
    !
  end function ddist
  !
  subroutine housg(x,tau)
    !
    use def_kinds, ONLY : RP => RPX
    implicit none
    intrinsic abs,sqrt,size,tiny
    !
    real(kind=RP), dimension(0:), intent(inout) :: x
    real(kind=RP),                intent(out)   :: tau
    !
    ! local variables and parameters
    !
    integer       :: i,kn,kred
    real(kind=RP) :: a1,a2,alpha,beta,rsafmin,tmp,xnrm,xscal
    !
    real(kind=RP), parameter :: zer=0.0_RP, one=1.0_RP
    real(kind=RP), parameter :: safmin=tiny(one)
    !--------------------------------------------------------
    ! For a given vector x stored in the array x(0:kn), an
    ! n=kn+1 dimensional Householder reflector
    !
    !      H = I - tau * y^T y,     H^T * H = I
    ! 
    ! is generated where tau is a scalar and y is an
    ! n-dimensional vector stored in the array y(0:kn)
    ! such that y(0) = 1.0 and
    !
    !      H * x = ( beta ),  beta = scalar
    !              (   0  )
    !
    ! Because of H^T * H = I we have x^T x = beta^2 and
    ! thus
    !
    !   beta = sqrt(x^T x),     H^T (beta) = x
    !                               ( 0  ) 
    ! whence
    !     tau = (beta - x(0))/beta
    !     y(0) = 1.0
    !     y(1:kn) = xscal*x(1:kn), xscal = 1/(beta - x(1))
    !
    ! If x = 0, then tau = 0 and H = I, else  1 <= tau <= 2.
    !
    ! The scalar beta is returned in x(0), and x(1:n) will
    ! contain y(1:n).
    !
    ! This is an edited version of the lapack routine DLARFG
    !------------------------------------------------------------
    kn = size(x)-1                        ! dimension index
    !
    xnrm = zer
    if(kn > 0) xnrm = ddist(x(1:kn))
    !
    if(kn <= 0 .or. xnrm == zer) then
      tau = zer                               ! quick return
      return
    endif
    !
    alpha = x(0)
    a1 = xnrm
    a2 = abs(alpha)
    if(a1 < a2) then
      a1 = a2
      a2 = xnrm
    endif
    if(a2 == zer) then
      beta = a1
    else
      tmp = a1/a2
      beta = a2*sqrt(one + tmp*tmp)
    endif
    if(alpha > zer) beta = -beta
    if(abs(beta) >= safmin) then
      tau = (beta - alpha)/beta
      xscal = one/(alpha-beta)
      alpha = beta
      !
    else
      !
      ! xnrm, beta may be inaccurate, scale x and recompute
      !
      rsafmin = one/safmin
      kred = 0
      do while(abs(beta) < safmin)
        !
        kred = kred + 1
        x(1:kn) = rsafmin*x(1:kn)
        beta    = rsafmin*beta
        alpha   = rsafmin*alpha
        !
      enddo
      !
      ! The new beta satisfies safmin <= beta <= 1.0
      !
      xnrm = ddist(x(1:kn))
      a1 = xnrm
      a2 = abs(alpha)
      if(a1 < a2) then
        a1 = a2
        a2 = xnrm
      endif
      if(a2 == zer) then
        beta = a1
      else
        tmp = a1/a2
        beta = a2*sqrt(one + tmp*tmp)
      endif
      if(alpha > zer) beta = -beta
      tau   = (beta - alpha)/beta
      xscal = one/(alpha-beta)
      alpha = beta
      do i = 1,kred
        alpha = alpha*safmin
      enddo
    endif
    x(0) = alpha
    x(1:kn) = xscal*x(1:kn)
    !
  end subroutine housg
  !
  subroutine housl(a,x,tau,info)
    !
    use def_kinds, ONLY : RP => RPX
    implicit none
    intrinsic dot_product,size
    !
    real(kind=RP),                   intent(in)    :: tau
    real(kind=RP), dimension(0:),    intent(in)    :: x
    real(kind=RP), dimension(0:,0:), intent(inout) :: a
    integer,                         intent(out)   :: info
    !
    integer                  :: km,kn,j    ! local variables
    real(kind=RP)            :: sum
    real(kind=RP), parameter :: zer=0.0_RP
    !--------------------------------------------------------
    ! Multiplies a given m x n matrix A from the left by an 
    ! m-dimensional Householder reflector to obtain
    !
    !  A := (I - tau*x*x^T)*A ,  dim x = m
    !
    ! A is stored in the array a(0:km,0:kn)
    ! x is stored in the array x(0:km) but only the 
    ! components x(1),...,x(km) are used and instead
    ! of x(0) the value 1.0 is enforced.
    !--------------------------------------------------------
    km = size(x)-1
    kn = size(a,2)-1
    !
    if(km+1 /= size(a,1)) then                ! dimension error
      info = 1
      return
    endif
    !
    info = 0
    if(km < 0 .or. kn < 0 .or. tau == zer) return ! zero data
    !
    ! form w := A^T*x and  A := A - tau*x*w^T
    !
    do j = 0,kn,1
      sum = a(0,j)                     ! x(0) = 1 is enforced
      if(km > 0) sum = sum + dot_product(a(1:km,j),x(1:km))
      if(sum /= zer) then
        sum = -tau*sum
        a(0,j) = a(0,j) + sum
        if(km > 0) a(1:km,j) = a(1:km,j) + sum*x(1:km)
      endif
    enddo
    !
  end subroutine housl
  !
  subroutine housr(a,x,tau,info)
    !
    use def_kinds, ONLY : RP => RPX
    implicit none
    intrinsic dot_product,size
    !
    real(kind=RP),                   intent(in)    :: tau
    real(kind=RP), dimension(0:),    intent(in)    :: x
    real(kind=RP), dimension(0:,0:), intent(inout) :: a
    integer,                         intent(out)   :: info
    !
    integer       :: km,kn,j            ! local variables
    real(kind=RP) :: sum
    real(kind=RP), parameter :: zer=0.0_RP
    !--------------------------------------------------------
    ! Multiplies a given m x n matrix A from the right by an 
    ! m-dimensional Householder reflector to obtain
    !
    !    A := A*(I - tau*x*x^T) ,  dim x = n 
    !
    ! A is stored in the array a(0:km,0:kn)
    ! x is stored in the array x(0:km) but only the 
    ! components x(1),...,x(km) are used and instead
    ! of x(0) the value 1.0 is enforced.
    !--------------------------------------------------------
    km = size(a,1)-1
    kn = size(x)-1
    !
    if(kn+1 /= size(a,2)) then             ! dimension error
      info = 1
      return
    endif
    !
    info = 0
    if(km < 0 .or. kn < 0 .or. tau == zer) return  !zero data
    !
    ! form w := A^T*x and  A := A - tau*x*w^T
    !
    do j = 0,km,1
      sum = a(j,0)          ! note that x(0) = 1 is enforced
      if(kn > 0) sum = sum + dot_product(a(j,1:kn),x(1:kn)) 
      if(sum /= zer) then
        sum = -tau*sum
        a(j,0) = a(j,0) + sum
        if(kn > 0) a(j,1:kn) = a(j,1:kn) + sum*x(1:kn)
      endif
    enddo
    !
  end subroutine housr
  !
  subroutine lqas(a,ipiv,tau,y,x,info)
    !
    use def_kinds, ONLY : RP => RPX
    implicit none
    intrinsic dot_product, size
    !
    real(kind=RP), dimension(0:,0:), intent(in)  :: a
    integer, dimension(0:),          intent(in)  :: ipiv
    real(kind=RP), dimension(0:),    intent(in)  :: tau,y
    integer,                         intent(out) :: info
    real(kind=RP), dimension(0:),    intent(out) :: x
    !
    integer       :: i,ip,j,km,kn     ! local variables
    real(kind=RP) :: sum
    real(kind=RP), parameter ::  zer=0.0_RP
    !----------------------------------------------------------------
    ! For an  m x n matrix A with m <= n and rank A = m and a given
    ! given m-vector y, compute the unique solution x of Ax = y
    ! which is orthogonal to ker A. The algorithm 
    !
    !     1.  solve Ly := Py
    !     2.  x := (y)
    !              (0)
    !     2.  x := Q^t x
    !
    ! is used under the assumption that the arrays a, tau, and ipiv 
    ! contain the LQ-factorization A = P^T L Q of the matrix A.
    !
    ! The return indicator has the values
    !     info < 0  zero pivot: info = -j-1 signifies that the
    !               j-th diagonal element of the triangular matrix 
    !               r is zero.
    !----------------------------------------------------------------
    km = size(a,1)-1
    kn = size(a,2)-1
    !
    info = 0
    if((km < 0) .or. (kn < 0) .or. km > kn) return  ! dimension error
    !
    if(kn == 0) then
      if(a(0,0) == zer) then                    ! scalar case
        info = -1                               ! zero pivot
      else
        x(0) = y(0)/a(0,0)
      endif
      return
    endif
    !
    ! solve Lz = Py for z and extend to an m-vector
    !
    do i = 0,kn
      if(i <= km) then
        if(a(i,i) == zer) then
          info = -i-1                           ! zero pivot
          return
        endif
        ip = ipiv(i)
        sum = y(ip)
        if(i > 0) then
          sum = sum - dot_product(a(i,0:i-1),x(0:i-1))
        endif
        x(i) = sum/a(i,i)
      else
        x(i) = zer
      endif
    enddo
    !
    ! form q^t*x
    !
    do j = km,0,-1
      sum = x(j)
      if(j < kn) sum = sum + dot_product(a(j,j+1:kn),x(j+1:kn))
      if(sum /= zer) then
        sum = - tau(j)*sum
        x(j) = x(j) + sum
        if(j < kn) x(j+1:kn) = x(j+1:kn) + sum*a(j,j+1:kn)
      endif
    enddo
    !
  end subroutine lqas
  !
  subroutine lqf(a,ipiv,tau,cnrm,info)
    !
    use def_kinds, ONLY : RP => RPX
    implicit none
    intrinsic abs,sqrt,size,tiny
    !
    real(kind=RP), dimension(0:,0:), intent(inout) :: a
    integer,       dimension(0:),    intent(out)   :: ipiv
    real(kind=RP), dimension(0:),    intent(out)   :: tau,cnrm
    integer,                         intent(out)   :: info
    !
    ! local variables and parameters
    !
    integer                  :: i,imax,itmp,j,km,kn
    real(kind=RP)            :: cnr,cnrj,taui,tmp1,tmp2,tmp3
    real(kind=RP), parameter :: fact=0.05_RP,zer=0.0_RP,one=1.0_RP
    !-------------------------------------------------------------
    ! Computes the LQ factorization  PA = LQ  with column 
    ! pivoting  on the columns of an m by n matrix A under the 
    ! assumption that  m <= n.
    ! 
    ! The matrix Q is represented as a product of m elementary 
    ! reflectors Q = H(0)*H(1)*...*H(km), H(i) = I - tau*v*v^T, 
    ! where tau is a scalar, and v = v(0:n-1) a vector with
    ! v(i) = 1 and v(0:i-1) = 0 for i > 0.
    !
    ! Upon return, the components v(i+1:n-1) are stored in 
    ! a(i,i+1:n-1) and tau in the array tau(0:m-1). The matrix 
    ! P is stored as follows: If jpvt(j) = i then the jth row
    ! of p is the (i+1)-st canonical unit vector.
    !
    ! The return indicator has the values 
    !     info = 0 :  successful exit
    !     info = 1 :  dimension error
    !     info < 0 :  for info = -i the i-th pivot is zero
    !-------------------------------------------------------------
    km = size(a,1)-1
    kn = size(a,2)-1
    !
    if(kn < km) then
      info = 1                           ! dimension error
      return
    endif
    !
    info = 0
    if(km < 0 .or. kn < 0) return        ! zero data
    !
    do i = 0,km,1
      cnr = ddist(a(i,0:kn))             ! Initialize row norms    
      tau(i)  = cnr
      cnrm(i)  = cnr
      ipiv(i) = i
    enddo
    !
    do i = 0,km,1                        ! Main factorization loop
      !
      cnr = tau(i)
      imax = i
      if(i < km) then
        do j=i+1,km,1                    ! determine pivot row
          if(tau(j) > cnr) then
            imax = j
            cnr = tau(j)
          endif
        enddo
        !
        if(imax /= i) then               ! swap rows if needed
          do j = 0,kn,1
            tmp1 = a(imax,j)
            a(imax,j) = a(i,j)
            a(i,j) = tmp1
          enddo
          tau(imax) = tau(i)
          cnrm(imax) = cnrm(i)
          itmp = ipiv(imax)
          ipiv(imax) = ipiv(i)
          ipiv(i) = itmp
        endif
      endif
      if(cnr <= zer) then
        info = -i-1                       ! zero pivot
      else
        !
        taui = zer
        if(i < kn) call housg(a(i,i:kn),taui)
                                         ! generate reflector h(i)
        if(i < km) then
          call housr(a(i+1:km,i:kn),a(i,i:kn),taui,info)
                   ! apply h(i) to a(i+1:km,i:kn) from the right
          !
          do j = i+1,km,1                ! loop to update row norms
            cnrj = tau(j)
            if(cnrj /= zer) then
              tmp1 = abs(a(j,i))/cnrj
              tmp2 = one - tmp1*tmp1
              if(tmp2 <= zer) then
                tmp2 = zer
                tmp3 = zer
              else
                tmp3 = sqrt(tmp2)
              endif
              tmp1 = cnrj/cnrm(j)
              tmp2 = one + fact*tmp2*tmp1*tmp1
              if(tmp2 == one) then
                tau(j) = ddist(a(j,i+1:kn))
                cnrm(j) = tau(j)
              else
                tau(j) = tau(j)*tmp3
              endif
            endif
          enddo                        ! end row-norm update
        endif
        tau(i) = taui
      endif
    enddo                              ! end main loop
    !
  end subroutine lqf
  !
  subroutine luf(a,ipiv,info)
    !
    use def_kinds, ONLY : RP => RPX
    implicit none
    intrinsic abs, size
    !
    real(kind=RP), dimension(0:,0:), intent(inout) :: a
    integer,       dimension(0:),    intent(out)   :: ipiv
    integer,                         intent(out)   :: info
    !
    integer       :: j,k,kpiv,km,kn        ! local variables
    real(kind=RP) :: apiv, tmp
    real(kind=RP), parameter :: zer=0.0_RP
    !-------------------------------------------------------------
    ! Computes an LU factorization of a general n x n matrix A
    ! using partial pivoting with row interchanges:
    !
    !   A = P * L * U
    !
    !   P = permutation matrix, 
    !   L = lower triangular matrix with unit diagonal elements
    !   U = upper triangular matrix
    !
    ! The factors L and U are returned in the corresponding parts 
    ! of A, except that the diagonal of L is not stored. The 
    ! permutation P is returned in the pivot array ipiv. More 
    ! specifically, row a(i,:) of the matrix was interchanged 
    ! with row a(ipiv(i),:).
    !
    ! The values of the return indicator are
    !      info = 0 :  successful exit
    !      info = 1 : data error
    !      info < 0 :  if info = i, then u(i,i) is exactly zero 
    !                  but the factorization has been completed.
    !-------------------------------------------------------------
    km = size(a,1)-1
    kn = size(a,2)-1
    !
    info = 0
    if(kn < 0 .or. km < 0) then
      info = 1
      return                                ! data error
    endif
    !
    do k = 0,kn,1                           ! main loop over columns
      !
      kpiv = k                              ! determine pivot 
      apiv = abs(a(k,k))
      if(k < kn) then      
        do j = k+1,kn,1
          if(abs(a(j,k)) > apiv) then
            kpiv = j
            apiv = abs(a(j,k))
          endif
        enddo
      endif
      ipiv(k) = kpiv
      if(apiv == zer) then
        info = -k-1                         ! zero pivot
      else
        !
        if(kpiv /= k) then                  ! apply interchange    
          do j = 0,kn,1
            tmp = a(kpiv,j)
            a(kpiv,j) = a(k,j)
            a(k,j) = tmp
          enddo
        endif
        !
        if(k < kn) then
          a(k+1:kn,k) = a(k+1:kn,k)/a(k,k)    ! compute multipliers 
          !
          do j = k+1,kn                ! update trailing submatrix
            a(k+1:kn,j) = a(k+1:kn,j) - a(k,j)*a(k+1:kn,k)
          enddo
        endif
      endif
    enddo                                   ! end main loop
    !
  end subroutine luf
  !
  subroutine lus1(a,ipiv,x,info)
    !
    use def_kinds, ONLY : RP => RPX
    implicit none
    intrinsic size
    !
    real(kind=RP), dimension(0:,0:), intent(in)    :: a
    integer,       dimension(0:),    intent(in)    :: ipiv
    real(kind=RP), dimension(0:),    intent(inout) :: x
    integer,                         intent(out)   :: info
    !
    integer       :: i,ip,j,kn               ! local variables
    real(kind=RP) :: tmp
    real(kind=RP), parameter :: zer=0.0_RP
    !----------------------------------------------------------------
    ! Solves a system of linear equations A * z = x for an n x n 
    ! matrix A and an n-vector x and returns the result in x. It is 
    ! assumed that the arrays a and ipiv contain the LU factorization 
    !  returned by subroutine luf.
    !
    ! The values of the return indicator are
    !      info = 0:  successful exit
    !      info = 1:  input-data error
    !      info < 0:  if iret = -k, then the k-th pivot is zero.
    !-------------------------------------------------------------
    kn = size(a,1)-1
    !
    if(kn+1 /= size(a,2) .or. kn+1 /= size(x)) then
      info = 1                               ! dimension error
      return
    endif
    !
    info = 0
    if(kn == 0) then                         ! scalar case
      if(a(0,0) == zer) then
        info = -1                            ! zero pivot
      else         
        x(0) = x(0)/a(0,0)
      endif
      return
    endif
    !
    do i = 0,kn,1
      ip = ipiv(i)
      if(ip /= i)then
        tmp  = x(ip)            ! apply interchanges to right side.
        x(ip) = x(i)
        x(i) = tmp
      endif
    enddo
    !
    ! solve l*z = x, overwriting b with z.
    !
    do j = 0,kn-1
      tmp = x(j)
      if(tmp /= zer .and. j < kn) &
                x(j+1:kn) = x(j+1:kn) - tmp*a(j+1:kn,j)
    enddo
    !
    ! solve u*z = x, overwriting x with z.
    !
    do j = kn,0,-1
      if(x(j) /= zer) then
        if(a(j,j) == zer) then
          info = -j-1                      ! zero pivot
        else         
          x(j) = x(j)/a(j,j)
          if(j > 0) x(0:j-1) = x(0:j-1) - x(j)*a(0:j-1,j) 
        endif
      endif
    enddo
    !
  end subroutine lus1
  !
  subroutine lusk(a,ipiv,b,info)
    !
    use def_kinds, ONLY : RP => RPX
    implicit none
    intrinsic size
    !
    real(kind=RP), dimension(0:,0:), intent(in)    :: a
    integer,       dimension(0:),    intent(in)    :: ipiv
    real(kind=RP), dimension(0:,0:), intent(inout) :: b
    integer,                         intent(out)   :: info
    !
    integer :: i,ip,j,k1,k2                  ! local variables
    real(kind=RP) :: tmp
    real(kind=RP), parameter ::  zer=0.0_RP
    !----------------------------------------------------------------
    ! Solves a system of linear equations a * z = b for an n x n 
    ! matrix a and an n x k matrix b on the right side. It is 
    ! assumed that the arrays a and ipiv contain the LU
    ! factorization returned by subroutine luf.
    !
    ! The solutions are returned in the array b.
    !
    ! The values of the return indicator are
    !      info = 0 :  successful exit
    !      info = 1 :  input-data error
    !      info < 0 :  if info = -j, then the j-th pivot is zero.
    !-------------------------------------------------------------
    k1 = size(a,1)-1
    k2 = size(b,2)-1
    !
    if(k1+1 /= size(a,2) .or. k1+1 /= size(b,1)) then
      info = 1                               ! dimension error
      return
    endif
    !
    info = 0
    if(k1 == 0) then                         ! scalar case
      if(a(0,0) == zer) then
        info = -1                                ! zero pivot
      else
        b(0,0:k2) = b(0,0:k2)/a(0,0)
      endif
      return
    endif
    !
    do i = 0,k1,1            ! apply row interchanges to right side.
      ip = ipiv(i)
      if(ip /= i) then
      do j = 0,k2,1
        tmp  = b(ip,j)
        b(ip,j) = b(i,j)
        b(i,j) = tmp
      enddo
      endif
    enddo
    !
    ! solve l*x = b, overwriting b with x.
    !
    do j = 0,k2,1
      do i = 0,k1-1
        if(b(i,j) /= zer .and. i < k1) &
          b(i+1:k1,j) = b(i+1:k1,j) - b(i,j)*a(i+1:k1,i)
      enddo
    enddo
    !
    ! solve u*x = b, overwriting b with the solution.
    !
    do j = 0,k2,1
      do i = k1,0,-1
        if(b(i,j) /= zer ) then
          if(a(i,i) == zer ) then
            info = -i-1                   ! zero pivot
            return
          else
            b(i,j) = b(i,j)/a(i,i)
          endif
          if(b(i,j) /= zer .and. i > 0) &
             b(0:i-1,j) = b(0:i-1,j) - b(i,j)*a(0:i-1,i)
        endif
      enddo
    enddo
    !
  end subroutine lusk
  !
  subroutine qrf(a,ipiv,tau,cnrm,info)
    !
    use def_kinds, ONLY : RP => RPX
    implicit none
    intrinsic abs, sqrt, size,tiny
    !
    real(kind=RP), dimension(0:,0:), intent(inout) :: a
    integer,       dimension(0:),    intent(out)   :: ipiv
    real(kind=RP), dimension(0:),    intent(out)   :: tau, cnrm
    integer,                         intent(out)   :: info
    !
    ! local variables and parameters
    !
    integer       :: i,imax,itmp,j,km,kn,kmn
    real(kind=RP) :: cnr,cnrj,taui,tmp1,tmp2,tmp3
    !
    real(kind=RP), parameter :: fact=0.05_RP,zer=0.0_RP,one=1.0_RP
    !-------------------------------------------------------------
    ! Computes the QR factorization  a*p = q*r  with pivoting
    ! on the columns of an m by n matrix a.
    ! 
    ! The matrix q is represented as a product of m elementary 
    ! reflectors q = h(1)*h(2)*...*h(k), h(i) = i - tau * v * v^T, 
    ! where k = min(m,n), tau is a real scalar, and v is a real 
    ! vector with v(1:i-1) = 0 and v(i) = 1.
    !
    ! Upon return, the components v(i+1:n) are stored in 
    ! a(i,i+1:n) and tau in the array tau(i). The matrix p is 
    ! stored as follows: If ipiv(j) = i then the jth column of 
    ! p is the i-th canonical unit vector.
    !
    ! The return indicator has the values 
    !     info = 0 :  successful exit
    !     info < 0 :  if info = -i then the i-th pivot
    !                 was found to be zero
    !-------------------------------------------------------------
    km = size(a,1)-1
    kn = size(a,2)-1
    kmn = min(km,kn)
    !
    info = 0
    !
    do i = 0,kn,1
      cnr = ddist(a(0:km,i))     
      tau(i)  = cnr                     ! Initialize row norms
      cnrm(i)  = cnr
      ipiv(i) = i                       ! Initialize pivot array
    enddo
    !
    do i = 0,kmn,1                      ! Main factorization loop
      !
      cnr = tau(i)
      imax = i
      if(i < kmn) then
        do j=i+1,kn,1                   ! determine pivot column
          if(tau(j) > cnr) then
            imax = j
            cnrm = tau(j)
          endif
        enddo          
      endif
      if(imax /= i) then                ! swap columns if needed
        do j = 1,km
          tmp1 = a(j,imax)
          a(j,imax) = a(j,i)
          a(j,i) = tmp1
        enddo
        itmp = ipiv(imax)
        ipiv(imax) = ipiv(i)
        ipiv(i) = itmp
        tau(imax) = tau(i)
        cnrm(imax)  = cnrm(i)
      endif                             ! end swap
      !
      if(cnr == zer) then
        info = -i-1                     ! zero pivot
        !
      else
        !
        taui = zer
        if(i < km) call housg(a(i:km,i),taui)
                             ! generate elementary reflector h(i)
        if(i < kn) then
          call housl(a(i:km,i+1:kn),a(i:km,i),taui,info)
                    ! apply h(i) to a(i:km,i+1:kn) from the left
          !
          do j = i+1,kn                    ! update column norms
            cnrj = tau(j)
            if(cnrj /= zer) then
              tmp1 = abs(a(i,j))/cnrj
              tmp2 = one - tmp1*tmp1
              if(tmp2 <= zer) then
                tmp2 = zer
                tmp3 = zer
              else
                tmp3 = sqrt(tmp2)
              endif
              tmp1 = cnrj/cnrm(j)
              tmp2 = one + fact*tmp2*tmp1*tmp1
              if(tmp2 == one) then
                tau(j) = ddist(a(i+1:km,j))
                cnrm(j) = tau(j)
              else
                tau(j) = tau(j)*tmp3
              endif
            endif
          enddo                           ! end norm update
        endif
        tau(i) = taui
      endif
      !
    enddo                                 ! end main loop
    !
  end subroutine qrf
  !
  subroutine qrs(a,tau,y,x,info)
    !
    use def_kinds, ONLY : RP => RPX
    implicit none
    intrinsic dot_product,size
    !
    real(kind=RP), dimension(0:,0:), intent(in)    :: a
    real(kind=RP), dimension(0:),    intent(in)    :: tau
    real(kind=RP), dimension(0:),    intent(inout) :: y
    real(kind=RP), dimension(0:),    intent(out)   :: x
    integer,                         intent(out)   :: info
    !
    integer       :: j,km,kn             ! local variables
    real(kind=RP) :: sum, tmp
    real(kind=RP), parameter :: zer=0.0_RP
    !----------------------------------------------------------------
    ! For an m x n matrix a with m >= n and rank a = n, and a given
    ! m-vector y, compute the least squares solution
    !
    !      min { || a*x - y || ; x in R^n }
    !
    ! by the algorithm  
    !
    !      1.  y := q^t y
    !      2.  solve r*x := (q^T y)(1:n)
    !
    ! under the assumption that the arrays a and tau contain the
    ! QR-factorization of an m x n matrix a.
    !
    ! The routine returns q^T y in the array y and the least squares
    ! solution in the array x. The array x may be identified with y 
    ! if q^T y is not needed.
    !
    ! The return indicator has the values
    !     info = 0 successful return
    !     info = 1 dimension error
    !     info < 0 if info = -i then the i-th element of the
    !              triangular matrix r is zero
    !----------------------------------------------------------------
    km = size(a,1)-1
    kn = size(x)-1
    !
    if(kn+1 /= size(a,2) .or. km+1 /= size(y)) then
      info = 1                               ! dimension error
      return
    endif
    !
    info = 0
    if(km == 0) then
      if(a(0,0) == zer) then                 ! scalar case
        info = -1                            ! zero pivot
      else
        x(0) = y(0)/a(0,0)
      endif
      return
    endif
    !
    ! compute q^t*y
    !
    do j = 0,kn
      sum = y(j)
      if(j < km) sum = sum + dot_product(a(j+1:km,j),y(j+1:km))
      if(sum /= zer) then
        tmp = - tau(j)*sum
        y(j) = y(j) + tmp
        if(j < km) y(j+1:km) = y(j+1:km) + tmp*a(j+1:km,j)
      endif
    enddo
    !
    x(0:kn) = y(0:kn)                         ! copy y into x
    !
    ! solve r*x := y(1:n)
    !
    do j = kn,0,-1
      if(a(j,j) == zer) then
        info = -j-1                           ! zero pivot
        return
      endif
      x(j) = x(j)/a(j,j)
      if(j > 0) &
        x(1:j-1) = x(1:j-1) - x(j)*a(0:j-1,j)
    enddo
    !
  end subroutine qrs
  !
end module lin_lib
