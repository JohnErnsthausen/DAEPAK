DAEPAK
======

DAE solver based on geometric techniques

W.C. Rheinboldt wrote the original solvers and supporting libraries in Fortran 95 released on 11/06/2000.
The initial commit is exactly these files unpacked into a directory structure for development.

No one is actively maintaining DAEPAK. I value this technology and created this git repository
to bring awareness, for continue maintainence, and for free distribution of this DAE solver library.


Differential Algebraic Equation Solvers
=======================================

Solver for quasilinear DAEs of the following form

### DAEQ1
    A(t,u)u' = H(t,u)
    F(t,u)   = 0
    u(t0) = u0, such that F(t0,u0) = 0
    
    dim(rge A)= dim(rge H) = nd, dim(rge F) = na, 
    nu = nd+na, nw = 0, dim(u) = nu,
    
    requires subroutine DAEFCT with 
      fname = 'amt ' for evaluating A(t,u)
      fname = 'hf  ' for evaluating H(t,u)
      fname = 'ff  ' for evaluating F(t,u)
      fname = 'dff ' for evaluating DF(t,u)

### DAEN1
    u' = v
    F(t,u,v,w) = 0,
    u(t0) = u0, v(t0) = up0, w(t0) = w0, 
                such that F(t0,u0,up0,w0) = 0
    
    dim(u) = nu, dim(w) = nw,
    nd = nu, na = nu + nw, dim(rge F)= na
    
    requires subroutine DAEFCT with 
      fname = 'ff  ' for evaluating F(t,u)
      fname = 'dff ' for evaluating DF(t,u)

### DAEQ2
    A(t,u)u' + B(t,u)w = H(t,u)
    F(t,u)             = 0
    u(t0) = u0, such that F(t0,u0) = 0
    
    dim(rge A)= dim(rge H) = nd, dim(rge F)= na, 
    dim(w) = nw, dim(u) = nu, nu + nw = nd + na
    
    requires subroutine DAEFCT with 
      fname = 'amt ' for evaluating A(t,u)
      fname = 'bmt ' for evaluating B(t,u)
      fname = 'hf  ' for evaluating H(t,u)
      fname = 'ff  ' for evaluating F(t,u)
      fname = 'dff ' for evaluating DF(t,u)

### NOHOL
    A(u) u" + d_{u'}F(t,u,u')^T w = H(t,u,u')
    F(t,u,u')                     = 0
    u(t0) = u0, u'(t0) = up0
                such that Ft0,u0,up0) = 0
    
    dim(rge A)= dim(rge H) = nd, dim(rge F)= na, 
    dim(w) = na, dim(u) = nd,
    
    requires subroutine DAEFCT with 
      fname = 'amt ' for evaluating A(u)
      fname = 'hf  ' for evaluating H(t,u,u')
      fname = 'ff  ' for evaluating F(t,u,u')
      fname = 'dff ' for evaluating DF(t,u,u')

### DAEQ3
    A(t,u,u')u" + B(t,u,u')w = H(t,u,u')
    F(t,u)                   = 0
    u(t0)=u0, u'(t0)=up0, such that F(t0,u0) = 0
    
    dim(rge F)= na, dim(w) = nw, 
    dim(rge A)= dim(rge H) = nd, dim(u) = nd+na-nw
    
    requires subroutine daefct with 
      fname = 'amt ' for evaluating A(t,u,u')
      fname = 'bmt ' for evaluating B(t,u,u')
      fname = 'hf  ' for evaluating H(t,u,u')
      fname = 'ff  ' for evaluating F(t,u)
      fname = 'dff ' for evaluating DF(t,u)
      fname = 'd2ff' for evaluating D2F(t,u)(v,v)

### EULAG
    A(u) u" + DF(u)^T w = H(u,u')
    F(u)                = 0
    u(t0) = u0, such that F(u0) = 0
    
    dim(rge A) = dim(rge H) = nd, dim(rge F) = na
    dim(w) = na, dim(u) = nd
    
    requires subroutine DAEFCT with 
      fname = 'amt ' for evaluating A(u)
      fname = 'hf  ' for evaluating H(u,u')
      fname = 'ff  ' for evaluating F(u)
      fname = 'dff ' for evaluating DF(u)
      fname = 'd2ff' for evaluating D2F(u)(v,v)

Getting Started with DAEPAK
===========================
>
> The main driver file has been renamed to *_prog.f95 to enable rake to pick it up as the main program.
>
> On CygWin version 1.7.14, the solvers partially work. All files compile under gfortran -std=f95. Examples
> ap_q2 and eight run but the example tskate fails due to step size and the examples ap_q3, catal, and
> sevbd fail on execution with a core dump. The linked list pattern in the module file_sup caused the program
> to core dump. A direct file open was adopted. The initial condition in tskate was corrected so that x(4)=1.
> The library now compiles and runs the examples without coredump and produces the runtime results Prof. Rheinboldt
> appened to the main example program.


Roadmap for the Future
======================
-  Incorporate unit testing
-  Move examples outside the library
-  Move linear algebra support to its own library
-  Move ordinary differential equation solver support to its own library
-  Move manifold support to its own library
-  Incorporate the examples from prior F77 release

