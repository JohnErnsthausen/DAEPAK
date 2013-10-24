DAEPAK-2.0.0
============

DAE solver based on geometric techniques

> W.C. Rheinboldt wrote the original solvers and supporting libraries in Fortran 95 which he released
> on 11/06/2000. The original tar files are available from various repositories on the internet. The
> initial commit is exactly these files unpacked into the standard directory structure for gem
> development. The main driver file has been renamed to *_prog.f95 to enable rake to pick is up. On
> CygWin version 1.7.14, the solvers partially work. All files compile under gfortran -std=f95. Examples
> ap_q2 and eight run but the example tskate fails due to step size and the examples ap_q3, catal, and
> sevbd fail on execution with a core dump.

The author believes no one is actively maintaining DAEPAK and has created this git repository to continue maintaining this DAE solver library.

Differential Algebraic Equation Solvers
=======================================

Solver for quasilinear DAEs of the following form

### DAEQ1
    a(t,u)u' = h(t,u)
    f(t,u)   = 0
    u(t0) = u0, such that f(t0,u0) = 0
    
    dim(rge a)= dim(rge h) = nd, dim(rge f) = na, 
    nu = nd na, nw = 0, dim(u) = nu,
    
    requires subroutine daefct with 
      fname = 'amt ' for evaluating a(t,u)
      fname = 'hf  ' for evaluating h(t,u)
      fname = 'ff  ' for evaluating f(t,u)
      fname = 'dff ' for evaluating df(t,u)

### DAEN1
    u' = v 
    f(t,u,v,w) = 0,
    u(t0) = u0, v(t0) = up0, w(t0) = w0, 
                such that f(t0,u0,up0,w0) = 0
    
    dim(u) = nu, dim(w) = nw,
    nd = nu, na = nu + nw, dim(rge f)= na
    
    requires subroutine daefct with 
      fname = 'ff  ' for evaluating f(t,u)
      fname = 'dff ' for evaluating df(t,u)

### DAEQ2
    a(t,u)u' + b(t,u)w = h(t,u)
    f(t,u)             = 0
    u(t0) = u0, such that f(t0,u0) = 0
    
    dim(rge a)= dim(rge h) = nd, dim(rge f)= na, 
    dim(w) = nw, dim(u) = nu, nu + nw = nd + na
    
    requires subroutine daefct with 
      fname = 'amt ' for evaluating a(t,u)
      fname = 'bmt ' for evaluating b(t,u)
      fname = 'hf  ' for evaluating h(t,u)
      fname = 'ff  ' for evaluating f(t,u)
      fname = 'dff ' for evaluating df(t,u)

### NOHOL
    a(u) u" + d_{u'}f(t,u,u')^T w = h(t,u,u')
    f(t,u,u')                     = 0
    u(t0) = u0, u'(t0) = up0
                such that f(t0,u0,up0) = 0
    
    dim(rge a)= dim(rge h) = nd, dim(rge f)= na, 
    dim(w) = na, dim(u) = nd,
    
    requires subroutine daefct with 
      fname = 'amt ' for evaluating a(u)
      fname = 'hf  ' for evaluating h(t,u,u')
      fname = 'ff  ' for evaluating f(t,u,u')
      fname = 'dff ' for evaluating df(t,u,u')

### DAEQ3
    a(t,u,u')u" + b(t,u,u')w = h(t,u,u')
    f(t,u)                   = 0
    u(t0)=u0, u'(t0)=up0, such that f(t0,u0) = 0
    
    dim(rge f)= na, dim(w) = nw, 
    dim(rge a)= dim(rge h) = nd, dim(u) = nd+na-nw
    
    requires subroutine daefct with 
      fname = 'amt ' for evaluating a(t,u,u')
      fname = 'bmt ' for evaluating b(t,u,u')
      fname = 'hf  ' for evaluating h(t,u,u')
      fname = 'ff  ' for evaluating f(t,u)
      fname = 'dff ' for evaluating df(t,u)
      fname = 'd2ff' for evaluating d2f(t,u)(v,v)

### EULAG
    a(u) u" + df(u)^T w = h(u,u')
    f(u)                = 0
    u(t0) = u0, such that f(u0) = 0
    
    dim(rge a) = dim(rge h) = nd, dim(rge f) = na
    dim(w) = na, dim(u) = nd
    
    requires subroutine daefct with 
      fname = 'amt ' for evaluating a(u)
      fname = 'hf  ' for evaluating h(u,u')
      fname = 'ff  ' for evaluating f(u)
      fname = 'dff ' for evaluating df(u)
      fname = 'd2ff' for evaluating d2f(u)(v,v)

Version:  W. Rheinboldt, Univ. of Pittsburgh, Nov. 6, 2000

