DAEPAK
======

DAE solver based on geometric techniques

W.C. Rheinboldt wrote the original solvers and supporting libraries in Fortran 95 released on 11/06/2000.
The initial commit is exactly these files unpacked into a directory structure for development.

No one is actively maintaining DAEPAK. I value this technology and created this git repository
to bring awareness, for continue maintenance, and for free distribution of this DAE solver library.


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
DAEPAK is written is Fortran 95 and requires a fortran compiler.

A choice on Windows in the MinGW environment is the [MinGW fortran compiler](http://www.mingw.org)
and there is (video) support to install this compiler. A recent installer installed GNU Fortran (GCC) version 4.8.1.
Another choice on Windows in [CygWin](https://www.cygwin.com) which has a package manager.

A [Ruby](http://www.ruby-doc.org) installation which includes the rake gem is required to utilize and execute the
rakefile rakefile.rb. This is worth doing.

A choice for unit testing is [PFunit](http://en.wikipedia.org/wiki/PFUnit). I plan to create unit tests in the future.

You are ready to build DAEPAK assuming your system is running GNU Fortran (CGG) version 4.8.1 and
Ruby 1.9.3p429 and Rake.
### Build the library
-  $rake clobber
-  $rake sources

Check that a library file libdaepak.a is in the top-level directory

### Build and run the library examples
-  $rake examples

Check out the resulting .LOG and .TXT files in the examples/EXAMPLE directory.

### Build your example
-  Follow a prototype --- if need be --- to setup your problem
-  Copy libdaepak.a to your examples directory
-  $ar x libdaepak.a
-  $gfortran -Wall -pedantic -c YourExamplesFile.f95 -o YourExamplesFile.o
-  $gfortran -Wall -pedantic ALLSOURCEFILES.o -o YourExamplesFile.x
-  $./YourExamplesFile.x

Roadmap for the Future
======================
-  Incorporate unit testing
-  Move examples outside the library
-  Move linear algebra support to its own library
-  Move ordinary differential equation solver support to its own library
-  Move manifold support to its own library
-  Incorporate the examples from prior F77 release

