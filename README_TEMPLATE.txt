#
#  = NAME
#
#  DAEPAK - A package of DAE solvers written by W.C. Rheinboldt
#
#  = VERSION #VERSION#
#
#   Released 11/06/2000
#
#  = DIRECTORY STRUCTURE
#
#  daepak-#VERSION#/lib holds solvers
#  daepak-#VERSION#/test holds drivers
#  daepak-#VERSION#/README_TEMPLATE.txt (this file)
#
#  = LIBRARIES
#
#  daepak-#VERSION#/cygdaepak-#VERSION#.dll compiled solver shared library
#  daepak-#VERSION#/libdaepak-#VERSION#.a compiled solver archive (ar) file
#  daepak-#VERSION#/libdaepak-#VERSION#G.a compiled solver archive (ar) file compiled with (-pg -g) compiler flags for gprof
#
#  = DEPENDENCIES
#
#DEPENDS#
#
#  = COMMENTS
#  == Fortran 95
#
#  == RunTime determined from unix date command before and after executing the driver. This timing is quite inaccurate.
#
#  == GPROF timings. See driver results in daepak-#VERSION#/test/{driver}_gmon.txt. gprof actually says no time spent
#     executing driver.
#
#  == Output results in daepak-#VERSION#/test/../{driver}_output.txt.
#
#  == Commented output. See driver results in daepak-#VERSION#/test/../{driver}.txt.
#
#  = OUTPUT at installation
#
#  == Command Line
#COMMAND#
#
#  == Output
#DRIVERS#
