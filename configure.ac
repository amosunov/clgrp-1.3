#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.68])
AC_INIT(clgrp, 1.3, amosunov@uwaterloo.ca)
AM_INIT_AUTOMAKE(clgrp, 1.3)
AC_CONFIG_SRCDIR([clgrp.c])
AC_CONFIG_HEADERS([config.h])
AC_PREFIX_DEFAULT (/usr/local)
AC_PROG_CC([mpicc])

# Checks for programs.
AC_PROG_CC

# Checks for libraries.
AC_PROG_RANLIB

# Checks for header files.
AC_CHECK_HEADERS([stdio.h fcntl.h stdlib.h string.h sys/statvfs.h sys/time.h unistd.h])
AC_CHECK_HEADER(mpi.h, ,
  [AC_MSG_ERROR([Open MPI not found, see https://www.open-mpi.org/])],[[#if HAVE_MPI_H
# include <mpi.h>
#endif
]])
#AC_CHECK_HEADER(omp.h, ,
#  [AC_MSG_ERROR([OpenMP not found, see https://www.openmp.org/])],[[#if HAVE_OMP_H
# include <omp.h>
#endif
#]])

# Checks for typedefs, structures, and compiler characteristics.
AC_TYPE_OFF_T
AC_TYPE_SIZE_T

# Checks for library functions.
AC_FUNC_MALLOC
AC_FUNC_MMAP
AC_CHECK_FUNCS([floor ftruncate getpagesize gettimeofday memset mkdir munmap sqrt])
AC_CHECK_LIB(gmp, __gmpz_init, ,
  [AC_MSG_ERROR([GNU MP not found, see https://gmplib.org/])])
AC_CHECK_LIB(optarith, group_pow_init, ,
  [AC_MSG_ERROR([GNU MP not found, see https://github.com/maxwellsayles/liboptarith/])])
AC_CHECK_LIB(qform, qform_random_primeform, ,
  [AC_MSG_ERROR([GNU MP not found, see https://github.com/maxwellsayles/libqform/])])

AC_CONFIG_FILES([Makefile])
AC_OUTPUT
