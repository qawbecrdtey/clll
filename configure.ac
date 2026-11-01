AC_PREREQ([2.71])
AC_INIT([clll], [1.0.0], [qawbecrdtey@kaist.ac.kr], [clll])
AC_CONFIG_SRCDIR([src/lll.c])
AC_CONFIG_AUX_DIR([build-aux])
AC_CONFIG_MACRO_DIR([m4])

AM_INIT_AUTOMAKE([1.16 foreign subdir-objects -Wall -Werror])
AM_SILENT_RULES([yes])

dnl Libtool version-info (current:revision:age). It tracks the ABI, not the
dnl release number.
AC_SUBST([CLLL_VERSION_INFO], [0:0:0])

AC_PROG_CC
AC_PROG_SED
AM_PROG_AR
LT_INIT

dnl ---------------------------------------------------------------------------
dnl C99
dnl ---------------------------------------------------------------------------

AC_MSG_CHECKING([for a C99 mode])
clll_std_cflags=
clll_saved_cflags=$CFLAGS
CFLAGS="$CFLAGS -std=c99"
AC_COMPILE_IFELSE(
  [AC_LANG_PROGRAM(
     [[#include <stdbool.h>
       #include <stdint.h>
       static inline int64_t clll_probe_twice(int64_t x) { return 2 * x; }]],
     [[bool ok = true;
       for (int i = 0; i < 2; i++) {
           ok = ok && clll_probe_twice(i) == 2 * i;
       }
       return !ok;]])],
  [clll_std_cflags=-std=c99])
CFLAGS=$clll_saved_cflags
AS_IF([test -z "$clll_std_cflags"],
      [AC_MSG_RESULT([none])
       AC_MSG_ERROR([a C99 compiler is required])],
      [AC_MSG_RESULT([$clll_std_cflags])])
AC_SUBST([CLLL_STD_CFLAGS], [$clll_std_cflags])

dnl ---------------------------------------------------------------------------
dnl Warnings
dnl ---------------------------------------------------------------------------

dnl Warnings never fail the build, since a newer compiler may add new ones.
dnl Use CFLAGS='-g -O2 -Werror' to make them errors, as make distcheck does.
CLLL_WARN_CFLAGS="-Wall -Wextra -pedantic"
AC_SUBST([CLLL_WARN_CFLAGS])

dnl ---------------------------------------------------------------------------
dnl Calls into GMP
dnl ---------------------------------------------------------------------------

dnl -fno-plt calls GMP through the GOT instead of a PLT stub, which is a little
dnl faster. CFLAGS='-g -O2 -fplt' turns it off again.
clll_saved_cflags=$CFLAGS
AC_MSG_CHECKING([whether $CC accepts -fno-plt])
CFLAGS="$clll_saved_cflags -Werror -fno-plt"
AC_LINK_IFELSE([AC_LANG_PROGRAM([], [])],
  [AC_MSG_RESULT([yes])
   clll_plt_cflags=-fno-plt],
  [AC_MSG_RESULT([no])
   clll_plt_cflags=])
CFLAGS=$clll_saved_cflags
AC_SUBST([CLLL_PLT_CFLAGS], [$clll_plt_cflags])

dnl Assertions are on unless NDEBUG is defined, for example with
dnl CPPFLAGS=-DNDEBUG. This check only reports which one it is.
AC_MSG_CHECKING([whether assertions are checked])
AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[#ifdef NDEBUG
#error NDEBUG is defined
#endif]], [])],
  [clll_assert=on],
  [clll_assert=off])
AC_MSG_RESULT([$clll_assert])

dnl ---------------------------------------------------------------------------
dnl GMP
dnl ---------------------------------------------------------------------------

AC_CHECK_HEADER([gmp.h], [],
  [AC_MSG_ERROR([gmp.h was not found; install GMP's development package (libgmp-dev, gmp-devel or dev-libs/gmp)])])

dnl mpz_init is a macro, so check for the real symbol, __gmpz_init.
AC_CHECK_LIB([gmp], [__gmpz_init], [GMP_LIBS=-lgmp],
  [AC_MSG_ERROR([libgmp was not found; install GMP's development package (libgmp-dev, gmp-devel or dev-libs/gmp)])])
AC_SUBST([GMP_LIBS])

dnl The public headers use mpz_t, so clll.pc lists GMP as a requirement. GMP has
dnl shipped gmp.pc only since 6.2, so older versions are linked with -lgmp.
AC_PATH_PROG([PKG_CONFIG], [pkg-config])
AS_IF([test -n "$PKG_CONFIG" && $PKG_CONFIG --exists gmp],
      [CLLL_PC_REQUIRES=gmp
       CLLL_PC_LIBS=],
      [CLLL_PC_REQUIRES=
       CLLL_PC_LIBS=-lgmp])
AC_SUBST([CLLL_PC_REQUIRES])
AC_SUBST([CLLL_PC_LIBS])

AC_CONFIG_FILES([Makefile clll.pc cmake/clllConfigVersion.cmake])
AC_OUTPUT

AS_IF([test -n "$CLLL_PC_REQUIRES"], [clll_pc_gmp="Requires: gmp"], [clll_pc_gmp="Libs: -lgmp"])

AC_MSG_NOTICE([

  $PACKAGE_NAME $PACKAGE_VERSION

    compiler          $CC $CLLL_STD_CFLAGS
    warnings          $CLLL_WARN_CFLAGS
    optimization      ${CLLL_PLT_CFLAGS:+$CLLL_PLT_CFLAGS }$CFLAGS
    assertions        $clll_assert
    shared / static   $enable_shared / $enable_static
    GMP in clll.pc    $clll_pc_gmp
    prefix            $prefix
])
