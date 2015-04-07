AC_INIT([scrm], [1.5.0],[https://github.com/paulstaab/scrm/issues])
AM_INIT_AUTOMAKE([subdir-objects -Wall -Werror foreign])

# Use -O3 as default optimization level
: ${CXXFLAGS="-O3"}

# Suppress output at compilation
m4_ifndef([AM_SILENT_RULES], [m4_define([AM_SILENT_RULES],[])])
AM_SILENT_RULES([yes])

# Load macros in 'm4'-dir
AC_CONFIG_MACRO_DIR([m4])

AC_PROG_INSTALL
AC_PREREQ
AC_CANONICAL_HOST

# Checks for programs.
AC_PROG_RANLIB

# Check for C++11
AX_CXX_COMPILE_STDCXX_11(,mandatory)

# Checks for libraries
AC_CHECK_LIB(cppunit,TestCase,[])

# Checks for header files for scrm.
AC_HEADER_STDC
AC_LANG(C++) 

# Checks for typedefs, structures, and compiler characteristics.
AC_C_CONST
AC_C_INLINE
AC_TYPE_SIZE_T

# Enable Link-time optimization if supported (gcc only)
if test x$CXX = xg++; then
  AX_CHECK_COMPILE_FLAG([-flto], [OPT_CXXFLAGS="-flto"], [], [-Werror])
  AC_SUBST(OPT_CXXFLAGS)
fi

# Enable Warnings
AX_CXXFLAGS_WARN_ALL

AC_CONFIG_FILES([Makefile])
AC_OUTPUT
