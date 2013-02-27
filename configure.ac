#AC_INIT([scrm],esyscmd([sh -c "date \"+%Y%m%d\" | tr -d '\n'"]),[staab@biologie.uni-muenchen.de])
AC_INIT([scrm],[beta],[staab@biologie.uni-muenchen.de])
dnl AC_CONFIG_SRCDIR([simu.c])
#AC_CONFIG_SRCDIR([main.cpp])
#AC_CONFIG_HEADER([config.h])

AM_INIT_AUTOMAKE
AC_PROG_INSTALL
AC_PREREQ
AC_CANONICAL_HOST

# Checks for programs.
#AC_PROG_CC
AC_PROG_CXX
AC_PROG_RANLIB

# Checks for libraries.
#AC_CHECK_LIB([m], [log])

# Checks for header files.
AC_HEADER_STDC
#AC_CHECK_HEADERS([float.h stdlib.h string.h unistd.h])
AC_LANG(C++) 
AC_CHECK_HEADERS([iostream string sstream fstream vector iomanip valarray math.h stdio.h ctime boost])

# Checks for typedefs, structures, and compiler characteristics.
AC_C_CONST
AC_C_INLINE
AC_TYPE_SIZE_T

# Checks for library functions.
AC_FUNC_MALLOC
AC_FUNC_REALLOC

AC_CHECK_FUNCS([floor])
AC_CHECK_FUNCS([pow])
AC_HEADER_STDBOOL




# Check for presence of pdfLaTeX
AC_CHECK_PROG(PDFLATEX, pdflatex, pdflatex)
AC_CHECK_PROG(BIBTEX, bibtex,bibtex)
if test -z "$PDFLATEX"; then
  AC_MSG_WARN([Unable to create PDF version of the user manual.])
fi

AM_CONDITIONAL([HAVE_PDFLATEX], test -n "$PDFLATEX")

AC_CONFIG_FILES([Makefile src/Makefile doc/Makefile])
AC_OUTPUT

