AC_INIT([scrm], [1.0-beta2],[https://github.com/paulstaab/scrm/issues])
AM_INIT_AUTOMAKE([subdir-objects -Wall -Werror foreign])
dnl AC_CONFIG_SRCDIR([simu.c])
#AC_CONFIG_SRCDIR([main.cpp])
#AC_CONFIG_HEADER([config.h])

m4_ifndef([AM_SILENT_RULES], [m4_define([AM_SILENT_RULES],[])])
AM_SILENT_RULES([yes])

AC_PROG_INSTALL
AC_PREREQ
AC_CANONICAL_HOST

# Checks for programs.
AC_PROG_RANLIB

# Checks for libraries
#AC_CHECK_LIB(cppunit,TestCase,[],AC_MSG_WARN([cppunit is not installed]))
#AC_CHECK_LIB(cppunit,TestCase,[],AC_MSG_WARN([User define 'cppunit' location. e.g. LDFLAGS=-L/opt/local/lib]))
AC_CHECK_LIB(cppunit,TestCase,[],LDFLAGS=-L/opt/local/lib)

# Checks for header files for scrm.
AC_HEADER_STDC
AC_LANG(C++) 


#
# Platform specific setup
#
#############################
AC_CANONICAL_HOST
# Check for which host we are on and setup a few things
# specifically based on the host
case $host_os in
  darwin* )
    AC_PROG_CXX([clang])
    export CXX=clang++
    export CXXFLAGS=-I/opt/local/include
        ;;
  linux*)
    AC_PROG_CXX
    # Checks for the compiler viersion
    GCC_VERSION_MAJOR=$(g++ -dumpversion | cut -d'.' -f1)
    GCC_VERSION_MINOR=$(g++ -dumpversion | cut -d'.' -f2)
    
    if test "$GCC_VERSION_MAJOR" -ge 4
    then
        if test "$GCC_VERSION_MINOR" -lt 6
        then
            AC_MSG_ERROR([Require g++ compiler 4.6+])
        fi
    else
        AC_MSG_ERROR([Require g++ compiler 4.6+])
    fi 

    ;;
   *BSD*)
        # Do something specific for BSD
        ;;
    *)
        #Default Case
        AC_MSG_ERROR([Your platform is not currently supported])
        ;;
esac

#AM_CONDITIONAL([DARWIN_HOST], test $host_os = darwin)
#AM_CONDITIONAL([LINUX_HOST], test $host_os = linux)



# Checks for the std=c++0x option
#AX_CXX_COMPILE_STDCXX_0X

# Checks for the std=c++11 option
#AX_CXX_COMPILE_STDCXX_11([ext|noext],[mandatory|optional])
#AX_CXX_COMPILE_STDCXX_11()

# Checks for typedefs, structures, and compiler characteristics.
AC_C_CONST
AC_C_INLINE
AC_TYPE_SIZE_T

# Checks for libraries.
#AC_CHECK_LIB([m], [log])

# Checks for library functions.
AC_FUNC_MALLOC
AC_FUNC_REALLOC

#AC_CHECK_FUNCS([floor])
#AC_CHECK_FUNCS([pow])
#AC_HEADER_STDBOOL

m4_include([m4/ac_doxygen.m4])
DX_INIT_DOXYGEN(scrm, Doxyfile, doc/doxygen)
DX_HTML_FEATURE(ON)
DX_CHM_FEATURE(OFF)
DX_CHI_FEATURE(OFF)
DX_MAN_FEATURE(OFF)
DX_RTF_FEATURE(OFF)
DX_XML_FEATURE(OFF)
DX_PDF_FEATURE(OFF)
DX_PS_FEATURE(OFF)


# Check for presence of pdfLaTeX
#AC_CHECK_PROG(PDFLATEX, pdflatex, pdflatex)
#AC_CHECK_PROG(BIBTEX, bibtex,bibtex)
#if test -z "$PDFLATEX"; then
#  AC_MSG_WARN([Unable to create PDF version of the user manual.])
#fi

#AM_CONDITIONAL([HAVE_PDFLATEX], test -n "$PDFLATEX")

# Check for presence of pandoc
AC_CHECK_PROG(PANDOC, pandoc, pandoc)
if test -z "$PANDOC"; then
  AC_MSG_WARN([Unable to create the user manual.])
fi

AM_CONDITIONAL([HAVE_PANDOC], test -n "$PANDOC")

AC_CONFIG_FILES([Makefile])
AC_OUTPUT
