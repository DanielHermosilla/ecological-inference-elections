# Initialize Autoconf
AC_INIT([ecological-inference-elections], [1.0], [daniel.hermosilla.r@ug.uchile.cl])
AC_PREREQ([2.69]) # Minimum Autoconf version
AM_INIT_AUTOMAKE([foreign])

# Compiler and Build Options
AC_PROG_CC                    # Find a C compiler
AC_CONFIG_MACRO_DIR([m4])     # Macro directory
AC_CONFIG_HEADERS([config.h]) # Generate a config.h file

# Set C Standard
AX_PROG_CC_C99                # Use C99/C11 standard if available

# Set Release Build Flags
CFLAGS="-Wall -Wextra -O3 -march=native -DNDEBUG -ffast-math -mtune=native -Xpreprocessor -fopenmp"
AC_SUBST(CFLAGS)

# Library and Include Checks
AC_CHECK_HEADERS([math.h stdlib.h string.h, stdio.h]) # Ensure standard headers exist
AC_CHECK_LIB([m], [main])                    # Math library (-lm)

# OpenBLAS
AC_CHECK_LIB([openblas], [cblas_dgemm], [], [AC_MSG_ERROR([OpenBLAS not found.])])
AC_CHECK_HEADERS([cblas.h], [], [AC_MSG_ERROR([OpenBLAS headers not found.])])

# GSL
AC_CHECK_LIB([gsl], [gsl_sf_bessel_J0], [], [AC_MSG_ERROR([GSL not found.])])
AC_CHECK_HEADERS([gsl/gsl_sf_bessel.h], [], [AC_MSG_ERROR([GSL headers not found.])])

# OpenMP
AC_CHECK_LIB([omp], [omp_get_num_threads], [], [AC_MSG_ERROR([OpenMP not found.])])
AC_CHECK_HEADERS([omp.h], [], [AC_MSG_ERROR([OpenMP headers not found.])])

# cJSON
AC_CHECK_LIB([cjson], [cJSON_Parse], [], [AC_MSG_ERROR([cJSON library not found.])])
AC_CHECK_HEADERS([cjson/cJSON.h], [], [AC_MSG_ERROR([cJSON headers not found.])])

# Source Files
AC_CONFIG_FILES([Makefile src/Makefile])
AC_OUTPUT

