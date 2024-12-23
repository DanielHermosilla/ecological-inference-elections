#include <R.h>
#include <R_ext/Rdynload.h>
#include <Rinternals.h>
#include <gsl/gsl_sf_bessel.h>
#include <stdio.h>

SEXP hello_gsl()
{
    printf("Hello, World from GSL!\n");

    // Ejemplo: calcular la función de Bessel de primer tipo J0(5.0)
    double x = 5.0;
    double result = gsl_sf_bessel_J0(x);

    printf("El valor de la función de Bessel J0(%.1f) es: %.5f\n", x, result);
    return R_NilValue;
}
