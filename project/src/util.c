#include <R.h>
#include <R_ext/Rdynload.h>
#include <Rinternals.h>
#include <cblas.h>
#include <stdio.h>

// note: Refer to https://www.seehuhn.de/pages/linear.html#blas
// double m[] = {
// 3, 1, 3,
// 1, 5, 9,
// 2, 6, 5
// };

double getInitialP()
{

    /**
     * @brief Computes the initial probability of the EM algoritm.
     *
     * Given the observables results, it computes a convenient initial "p" value for initiating the
     * algorithm. Currently it supports the "uniform", "group_proportional" and "proportional" methods.
     *
     * @param[in] X Matrix of dimension (cxb) that stores the results of candidate "c" on ballot box "b".
     * @param[in] w Matrix of dimension (bxg) that stores the amount of votes from the demographic group "g".
     * @param[in] p_method The method for calculating the initial parameter. Currently it supports "uniform",
     * "group_proportional" and "proportional" methods.
     * @return Matrix of dimension (gxc) with the initial probability for each demographic group "g" voting for a given
     * candidate "c".
     * @note This should be used only that the first iteration of the EM-algorithm.
     * @warning Warnings about misuse or side effects.
     * @see Reference to other relevant functions or documentation.
     */
    return R_NilValue;
}

SEXP hello_gsl()
{
    printf("Hello, World from GSL!\n");

    // Ejemplo: calcular la función de Bessel de primer tipo J0(5.0)
    double x = 5.0;
    double result = gsl_sf_bessel_J0(x);

    printf("El valor de la función de Bessel J0(%.1f) es: %.5f\n", x, result);
    return R_NilValue;
}
