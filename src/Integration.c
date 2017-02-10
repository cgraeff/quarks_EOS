
// This function should be removed as it simplify very little
#include <gsl/gsl_integration.h>

#include "Parameters.h"

double OnedimensionalIntegrator(gsl_function * F, double lower_limit, double upper_limit)
{
    gsl_integration_workspace * workspace =
    gsl_integration_workspace_alloc(parameters.fermi_dirac_integrals.max_interval_num);

    double integral = 0;
    double abserr;
    gsl_integration_qag(F,
                        lower_limit,
                        upper_limit,
                        parameters.fermi_dirac_integrals.abs_error,
                        parameters.fermi_dirac_integrals.rel_error,
                        parameters.fermi_dirac_integrals.max_sub_interval,
                        parameters.fermi_dirac_integrals.integration_key,
                        workspace,
                        &integral,
                        &abserr);

    gsl_integration_workspace_free(workspace);

    return integral;
}
