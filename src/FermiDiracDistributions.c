//
//  FermiDiracDistributions.c
//  quarks EOS
//
//  Created by Clebson Graeff on 2017-02-05.
//  Copyright © 2016 Clebson Graeff. All rights reserved.
//

#include <math.h>
#include <gsl/gsl_integration.h>
#include "Parameters.h"
#include "FermiDiracDistributions.h"
#include "Integration.h"

typedef struct _fermi_dirac_distrib_integrand{
    double mass;
    double chemical_potential;
    double temperature;
} fermi_dirac_distrib_integrand;


double FermiDiracDistributionFromDensityIntegralIntegrand(double momentum, void * parameters);
double FermiDiracDistributionIntegralFromGapEquationIntegrand(double momentum,
                                                              void * params);


double FermiDiracDistributionForParticles(double energy,
                                          double chemical_potential,
                                          double temperature)
{
    return 1.0 / (1.0 + exp((energy - chemical_potential)/temperature));
}

double FermiDiracDistributionForAntiparticles(double energy,
                                              double chemical_potential,
                                              double temperature)
{
    return 1.0 / (1.0 + exp((energy + chemical_potential)/temperature));
}

double FermiDiracDistributionFromDensityIntegral(double mass,
                                                 double renormalized_chemical_potential)
{
    fermi_dirac_distrib_integrand p;
    p.mass = mass;
    p.chemical_potential = renormalized_chemical_potential;
    p.temperature = parameters.finite_temperature.temperature;

    gsl_function F;
    F.function = &FermiDiracDistributionFromDensityIntegralIntegrand;
    F.params = &p;

  	double integral = OnedimensionalIntegrator(&F, 0.0, parameters.model.cutoff);

  	return integral;
}

double FermiDiracDistributionFromDensityIntegralIntegrand(double momentum, void * params)
{
    fermi_dirac_distrib_integrand * p = (fermi_dirac_distrib_integrand *) params;

    double E = sqrt(pow(momentum, 2.0) + pow(p->mass, 2.0));

    double particle_term = FermiDiracDistributionForParticles(E,
                                                              p->chemical_potential,
                                                              p->temperature);
    double antiparticle_term = FermiDiracDistributionForAntiparticles(E,
                                                                      p->chemical_potential,
                                                                      p->temperature);
    return (particle_term - antiparticle_term) * pow(momentum, 2.0);
}


double FermiDiracDistributionIntegralFromGapEquation(double mass,
                                                     double renormalized_chemical_potential)
{
    // Maybe it's possible to write the first term as the scalar density
    // at zero temperature

    fermi_dirac_distrib_integrand p;
    p.mass = mass;
    p.chemical_potential = renormalized_chemical_potential;
    p.temperature = parameters.finite_temperature.temperature;

    gsl_function F;
    F.function = &FermiDiracDistributionIntegralFromGapEquationIntegrand;
    F.params = &p;

	double integral = OnedimensionalIntegrator(&F, 0.0, parameters.model.cutoff);

    return integral;
}

double FermiDiracDistributionIntegralFromGapEquationIntegrand(double momentum,
                                                        void * params)
{
    fermi_dirac_distrib_integrand * p = (fermi_dirac_distrib_integrand *) params;

    double E = sqrt(pow(momentum, 2.0) + pow(p->mass, 2.0));

    double particle_term = FermiDiracDistributionForParticles(E,
                                                              p->chemical_potential,
                                                              p->temperature);
    double antiparticle_term = FermiDiracDistributionForAntiparticles(E,
                                                                      p->chemical_potential,
                                                                      p->temperature);
    return (1.0 - particle_term - antiparticle_term) * pow(momentum, 2.0) / E;
}

