#ifndef FermiDiracDistributions_h
#define FermiDiracDistributions_h

typedef struct _FermiDiracIntegralsParameters{
  	int max_interval_num;
   	int integration_key;
   	int max_sub_interval;
   	double abs_error;
   	double rel_error;
} FermiDiracIntegralsParameters;

double FermiDiracDistributionFromDensityIntegral(double mass,
                                                 double renormalized_chemical_potential);
double FermiDiracDistributionIntegralFromGapEquation(double mass,
                                                     double renormalized_chemical_potential);

#endif /* FermiDiracDistributions_h */
