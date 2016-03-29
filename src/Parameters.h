//
//  Parameters.h
//  quarks EOS
//
//  Created by Clebson Graeff on 2016-03-17.
//  Copyright © 2016 Clebson Graeff. All rights reserved.
//

#ifndef Parameters_h
#define Parameters_h

typedef struct _parameters
{
    // Parameters for performing the calculation of the EOS
    int points_number;       // Number of points in which the above range will be divide into
    double minimum_density;     // (fm^-3)
    double maximum_density;     // (fm^-3)
    double proton_fraction;     // (no dimension)
    
    // Parameters used when writing the zeroed gap equation
    double gap_minimum_mass;
    double gap_maximum_mass;
    double gap_points_number;
    
/* Parameters for root finding in the solution of the Gap Equation and Vacuum Mass:
       - The root is assumed to be in the region bracketed by *_lower_bound and *_upper_bound;
 		- Two kinds of error are defined:
            - *_const_abs_error: An absolute error in MeV for the width of the region
                                 which contais the root
            - *_const_rel_error: Error as a fraction of the value of MIN(x_lower, x_upper),
                                 where x_lower and x_upper are the values that bracket the
                                 root after each iteration of the root findind algorithm
       - The root is accepted when
            |x_lower - x_upper| < *_const_abs_error
                                  + *_const_rel_error * MIN(|x_upper|, |x_lower|)
*/
    double gap_eq_solver_lower_bound;
    double gap_eq_solver_upper_bound;
    double gap_eq_solver_abs_error;
    double gap_eq_solver_rel_error;
    double gap_eq_solver_max_iterations;
    
    double vac_mass_det_lower_bound;
    double vac_mass_det_upper_bound;
    double vac_mass_det_abs_error;
    double vac_mass_det_rel_error;
    double vac_mass_det_max_iterations;
    
    // Parameterization variables
    double G_S;             // scalar-isoscalar coupling (fm^2)
    double cutoff;          // \Lambda (MeV)
    double bare_mass;       // (MeV)
	

} Parameters;

extern Parameters parameters;

void ParametersSetup(void);

#endif /* Parameters_h */
