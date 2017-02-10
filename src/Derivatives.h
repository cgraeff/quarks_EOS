//
//  Derivatives.h
//  quarks EOS
//
//  Created by Clebson Graeff on 2017-02-19.
//  Copyright © 2016 Clebson Graeff. All rights reserved.
//

/* DerivativeWithOrder:
 *
 * Calculates the derivative of a function described
 * by a collection of coordinated pairs. The abscissas
 * are assumed to be regularly spaced. The derivatives
 * are calculated only at the same abscissas as the function.
 *
 *  - The function must be known in at least six points.
 *  - Steps must be constant (or you will get undefined behavior)
 *  - Higher order derivatives have less precision than lower orders
 *  - Derivatives on first two and last two points have bigger errors
 *    than the other points (see the reference for numbers).
 *  - The vector should have the same size and must be ordered
 *    with ascending abscissas
 *
 * Arguments:
 *  - order: order of the derivative (from 1 for the first derivative to
 *           5 for the fifth)
 *  - abscissas: vector with the points abscissas
 *  - ordinates: vector with the points ordinates
 *  - point_index: index of the value of abscissa at which the
 *                 derivative will be calculated
 *
 * Reference: Abramowitz & Stegun
*/

#include <gsl/gsl_vector.h>

double DerivativeWithOrder(int          order,
                           gsl_vector * abscissas,
                           gsl_vector * ordinates, 
                           int          point_index);
