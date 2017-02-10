//
//  Derivatives.c
//  quarks EOS
//
//  Created by Clebson Graeff on 2017-02-19.
//  Copyright © 2016 Clebson Graeff. All rights reserved.
//

#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_vector.h>

double CoefficientForDerivativeOrder(int order, int point, int index);

double DerivativeWithOrder(int order, gsl_vector * abscissas, gsl_vector * ordinates, int point_index)
{
    // Supported derivatives range
    const int MIN_ORDER = 1;
    const int MAX_ORDER = 5;
    
    if (order > MAX_ORDER || order < MIN_ORDER){
        printf("DerivativeWithOrder: Derivative order %d unsupported. Acceptable values range from %d to %d.\n",
               order,
               MIN_ORDER,
               MAX_ORDER);
        abort();
    }

    // The method uses a set of six points.
    const int NUM_PTS_SET = 6;

    if (abscissas->size < NUM_PTS_SET){
        printf("DerivativeWithOrder: A minimum of six points is required.\n");
        abort();
    }
    
    int num_data_pts = abscissas->size;
    
    // Determine the position of the requested point in
    // the set. This is necessary to calculate proper
    // derivatives for the first two and the last two
    // points in the set. For all others we will treat
    // the point as the third in the set and use the
    // neighbouring points. This ensures that the error
    // is kept as low as possible.
    
    int j;
    
    if (point_index <= 2){
        j = point_index;
    }
    else if (point_index == (num_data_pts - 1)){
        j = 5;
    }
    else if (point_index == (num_data_pts - 2)){
        j = 4;
    }
    else if (point_index == (num_data_pts - 3)){
        j = 3;
    }
    else {
        j = 2;
    }
   
    double sum = 0.0;
    for (int i = 0; i < NUM_PTS_SET; i++)
        sum += CoefficientForDerivativeOrder(order, j, i)
               * gsl_vector_get(ordinates, point_index + i - j);
    
    double step = gsl_vector_get(abscissas, 1) - gsl_vector_get(abscissas, 0);

    return gsl_sf_fact(order) * sum
           / (gsl_sf_fact(NUM_PTS_SET - 1) * pow(step, order));
}

double CoefficientForDerivativeOrder(int order, int point, int index){
    
    if (order > 5 || order < 1){
        printf("CoefficientForDerivativeOrder: order out of bounds.\n");
        abort();
    }
    
    if (point > 5 || point < 0){
        printf("CoefficientForDerivativeOrder: point out of bounds.\n");
        abort();
    }
    
    if (index > 5 || index < 0){
        printf("CoefficientForDerivativeOrder: index out of bounds.\n");
        abort();
    }
        
    const double A[5][6][6] =
        {
            {
                {-274.0,  600.0, -600.0,  400.0, -150.0,  24.0},
                { -24.0, -130.0,  240.0, -120.0,   40.0,  -6.0},
                {   6.0,  -60.0,  -40.0,  120.0,  -30.0,   4.0},
                {  -4.0,   30.0, -120.0,   40.0,   60.0,  -6.0},
                {   6.0,  -40.0,  120.0, -240.0,  130.0,  24.0},
                { -24.0,  150.0, -400.0,  600.0, -600.0, 274.0}
            },

            {
                { 225.0, -770.0, 1070.0, -780.0,  305.0, -50.0},
                {  50.0,  -75.0,  -20.0,   70.0,  -30.0,   5.0},
                {  -5.0,   80.0, -150.0,   80.0,   -5.0,   0.0},
                {   0.0,   -5.0,   80.0, -150.0,   80.0,  -5.0},
                {   5.0,  -30.0,   70.0,  -20.0,  -75.0,  50.0},
                { -50.0,  305.0, -780.0, 1070.0, -770.0, 225.0}
            },

            {
                { -85.0,  355.0, -590.0,  490.0, -205.0,  35.0},
                { -35.0,  125.0, -170.0,  110.0,  -35.0,   5.0},
                {  -5.0,   -5.0,   50.0,  -70.0,   35.0,  -5.0},
                {   5.0,  -35.0,   70.0,  -50.0,    5.0,   5.0},
                {  -5.0,   35.0, -110.0,  170.0, -125.0,  35.0},
                { -35.0,  205.0, -490.0,  590.0, -355.0,  85.0}
            },

            {
                {  15.0,  -70.0,  130.0, -120.0,   55.0, -10.0},
                {  10.0,  -45.0,   80.0,  -70.0,   30.0,  -5.0},
                {   5.0,  -20.0,   30.0,  -20.0,    5.0,   0.0},
                {   0.0,    5.0,  -20.0,   30.0,  -20.0,   5.0},
                {  -5.0,   30.0,  -70.0,   80.0,  -45.0,  10.0},
                { -10.0,   55.0, -120.0,  130.0,  -70.0,  15.0}
            },

            {
                { -1.0,     5.0,  -10.0,   10.0,   -5.0,   1.0},
                { -1.0,     5.0,  -10.0,   10.0,   -5.0,   1.0},
                { -1.0,     5.0,  -10.0,   10.0,   -5.0,   1.0},
                { -1.0,     5.0,  -10.0,   10.0,   -5.0,   1.0},
                { -1.0,     5.0,  -10.0,   10.0,   -5.0,   1.0},
                { -1.0,     5.0,  -10.0,   10.0,   -5.0,   1.0}
            }
        };

    return A[order - 1][point][index];
}

