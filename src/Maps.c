//
//  Maps.c
//  quarks EOS
//
//  Created by Clebson Graeff on 2016-08-04.
//  Copyright © 2016 Clebson Graeff. All rights reserved.
//

#include <math.h>
#include <gsl/gsl_vector.h>

#include "Constants.h"
#include "Parameters.h"
#include "CommandlineOptions.h"

#include "Maps.h"

/*
 * Internal functions interface declarations
 */

// Takes a map and discard all points outside a given range.
// The output (clipped_map_*) should be pre-allocated with the
// same size as the input, as it's possible that all data is in
// the given range.
void ClipMapToXRange(gsl_vector * map_x,
                     gsl_vector * map_y,
                     int map_num_points,
                     double x_range_start,
                     double x_range_end,
                     gsl_vector * clipped_map_x,
                     gsl_vector * clipped_map_y,
                     int * clipped_map_num_points);

// Loop looking for repeated values of x. The multiple points
// will be reduced to one by taking the mean of the various values for y.
// The output (reduced_map_*) must be pre-allocated to the same size as
// the input as it may be possible that there are no points with the same x.
void ReduceMultiplicityOfMapPoints(gsl_vector * map_x,
                                   gsl_vector * map_y,
                                   int map_num_points,
                                   gsl_vector * reduced_map_x,
                                   gsl_vector * reduced_map_y,
                                   int * reduced_map_num_points);

// Given two reduced maps, determine the index in each map of the points
// that are closest to the intersection of the maps.
void IndexOfClosestPointsForTwoReducedMaps(gsl_vector * map1_x,
                                           gsl_vector * map1_y,
                                           int map1_num_points,
                                           gsl_vector * map2_x,
                                           gsl_vector * map2_y,
                                           int map2_num_points,
                                           int * map1_closest_index_return,
                                           int * map2_closest_index_return);

// Given three points (x_1, y_1), (x_2, y_2), and (x_3, y_3),
// calculate the coefficients of the y = A + B * x line by
// least squares fitting
void LinearFitOfThreePoints(double x_1, double y_1,
                            double x_2, double y_2,
                            double x_3, double y_3,
                            double * A, double * B);

/*
 * External functions implementation
 */

void MapFunction(double (*f)(double x, double y, void * params),
                 double min_x,
                 double max_x,
                 int x_num_pts,
                 double min_y,
                 double max_y,
                 int y_num_pts,
                 double tolerance,
                 void * params,
                 bool show_progress,
                 gsl_vector * map_x,
                 gsl_vector * map_y,
                 int * num_points)
{
    double x_step = (max_x - min_x) / (double)(x_num_pts - 1);
    double y_step = (max_y - min_y) / (double)(y_num_pts - 2);

    int progress_counter = 0;
    int pts_counter = 0;
    double x = min_x;
    double y = min_y;

    for (int i = 0; i < x_num_pts; i++) {
        for (int j = 0; j < y_num_pts; j++){

            double val = f(x, y, params);

            if (fabs(val) <= tolerance){
                gsl_vector_set(map_x, pts_counter, x);
                gsl_vector_set(map_y, pts_counter, y);
                pts_counter++;
            }

            if (show_progress)
                if ((100 * progress_counter) % (x_num_pts * y_num_pts) == 0){
                    printf("\rMap %d %%", (100 * progress_counter) / (x_num_pts * y_num_pts));
                    fflush(stdout);
                }

            progress_counter++;
            y += y_step;
        }

        x += x_step;
    }

    *num_points = pts_counter;
    return;
}

void IntersectionPointOfTwoMaps(gsl_vector * map1_x,
                                gsl_vector * map1_y,
                                int map1_num_points,
                                gsl_vector * map2_x,
                                gsl_vector * map2_y,
                                int map2_num_points,
                                double * x_intersection_return,
                                double * y_intersection_return)
{
    //
    // Determine common range in x axis and rebuild
	// the index
	//

    double common_interval_start = fmax(gsl_vector_get(map1_x, 1),
                                        gsl_vector_get(map2_x, 1));
    double common_interval_end = fmin(gsl_vector_get(map1_x, map1_num_points),
                                      gsl_vector_get(map2_x, map2_num_points));

    gsl_vector * common_map1_x = gsl_vector_alloc(map1_num_points);
    gsl_vector * common_map1_y = gsl_vector_alloc(map1_num_points);
    int common_map1_num_points;

    gsl_vector * common_map2_x = gsl_vector_alloc(map2_num_points);
    gsl_vector * common_map2_y = gsl_vector_alloc(map2_num_points);
    int common_map2_num_points;

    ClipMapToXRange(map1_x,
                    map1_y,
                    map1_num_points,
                    common_interval_start,
                    common_interval_end,
                    common_map1_x,
                    common_map1_y,
                    &common_map1_num_points);
    ClipMapToXRange(map2_x,
                    map2_y,
                    map2_num_points,
                    common_interval_start,
                    common_interval_end,
                    common_map2_x,
                    common_map2_y,
                    &common_map2_num_points);

    //
    // Reduce multiple points
    //

    gsl_vector * red_map1_x = gsl_vector_alloc(common_map1_num_points);
    gsl_vector * red_map1_y = gsl_vector_alloc(common_map1_num_points);
    int red_map1_num_points;

    gsl_vector * red_map2_x = gsl_vector_alloc(common_map2_num_points);
    gsl_vector * red_map2_y = gsl_vector_alloc(common_map2_num_points);
    int red_map2_num_points;

    ReduceMultiplicityOfMapPoints(common_map1_x,
                                  common_map1_y,
                                  common_map1_num_points,
                                  red_map1_x,
                                  red_map1_y,
                                  &red_map1_num_points);

    ReduceMultiplicityOfMapPoints(common_map2_x,
                                  common_map2_y,
                                  common_map2_num_points,
                                  red_map2_x,
                                  red_map2_y,
                                  &red_map2_num_points);

    // Discard vectors that will not be used anymore
    gsl_vector_free(common_map1_x);
    gsl_vector_free(common_map1_y);
    gsl_vector_free(common_map2_x);
    gsl_vector_free(common_map2_y);

  	//
	// Search for the indexes of the points of the maps that have the minimal
    // distance between them.
	//

    int map1_closest_index;
    int map2_closest_index;

    IndexOfClosestPointsForTwoReducedMaps(red_map1_x,
                                          red_map1_y,
                                          red_map1_num_points,
                                          red_map2_x,
                                          red_map2_y,
                                          red_map2_num_points,
                                          &map1_closest_index,
                                          &map2_closest_index);


    // Calculate a linear fit of three points (the closest and its neighbours)
    // for each map; The intersection of the maps is given by the intersection of the
    // linear fits.
    double x_1 = gsl_vector_get(red_map1_x, map1_closest_index - 1);
    double x_2 = gsl_vector_get(red_map1_x, map1_closest_index);
    double x_3 = gsl_vector_get(red_map1_x, map1_closest_index + 1);
    double y_1 = gsl_vector_get(red_map1_y, map1_closest_index - 1);
    double y_2 = gsl_vector_get(red_map1_y, map1_closest_index);
    double y_3 = gsl_vector_get(red_map1_y, map1_closest_index + 1);

    double A_1;
    double B_1;

    LinearFitOfThreePoints(x_1, y_1,
                           x_2, y_2,
                           x_3, y_3,
                           &A_1, &B_1);

    x_1 = gsl_vector_get(red_map2_x, map2_closest_index - 1);
    x_2 = gsl_vector_get(red_map2_x, map2_closest_index);
    x_3 = gsl_vector_get(red_map2_x, map2_closest_index + 1);
    y_1 = gsl_vector_get(red_map2_y, map2_closest_index - 1);
    y_2 = gsl_vector_get(red_map2_y, map2_closest_index);
    y_3 = gsl_vector_get(red_map2_y, map2_closest_index + 1);

    double A_2;
    double B_2;

    LinearFitOfThreePoints(x_1, y_1,
                           x_2, y_2,
                           x_3, y_3,
                           &A_2, &B_2);

    double x_intersection = (A_2 - A_1) / (B_1 - B_2);
    double y_intersection = A_1 + B_1 * x_intersection;

    // Free vectors
    gsl_vector_free(red_map1_x);
    gsl_vector_free(red_map1_y);
    gsl_vector_free(red_map2_x);
    gsl_vector_free(red_map2_y);

    // return values
    *x_intersection_return = x_intersection;
    *y_intersection_return = y_intersection;

    return;
}


/*
 * Internal functions implementations
 */

void ClipMapToXRange(gsl_vector * map_x,
                     gsl_vector * map_y,
                     int map_num_points,
                     double x_range_start,
                     double x_range_end,
                     gsl_vector * clipped_map_x,
                     gsl_vector * clipped_map_y,
                     int * clipped_map_num_points)
{
    int count = 0;

    for (int i = 0; i < map_num_points; i++){

        double x = gsl_vector_get(map_x, i);

        if (x >= x_range_start){
            gsl_vector_set(clipped_map_x, count, x);
            gsl_vector_set(clipped_map_y, count, gsl_vector_get(map_y, i));
            count++;
        }

        if (x > x_range_end)
            break;
    }

    *clipped_map_num_points = count;

    return;
}

void ReduceMultiplicityOfMapPoints(gsl_vector * map_x,
                                   gsl_vector * map_y,
                                   int map_num_points,
                                   gsl_vector * reduced_map_x,
                                   gsl_vector * reduced_map_y,
                                   int * reduced_map_num_points)
{
    int reduced_map_count = 0;

    // Use first point as a reference.
    double x = gsl_vector_get(map_x, 0);
    double y_sum = gsl_vector_get(map_y, 0);
    int multiplicity = 1;

    for (int i = 1; i < map_num_points; i++){

        double x_i = gsl_vector_get(map_x, i);

	  	// if the value of x is repeated, add the
		// value of y and increase the multiplicity
        if (x_i == x){
            y_sum += gsl_vector_get(map_y, i);
            multiplicity++;
        }
        else{
            // x has changed, save old x and corresponding mean y,
			// increase the number of reduced points
            gsl_vector_set(reduced_map_x, reduced_map_count, x);
            gsl_vector_set(reduced_map_y, reduced_map_count, y_sum / (double)multiplicity);
            reduced_map_count++;

            // prepare next iteration
            x = x_i;
            y_sum = gsl_vector_get(map_y, i);
            multiplicity = 1;
        }
    }
    // save last point
    gsl_vector_set(reduced_map_x, reduced_map_count, x);
    gsl_vector_set(reduced_map_y, reduced_map_count, y_sum / (double)multiplicity);
    reduced_map_count++;

    *reduced_map_num_points = reduced_map_count;

    return;
}

void IndexOfClosestPointsForTwoReducedMaps(gsl_vector * map1_x,
                                           gsl_vector * map1_y,
                                           int map1_num_points,
                                           gsl_vector * map2_x,
                                           gsl_vector * map2_y,
                                           int map2_num_points,
                                           int * map1_closest_index_return,
                                           int * map2_closest_index_return)
{
    // Start with first point
    int map1_closest_index = 0;
    int map2_closest_index = 0;

    double x_1 = gsl_vector_get(map1_x, 0);
    double y_1 = gsl_vector_get(map1_y, 0);
    double x_2 = gsl_vector_get(map2_x, 0);
    double y_2 = gsl_vector_get(map2_y, 0);

    double min_distance = sqrt(pow(y_2 - y_1, 2.0) + pow(x_2 - x_1, 2.0));

    // loop over all remaining points looking for the least distance pair
    for (int i = 1; i < map1_num_points; i++){

        x_1 = gsl_vector_get(map1_x, i);
        y_1 = gsl_vector_get(map1_y, i);

        for (int j = 1; j < map2_num_points; j++){

            x_2 = gsl_vector_get(map2_x, 0);
            y_2 = gsl_vector_get(map2_y, 0);

            double distance = sqrt(pow(y_2 - y_1, 2.0) + pow(x_2 - x_1, 2.0));

            if (distance < min_distance){
                min_distance = distance;
                map1_closest_index = i;
                map2_closest_index = j;
            }

        }
    }

    *map1_closest_index_return = map1_closest_index;
    *map2_closest_index_return = map2_closest_index;

    return;
}

void LinearFitOfThreePoints(double x_1, double y_1,
                            double x_2, double y_2,
                            double x_3, double y_3,
                            double * A, double * B)
{
    const double num_points = 3.0;

    double sum_x_sq = pow(x_1, 2.0) + pow(x_2, 2.0) + pow(x_3, 2.0);
    double sum_x = x_1 + x_2 + x_3;
    double sum_xy = x_1 * y_1 + x_2 * y_2 + x_3 * y_3;
    double sum_y = y_1 + y_2 + y_3;

    double delta = num_points * sum_x_sq - pow(sum_x, 2.0);

    *A = (sum_x_sq * sum_y - sum_xy) / delta;
    *B = (num_points * sum_xy - sum_x * sum_y) / delta;
}
