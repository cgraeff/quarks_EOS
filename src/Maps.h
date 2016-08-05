//
//  Maps.h
//  quarks EOS
//
//  Created by Clebson Graeff on 2016-08-04.
//  Copyright © 2016 Clebson Graeff. All rights reserved.
//

/*
 * The following functions are intended as a device to visualize the regions of
 * a two variables function f(x, y) that are close to zero and calculate the values
 * of the variables that can make two such functions (e.g, f(x, y) and g(x, y)) be
 * close to zero simultaneously (this is given by the intersection of two such regions).
 * This is useful as a way to approximate the values of the two variables that are
 * roots to two functions simultaneously. This approximation can be used as a
 * starting point for a more sofisticated method of bi-dimensional root-finding.
 *
 */

#ifndef Maps_h
#define Maps_h

// The following function takes another function and verifies whether its absolute return value is
// less than a given tolerance, that is, whether |f(x,y)| <= tolerance is true
// INPUT:
//      double (*f)(double x, double y, void * params): a pointer to a function which takes two doubles and return a double
//      min_*, max_*, *_num_pts: minimum and maximum values of each variable, as well as the number of points
//                               in which the function will be probed in a given axis
//      tolerance: the maximum value the function may return (in absolute numbers) for the point to be
//                 registered
//      params: parameters to be passed for the function
//      show_progress: The execution of the map procedure may take a while if the map is too large
//                     (too many points), or if the provided function takes too long to calculate.
//                     In these ocasions, a percentual indicator may be shown to report progress.
// OUTPUT:
//      map_*: Return the values of x and y for each point that satisfies the condition |f(x,y)| <= tolerance
//             Those vectors must be pre-allocated and its sizes must be x_num_pts * y_num_pts (the vectors
//             must be able to hold all points). The i-th ordinate pair (x, y) is given by
//             (x, y)_i = (gsl_vector_get(output_x, i), gsl_vector_get(output_y, i))
//      num_points: returns the number of points that satisfies |f(x,y)| <= tolerance
//
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
                 int * num_points);

// The following routine takes two maps and calculates an approximation
// to the intersection of them.
//
// Before the approximation is made, the maps
// are looped over eliminating repeated values of x by taking the mean of the
// corresponding y values, combining multiple points in just one for a given
// value of x. After that, the pair of points (one from each map) that are
// separeted by the least distance is determined. From that point (and their
// corresponding first neighbours) two linear least squares fittings are made.
// The intersection of the two fittings is returned as an approximation for the
// interception of both maps.
//
// No steps are taken to verify the possibility of two intersections.
//
// INPUT:
//      Two maps calculated with MapFunction.
// OUTPUT:
//      x and y values of the (approximate) intersection point.
//
void IntersectionPointOfTwoMaps(gsl_vector * map1_x,
                                gsl_vector * map1_y,
                                int map1_num_points,
                                gsl_vector * map2_x,
                                gsl_vector * map2_y,
                                int map2_num_points,
                                double * x_intersection_return,
                                double * y_intersection_return);

#endif /* Maps_h */
