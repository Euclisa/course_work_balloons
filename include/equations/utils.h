#ifndef _EQUATIONS_UTILS_H
#define _EQUATIONS_UTILS_H

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <stdio.h>

#define MIN(a,b) (((a) < (b)) ? (a) : (b))


gsl_vector *vector_centred(const gsl_vector *v, const gsl_vector *center);

int center_from_points_and_radius(const gsl_vector *p1, const gsl_vector *p2, double r, gsl_vector *center);

int vectors_ang_clockwise(const gsl_vector *a, const gsl_vector *b, double *ang);

double area_between_vectors_triangle(const gsl_vector *a, const gsl_vector *b);

double area_segment(double ang, double r);

int point_from_alpha(double alpha, double norm, gsl_vector *center, gsl_vector *v);

void print_matrix(FILE *stream, const gsl_matrix *m);

void J_estimate(const gsl_vector *x0, void *params, gsl_matrix *J_est, int (*f)(const gsl_vector *x, void *params, gsl_vector *f));

void print_J_diff(FILE *stream, const gsl_vector *x, void *params, int (*f)(const gsl_vector *x, void *params, gsl_vector *f), int (*df)(const gsl_vector *x, void *params, gsl_matrix *df));

#endif