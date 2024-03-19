#include <equations/utils.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_sf.h>
#include <math.h>


inline gsl_vector *vector_centred(const gsl_vector *v, const gsl_vector *center)
{
    gsl_vector *res = gsl_vector_alloc(2);
    gsl_vector_memcpy(res,v);
    gsl_vector_sub(res,center);
    return res;
}


int center_from_points_and_radius(const gsl_vector *p1, const gsl_vector *p2, double r, gsl_vector *center)
{
    gsl_vector *v12 = vector_centred(p1,p2);
    gsl_vector_scale(v12,0.5);
    double x = gsl_vector_get(v12,0);
    double y = gsl_vector_get(v12,1);
    if(x > 0)
    {
        gsl_vector_set(center,0,y);
        gsl_vector_set(center,1,-x);
    }
    else
    {
        gsl_vector_set(center,0,-y);
        gsl_vector_set(center,1,x);
    }
    double rv_norm = gsl_blas_dnrm2(center);
    double target_norm = sqrt(gsl_pow_2(r) - gsl_pow_2(rv_norm));
    gsl_vector_scale(center,target_norm/rv_norm);
    gsl_vector_add(center,p2);
    gsl_vector_add(center,v12);

    gsl_vector_free(v12);

    return GSL_SUCCESS;
}


int vectors_ang_clockwise(const gsl_vector *a, const gsl_vector *b, double *ang)
{
    double dot;
    gsl_blas_ddot(a,b,&dot);

    double a1 = gsl_vector_get(a,0);
    double a2 = gsl_vector_get(a,1);
    double b1 = gsl_vector_get(b,0);
    double b2 = gsl_vector_get(b,1);
    double det = a1*b2 - a2*b1;

    *ang = -atan2(det,dot);
    if(*ang < 0)
        *ang = 2*M_PI + *ang;

    return GSL_SUCCESS;
}


inline double area_between_vectors_triangle(const gsl_vector *a, const gsl_vector *b)
{
    double a1 = gsl_vector_get(a,0);
    double a2 = gsl_vector_get(a,1);
    double b1 = gsl_vector_get(b,0);
    double b2 = gsl_vector_get(b,1);
    double area = (a1*b2 - a2*b1)/2;
    area = fabs(area);

    return area;
}


inline double area_segment(double ang, double r)
{
    return ang/2 * gsl_pow_2(r); 
}


inline int point_from_alpha(double alpha, double norm, gsl_vector *center, gsl_vector *v)
{
    gsl_vector_set(v,0,-gsl_sf_cos(alpha));
    gsl_vector_set(v,1,gsl_sf_sin(alpha));
    gsl_vector_scale(v,norm);
    gsl_vector_add(v,center);

    return GSL_SUCCESS;
}


void print_matrix(FILE *stream, const gsl_matrix *m)
{
    for(int i = 0; i < m->size1; ++i)
    {
        for(int j = 0; j < m->size2; ++j)
            fprintf(stream,"%3.3f ",gsl_matrix_get(m,i,j));
        fprintf(stream,"\n");
    }
}