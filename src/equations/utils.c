#include <equations/utils.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_sf.h>
#include <math.h>
#include <stdlib.h>


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


void J_estimate(const gsl_vector *x0, void *params, gsl_matrix *J_est, int (*f)(const gsl_vector *x, void *params, gsl_vector *f))
{
    size_t dim = x0->size;
    gsl_vector *f0 = gsl_vector_alloc(dim);
    f(x0,params,f0);
    gsl_vector *x0_new = gsl_vector_alloc(dim);
    gsl_vector_memcpy(x0_new,x0);
    const double eps = 1e-9;
    for(int i = 0; i < dim; ++i)
    {
        gsl_vector df_dxi = gsl_matrix_column(J_est,i).vector;
        double old_xi_v = gsl_vector_get(x0,i);
        gsl_vector_set(x0_new,i,old_xi_v+eps);
        f(x0_new,params,&df_dxi);
        gsl_vector_sub(&df_dxi,f0);
        gsl_vector_scale(&df_dxi,1.0/eps);
        gsl_vector_set(x0_new,i,old_xi_v);
    }

    gsl_vector_free(f0);
    gsl_vector_free(x0_new);
}


inline void print_J_diff(FILE *stream, const gsl_vector *x, void *params, int (*f)(const gsl_vector *x, void *params, gsl_vector *f), int (*df)(const gsl_vector *x, void *params, gsl_matrix *df))
{
    size_t dim = x->size;
    gsl_matrix *J_est = gsl_matrix_alloc(dim,dim);
    gsl_matrix *J_diff = gsl_matrix_alloc(dim,dim);
    gsl_matrix *J = gsl_matrix_alloc(dim,dim);
    df(x,params,J);
    J_estimate(x,params,J_est,f);
    gsl_matrix_memcpy(J_diff,J);
    gsl_matrix_sub(J_diff,J_est);
    fprintf(stream,"Jacobian true:\n");
    print_matrix(stream,J);
    fprintf(stream,"Jacobian est:\n");
    print_matrix(stream,J_est);
    fprintf(stream,"Jacobian diff:\n");
    print_matrix(stream,J_diff);

    fprintf(stream,"Errors:\n");
    for(size_t row_i = 0; row_i < dim; ++row_i)
    {
        for(size_t col_i = 0; col_i < dim; ++col_i)
        {
            double df_real = gsl_matrix_get(J,row_i,col_i);
            double df_est = gsl_matrix_get(J_est,row_i,col_i);
            double diff = gsl_matrix_get(J_diff,row_i,col_i);
            double err = fabs(diff/(MIN(df_real,df_est)+1e-9));
            if(err > 0.05)
                fprintf(stream,"J(%zd, %zd): real df: %f; estimated df: %f\n",row_i,col_i,df_real,df_est);
        }
    }

    gsl_matrix_free(J_est);
    gsl_matrix_free(J_diff);
    gsl_matrix_free(J);
}