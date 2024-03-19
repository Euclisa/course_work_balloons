#ifndef _EQUATIONS_2_LEVELS
#define _EQUATIONS_2_LEVELS

#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>


extern const size_t N_eq;

struct system_2_levels_params
{
    double phi_ad_0;
    double phi_cb_0;
    double phi_dc_0;
    double phi_ec_0;
    double phi_ed_0;

    double r_top_0,r_bot_0;
    double S_top_0,S_bot_0;
    double p_top_0,p_bot_0;

    double Ax, Ay;
    double Bx, By;

    double p,k;
    double p_ac;
};


int system_2_levels_f(const gsl_vector *x, void *p, gsl_vector *f);
int system_2_levels_df(const gsl_vector *x, void *p, gsl_matrix *J);
int system_2_levels_fdf(const gsl_vector *x, void *p, gsl_vector *f, gsl_matrix *J);
int system_2_levels_eval();

#endif