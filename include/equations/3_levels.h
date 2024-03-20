#ifndef _EQUATIONS_3_LEVELS_H
#define _EQUATIONS_3_LEVELS_H

#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>


struct system_3_levels_params
{
    double phi_ad_0;
    double phi_cb_0;
    double phi_dc_0;
    double phi_df_0;
    double phi_ec_0;
    double phi_fe_0;
    double phi_ge_0;
    double phi_gf_0;

    double r_top_0,r_mid_0,r_bot_0;
    double p_top_0,p_mid_0,p_bot_0;

    double Ax, Ay;
    double Bx, By;

    double p,k;
    double p_ac;
};


int system_3_levels_f(const gsl_vector *x, void *p, gsl_vector *f);
int system_3_levels_df(const gsl_vector *x, void *p, gsl_matrix *J);
int system_3_levels_fdf(const gsl_vector *x, void *p, gsl_vector *f, gsl_matrix *J);
int system_3_levels_eval();

#endif // _EQUATIONS_3_LEVELS_H