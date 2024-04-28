#ifndef _EQUATIONS_3_LEVELS_H
#define _EQUATIONS_3_LEVELS_H

#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>


struct system_3_levels_user_params
{
    double phi_ad_0;
    double phi_dc_0;
    double phi_df_0;
    double phi_fe_0;

    double r_top_0,r_mid_0,r_bot_0;
    double p_top_0,p_mid_0,p_bot_0;

    double Ax, Ay;
    double Bx, By;

    double p_atm;
    double p_ac;
};

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

    double p_atm,k;
    double p_ac;
};

struct system_3_levels_result
{
    double phi_ad, r_ad, x_ad, y_ad, a_ad;
    double phi_cb, r_cb, x_cb, y_cb, a_cb;
    double phi_dc, r_dc, x_dc, y_dc, a_dc;
    double phi_df, r_df, x_df, y_df, a_df;
    double phi_ec, r_ec, x_ec, y_ec, a_ec;
    double phi_fe, r_fe, x_fe, y_fe, a_fe;
    double phi_ge, r_ge, y_ge;
    double phi_gf, r_gf, y_gf;

    double x_bot;

    double p_bot, p_mid, p_top;
};


int system_3_levels_f(const gsl_vector *x, void *p, gsl_vector *f);
int system_3_levels_df(const gsl_vector *x, void *p, gsl_matrix *J);
int system_3_levels_fdf(const gsl_vector *x, void *p, gsl_vector *f, gsl_matrix *J);
int system_3_levels_eval_f();
int system_3_levels_eval(const struct system_3_levels_user_params *user_params, struct system_3_levels_result *result);

#endif // _EQUATIONS_3_LEVELS_H