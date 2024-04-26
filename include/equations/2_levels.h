#ifndef _EQUATIONS_2_LEVELS_H
#define _EQUATIONS_2_LEVELS_H

#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>


struct system_2_levels_user_params
{
    double phi_ad_0;
    double phi_dc_0;

    double r_top_0,r_bot_0;
    double p_top_0,p_bot_0;

    double Ax, Ay;
    double Bx, By;

    double p_ac, p_atm;

    double k;
};


struct system_2_levels_params
{
    double phi_ad_0;
    double phi_cb_0;
    double phi_dc_0;
    double phi_ec_0;
    double phi_ed_0;

    double r_top_0,r_bot_0;
    double p_top_0,p_bot_0;

    double Ax, Ay;
    double Bx, By;

    double p_ac, p_atm;

    double S_top_0, S_bot_0;

    double k;
};


struct system_2_levels_result
{
    double phi_ad, r_ad, x_ad, y_ad, a_ad;
    double phi_cb, r_cb, x_cb, y_cb, a_cb;
    double phi_dc, r_dc, x_dc, y_dc, a_dc;
    double phi_ed, r_ed, y_ed;
    double phi_ec, r_ec, y_ec;

    double x_bot;

    double p_bot, p_top;
};


int system_2_levels_f(const gsl_vector *x, void *p, gsl_vector *f);
int system_2_levels_df(const gsl_vector *x, void *p, gsl_matrix *J);
int system_2_levels_fdf(const gsl_vector *x, void *p, gsl_vector *f, gsl_matrix *J);
int system_2_levels_eval_f();
int system_2_levels_eval(const struct system_2_levels_user_params *user_params, struct system_2_levels_result *result);
int system_2_levels_adiabatic_eval(const struct system_2_levels_user_params *user_params, struct system_2_levels_result *result);

#endif // _EQUATIONS_2_LEVELS_H