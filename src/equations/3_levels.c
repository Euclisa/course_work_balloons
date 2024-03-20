#include <equations/3_levels.h>
#include <equations/utils.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>
#include <gsl/gsl_sf_trig.h>
#include <gsl/gsl_blas.h>

static const size_t N_eq = 40;


int system_3_levels_f(const gsl_vector *x, void *p, gsl_vector *f)
{
    struct system_3_levels_params *params = (struct system_3_levels_params*)p;
    
    const double phi_ad = gsl_vector_get(x,0);
    const double r_ad = gsl_vector_get(x,1);
    const double x_ad = gsl_vector_get(x,2);
    const double y_ad = gsl_vector_get(x,3);
    const double a_ad = gsl_vector_get(x,4);
    const double phi_cb = gsl_vector_get(x,5);
    const double r_cb = gsl_vector_get(x,6);
    const double x_cb = gsl_vector_get(x,7);
    const double y_cb = gsl_vector_get(x,8);
    const double a_cb = gsl_vector_get(x,9);
    const double phi_dc = gsl_vector_get(x,10);
    const double r_dc = gsl_vector_get(x,11);
    const double x_dc = gsl_vector_get(x,12);
    const double y_dc = gsl_vector_get(x,13);
    const double a_dc = gsl_vector_get(x,14);
    const double phi_df = gsl_vector_get(x,15);
    const double r_df = gsl_vector_get(x,16);
    const double x_df = gsl_vector_get(x,17);
    const double y_df = gsl_vector_get(x,18);
    const double a_df = gsl_vector_get(x,19);
    const double phi_ec = gsl_vector_get(x,20);
    const double r_ec = gsl_vector_get(x,21);
    const double x_ec = gsl_vector_get(x,22);
    const double y_ec = gsl_vector_get(x,23);
    const double a_ec = gsl_vector_get(x,24);
    const double phi_fe = gsl_vector_get(x,25);
    const double r_fe = gsl_vector_get(x,26);
    const double x_fe = gsl_vector_get(x,27);
    const double y_fe = gsl_vector_get(x,28);
    const double a_fe = gsl_vector_get(x,29);
    const double phi_ge = gsl_vector_get(x,30);
    const double r_ge = gsl_vector_get(x,31);
    const double y_ge = gsl_vector_get(x,32);
    const double phi_gf = gsl_vector_get(x,33);
    const double r_gf = gsl_vector_get(x,34);
    const double y_gf = gsl_vector_get(x,35);
    const double x_bot = gsl_vector_get(x,36);
    const double p_top = gsl_vector_get(x,37);
    const double p_mid = gsl_vector_get(x,38);
    const double p_bot = gsl_vector_get(x,39);

    const double a_ge = 3*M_PI_2;
    const double a_gf = 3*M_PI_2;

    // Preservation of length
    double eq1 = r_ad*phi_ad - params->r_top_0*params->phi_ad_0;
    gsl_vector_set(f,0,eq1);
    double eq2 = r_cb*phi_cb - params->r_top_0*params->phi_cb_0;
    gsl_vector_set(f,1,eq2);
    double eq3 = r_dc*phi_dc - params->r_top_0*params->phi_dc_0;
    gsl_vector_set(f,2,eq3);
    double eq4 = r_df*phi_df - params->r_mid_0*params->phi_df_0;
    gsl_vector_set(f,3,eq4);
    double eq5 = r_ec*phi_ec - params->r_mid_0*params->phi_ec_0;
    gsl_vector_set(f,4,eq5);
    double eq6 = r_fe*phi_fe - params->r_mid_0*params->phi_fe_0;
    gsl_vector_set(f,5,eq6);
    double eq7 = r_ge*phi_ge + r_gf*phi_gf - params->r_bot_0*(params->phi_ge_0 + params->phi_gf_0);
    gsl_vector_set(f,6,eq7);

    // Point A continuity
    double eq8 = -r_ad*gsl_sf_cos(a_ad) + x_ad - params->Ax;
    gsl_vector_set(f,7,eq8);
    double eq9 = r_ad*gsl_sf_sin(a_ad) + y_ad - params->Ay;
    gsl_vector_set(f,8,eq9);

    // Point B continuity
    double eq10 = -r_cb*gsl_sf_cos(a_cb+phi_cb) + x_cb - params->Bx;
    gsl_vector_set(f,9,eq10);
    double eq11 = r_cb*gsl_sf_sin(a_cb+phi_cb) + y_cb - params->By;
    gsl_vector_set(f,10,eq11);

    // Point C continuity
    double eq12 = -r_dc*gsl_sf_cos(a_dc+phi_dc) + x_dc - (-r_cb*gsl_sf_cos(a_cb) + x_cb);
    gsl_vector_set(f,11,eq12);
    double eq13 = r_dc*gsl_sf_sin(a_dc+phi_dc) + y_dc - (r_cb*gsl_sf_sin(a_cb) + y_cb);
    gsl_vector_set(f,12,eq13);
    double eq14 = -r_dc*gsl_sf_cos(a_dc+phi_dc) + x_dc - (-r_ec*gsl_sf_cos(a_ec+phi_ec) + x_ec);
    gsl_vector_set(f,13,eq14);
    double eq15 = r_dc*gsl_sf_sin(a_dc+phi_dc) + y_dc - (r_ec*gsl_sf_sin(a_ec+phi_ec) + y_ec);
    gsl_vector_set(f,14,eq15);

    // Point D continuity
    double eq16 = -r_ad*gsl_sf_cos(a_ad+phi_ad) + x_ad - (-r_dc*gsl_sf_cos(a_dc) + x_dc);
    gsl_vector_set(f,15,eq16);
    double eq17 = r_ad*gsl_sf_sin(a_ad+phi_ad) + y_ad - (r_dc*gsl_sf_sin(a_dc) + y_dc);
    gsl_vector_set(f,16,eq17);
    double eq18 = -r_ad*gsl_sf_cos(a_ad+phi_ad) + x_ad - (-r_df*gsl_sf_cos(a_df) + x_df);
    gsl_vector_set(f,17,eq18);
    double eq19 = r_ad*gsl_sf_sin(a_ad+phi_ad) + y_ad - (r_df*gsl_sf_sin(a_df) + y_df);
    gsl_vector_set(f,18,eq19);

    // Point E continuity
    double eq20 = -r_ec*gsl_sf_cos(a_ec) + x_ec - (-r_fe*gsl_sf_cos(a_fe+phi_fe) + x_fe);
    gsl_vector_set(f,19,eq20);
    double eq21 = r_ec*gsl_sf_sin(a_ec) + y_ec - (r_fe*gsl_sf_sin(a_fe+phi_fe) + y_fe);
    gsl_vector_set(f,20,eq21);
    double eq22 = -r_ec*gsl_sf_cos(a_ec) + x_ec - (-r_ge*gsl_sf_cos(a_ge+phi_ge) + x_bot);
    gsl_vector_set(f,21,eq22);
    double eq23 = r_ec*gsl_sf_sin(a_ec) + y_ec - (r_ge*gsl_sf_sin(a_ge+phi_ge) + y_ge);
    gsl_vector_set(f,22,eq23);

    // Point F continuity
    double eq24 = -r_fe*gsl_sf_cos(a_fe) + x_fe - (-r_df*gsl_sf_cos(a_df+phi_df) + x_df);
    gsl_vector_set(f,23,eq24);
    double eq25 = r_fe*gsl_sf_sin(a_fe) + y_fe - (r_df*gsl_sf_sin(a_df+phi_df) + y_df);
    gsl_vector_set(f,24,eq25);
    double eq26 = -r_fe*gsl_sf_cos(a_fe) + x_fe - (-r_gf*gsl_sf_cos(a_gf-phi_gf) + x_bot);
    gsl_vector_set(f,25,eq26);
    double eq27 = r_fe*gsl_sf_sin(a_fe) + y_fe - (r_gf*gsl_sf_sin(a_gf-phi_gf) + y_gf);
    gsl_vector_set(f,26,eq27);

    // Point G continuity
    double eq28 = y_ge - r_ge - (y_gf - r_gf);
    gsl_vector_set(f,27,eq28);

    // Point D steadiness
    double eq29 = -p_top*r_ad*gsl_sf_sin(a_ad+phi_ad) + (p_top-p_mid)*r_dc*gsl_sf_sin(a_dc) + p_mid*r_df*gsl_sf_sin(a_df);
    gsl_vector_set(f,28,eq29);
    double eq30 = -p_top*r_ad*gsl_sf_cos(a_ad+phi_ad) + (p_top-p_mid)*r_dc*gsl_sf_cos(a_dc) + p_mid*r_df*gsl_sf_cos(a_df);
    gsl_vector_set(f,29,eq30);

    // Point C steadiness
    double eq31 = (p_top-params->p_ac)*r_cb*gsl_sf_sin(a_cb) - (p_top-p_mid)*r_dc*gsl_sf_sin(a_dc+phi_dc) - (p_mid-params->p_ac)*r_ec*gsl_sf_sin(a_ec+phi_ec);
    gsl_vector_set(f,30,eq31);
    double eq32 = (p_top-params->p_ac)*r_cb*gsl_sf_cos(a_cb) - (p_top-p_mid)*r_dc*gsl_sf_cos(a_dc+phi_dc) - (p_mid-params->p_ac)*r_ec*gsl_sf_cos(a_ec+phi_ec);
    gsl_vector_set(f,31,eq32);

    // Point E steadiness
    double eq33 = (p_mid-params->p_ac)*r_ec*gsl_sf_sin(a_ec) - (p_mid-p_bot)*r_fe*gsl_sf_sin(a_fe+phi_fe) - (p_bot-params->p_ac)*r_ge*gsl_sf_sin(a_ge+phi_ge);
    gsl_vector_set(f,32,eq33);
    double eq34 = (p_mid-params->p_ac)*r_ec*gsl_sf_cos(a_ec) - (p_mid-p_bot)*r_fe*gsl_sf_cos(a_fe+phi_fe) - (p_bot-params->p_ac)*r_ge*gsl_sf_cos(a_ge+phi_ge);
    gsl_vector_set(f,33,eq34);

    // Point F steadiness
    double eq35 = -p_mid*r_df*gsl_sf_sin(a_df+phi_df) + (p_mid-p_bot)*r_fe*gsl_sf_sin(a_fe) + p_bot*r_gf*gsl_sf_sin(a_gf-phi_gf);
    gsl_vector_set(f,34,eq35);
    double eq36 = -p_mid*r_df*gsl_sf_cos(a_df+phi_df) + (p_mid-p_bot)*r_fe*gsl_sf_cos(a_fe) + p_bot*r_gf*gsl_sf_cos(a_gf-phi_gf);
    gsl_vector_set(f,35,eq36);

    // Point G steadiness
    double eq37 = (p_bot-params->p_ac)*r_ge - p_bot*r_gf;
    gsl_vector_set(f,36,eq37);


    // Balloons pressures
    double eq38 = p_top - params->p_top_0;
    gsl_vector_set(f,37,eq38);
    double eq39 = p_mid - params->p_mid_0;
    gsl_vector_set(f,38,eq39);
    double eq40 = p_bot - params->p_bot_0;
    gsl_vector_set(f,39,eq40);

    return GSL_SUCCESS;
}


int system_3_levels_df(const gsl_vector *x, void *p, gsl_matrix *J)
{
    struct system_3_levels_params *params = (struct system_3_levels_params*)p;
    
    const double phi_ad = gsl_vector_get(x,0);
    const double r_ad = gsl_vector_get(x,1);
    //const double x_ad = gsl_vector_get(x,2);
    //const double y_ad = gsl_vector_get(x,3);
    const double a_ad = gsl_vector_get(x,4);
    const double phi_cb = gsl_vector_get(x,5);
    const double r_cb = gsl_vector_get(x,6);
    //const double x_cb = gsl_vector_get(x,7);
    //const double y_cb = gsl_vector_get(x,8);
    const double a_cb = gsl_vector_get(x,9);
    const double phi_dc = gsl_vector_get(x,10);
    const double r_dc = gsl_vector_get(x,11);
    //const double x_dc = gsl_vector_get(x,12);
    //const double y_dc = gsl_vector_get(x,13);
    const double a_dc = gsl_vector_get(x,14);
    const double phi_df = gsl_vector_get(x,15);
    const double r_df = gsl_vector_get(x,16);
    //const double x_df = gsl_vector_get(x,17);
    //const double y_df = gsl_vector_get(x,18);
    const double a_df = gsl_vector_get(x,19);
    const double phi_ec = gsl_vector_get(x,20);
    const double r_ec = gsl_vector_get(x,21);
    //const double x_ec = gsl_vector_get(x,22);
    //const double y_ec = gsl_vector_get(x,23);
    const double a_ec = gsl_vector_get(x,24);
    const double phi_fe = gsl_vector_get(x,25);
    const double r_fe = gsl_vector_get(x,26);
    //const double x_fe = gsl_vector_get(x,27);
    //const double y_fe = gsl_vector_get(x,28);
    const double a_fe = gsl_vector_get(x,29);
    const double phi_ge = gsl_vector_get(x,30);
    const double r_ge = gsl_vector_get(x,31);
    //const double y_ge = gsl_vector_get(x,32);
    const double phi_gf = gsl_vector_get(x,33);
    const double r_gf = gsl_vector_get(x,34);
    //const double y_gf = gsl_vector_get(x,35);
    //const double x_bot = gsl_vector_get(x,36);
    const double p_top = gsl_vector_get(x,37);
    const double p_mid = gsl_vector_get(x,38);
    const double p_bot = gsl_vector_get(x,39);

    const double a_ge = 3*M_PI_2;
    const double a_gf = 3*M_PI_2;

    gsl_matrix_set_all(J,0);

    // Preservation of length
    gsl_matrix_set(J,0,0,r_ad);
    gsl_matrix_set(J,0,1,phi_ad);
    
    gsl_matrix_set(J,1,5,r_cb);
    gsl_matrix_set(J,1,6,phi_cb);

    gsl_matrix_set(J,2,10,r_dc);
    gsl_matrix_set(J,2,11,phi_dc);
    
    gsl_matrix_set(J,3,15,r_df);
    gsl_matrix_set(J,3,16,phi_df);

    gsl_matrix_set(J,4,20,r_ec);
    gsl_matrix_set(J,4,21,phi_ec);

    gsl_matrix_set(J,5,25,r_fe);
    gsl_matrix_set(J,5,26,phi_fe);

    gsl_matrix_set(J,6,30,r_ge);
    gsl_matrix_set(J,6,31,phi_ge);
    gsl_matrix_set(J,6,33,r_gf);
    gsl_matrix_set(J,6,34,phi_gf);


    // Point A continuity
    gsl_matrix_set(J,7,1,-gsl_sf_cos(a_ad));
    gsl_matrix_set(J,7,2,1);
    gsl_matrix_set(J,7,4,r_ad*gsl_sf_sin(a_ad));
    
    gsl_matrix_set(J,8,1,gsl_sf_sin(a_ad));
    gsl_matrix_set(J,8,3,1);
    gsl_matrix_set(J,8,4,r_ad*gsl_sf_cos(a_ad));
    

    // Point B continuity
    gsl_matrix_set(J,9,5,r_cb*gsl_sf_sin(phi_cb+a_cb));
    gsl_matrix_set(J,9,6,-gsl_sf_cos(phi_cb+a_cb));
    gsl_matrix_set(J,9,7,1);
    gsl_matrix_set(J,9,9,r_cb*gsl_sf_sin(phi_cb+a_cb));
    
    gsl_matrix_set(J,10,5,r_cb*gsl_sf_cos(phi_cb+a_cb));
    gsl_matrix_set(J,10,6,gsl_sf_sin(phi_cb+a_cb));
    gsl_matrix_set(J,10,8,1);
    gsl_matrix_set(J,10,9,r_cb*gsl_sf_cos(phi_cb+a_cb));
    

    // Point C continuity
    gsl_matrix_set(J,11,6,gsl_sf_cos(a_cb));
    gsl_matrix_set(J,11,7,-1);
    gsl_matrix_set(J,11,9,-r_cb*gsl_sf_sin(a_cb));
    gsl_matrix_set(J,11,10,r_dc*gsl_sf_sin(phi_dc+a_dc));
    gsl_matrix_set(J,11,11,-gsl_sf_cos(phi_dc+a_dc));
    gsl_matrix_set(J,11,12,1);
    gsl_matrix_set(J,11,14,r_dc*gsl_sf_sin(phi_dc+a_dc));

    gsl_matrix_set(J,12,6,-gsl_sf_sin(a_cb));
    gsl_matrix_set(J,12,8,-1);
    gsl_matrix_set(J,12,9,-r_cb*gsl_sf_cos(a_cb));
    gsl_matrix_set(J,12,10,r_dc*gsl_sf_cos(phi_dc+a_dc));
    gsl_matrix_set(J,12,11,gsl_sf_sin(phi_dc+a_dc));
    gsl_matrix_set(J,12,13,1);
    gsl_matrix_set(J,12,14,r_dc*gsl_sf_cos(phi_dc+a_dc));
    
    gsl_matrix_set(J,13,10,r_dc*gsl_sf_sin(phi_dc+a_dc));
    gsl_matrix_set(J,13,11,-gsl_sf_cos(phi_dc+a_dc));
    gsl_matrix_set(J,13,12,1);
    gsl_matrix_set(J,13,14,r_dc*gsl_sf_sin(phi_dc+a_dc));
    gsl_matrix_set(J,13,20,-r_ec*gsl_sf_sin(a_ec+phi_ec));
    gsl_matrix_set(J,13,21,gsl_sf_cos(a_ec+phi_ec));
    gsl_matrix_set(J,13,22,-1);
    gsl_matrix_set(J,13,24,-r_ec*gsl_sf_sin(a_ec+phi_ec));

    gsl_matrix_set(J,14,10,r_dc*gsl_sf_cos(phi_dc+a_dc));
    gsl_matrix_set(J,14,11,gsl_sf_sin(phi_dc+a_dc));
    gsl_matrix_set(J,14,13,1);
    gsl_matrix_set(J,14,14,r_dc*gsl_sf_cos(phi_dc+a_dc));
    gsl_matrix_set(J,14,20,-r_ec*gsl_sf_cos(a_ec+phi_ec));
    gsl_matrix_set(J,14,21,-gsl_sf_sin(a_ec+phi_ec));
    gsl_matrix_set(J,14,23,-1);
    gsl_matrix_set(J,14,24,-r_ec*gsl_sf_cos(a_ec+phi_ec));


    // Point D continuity
    gsl_matrix_set(J,15,0,r_ad*gsl_sf_sin(phi_ad+a_ad));
    gsl_matrix_set(J,15,1,-gsl_sf_cos(phi_ad+a_ad));
    gsl_matrix_set(J,15,2,1);
    gsl_matrix_set(J,15,4,r_ad*gsl_sf_sin(phi_ad+a_ad));
    gsl_matrix_set(J,15,11,gsl_sf_cos(a_dc));
    gsl_matrix_set(J,15,12,-1);
    gsl_matrix_set(J,15,14,-r_dc*gsl_sf_sin(a_dc));

    gsl_matrix_set(J,16,0,r_ad*gsl_sf_cos(phi_ad+a_ad));
    gsl_matrix_set(J,16,1,gsl_sf_sin(phi_ad+a_ad));
    gsl_matrix_set(J,16,3,1);
    gsl_matrix_set(J,16,4,r_ad*gsl_sf_cos(phi_ad+a_ad));
    gsl_matrix_set(J,16,11,-gsl_sf_sin(a_dc));
    gsl_matrix_set(J,16,13,-1);
    gsl_matrix_set(J,16,14,-r_dc*gsl_sf_cos(a_dc));

    gsl_matrix_set(J,17,0,r_ad*gsl_sf_sin(phi_ad+a_ad));
    gsl_matrix_set(J,17,1,-gsl_sf_cos(phi_ad+a_ad));
    gsl_matrix_set(J,17,2,1);
    gsl_matrix_set(J,17,4,r_ad*gsl_sf_sin(phi_ad+a_ad));
    gsl_matrix_set(J,17,16,gsl_sf_cos(a_df));
    gsl_matrix_set(J,17,17,-1);
    gsl_matrix_set(J,17,19,-r_df*gsl_sf_sin(a_df));

    gsl_matrix_set(J,18,0,r_ad*gsl_sf_cos(phi_ad+a_ad));
    gsl_matrix_set(J,18,1,gsl_sf_sin(phi_ad+a_ad));
    gsl_matrix_set(J,18,3,1);
    gsl_matrix_set(J,18,4,r_ad*gsl_sf_cos(phi_ad+a_ad));
    gsl_matrix_set(J,18,16,-gsl_sf_sin(a_df));
    gsl_matrix_set(J,18,18,-1);
    gsl_matrix_set(J,18,19,-r_df*gsl_sf_cos(a_df));


    // Point E continuity
    gsl_matrix_set(J,19,21,-gsl_sf_cos(a_ec));
    gsl_matrix_set(J,19,22,1);
    gsl_matrix_set(J,19,24,r_ec*gsl_sf_sin(a_ec));
    gsl_matrix_set(J,19,25,-r_fe*gsl_sf_sin(a_fe+phi_fe));
    gsl_matrix_set(J,19,26,gsl_sf_cos(a_fe+phi_fe));
    gsl_matrix_set(J,19,27,-1);
    gsl_matrix_set(J,19,29,-r_fe*gsl_sf_sin(a_fe+phi_fe));

    gsl_matrix_set(J,20,21,gsl_sf_sin(a_ec));
    gsl_matrix_set(J,20,23,1);
    gsl_matrix_set(J,20,24,r_ec*gsl_sf_cos(a_ec));
    gsl_matrix_set(J,20,25,-r_fe*gsl_sf_cos(a_fe+phi_fe));
    gsl_matrix_set(J,20,26,-gsl_sf_sin(a_fe+phi_fe));
    gsl_matrix_set(J,20,28,-1);
    gsl_matrix_set(J,20,29,-r_fe*gsl_sf_cos(a_fe+phi_fe));

    gsl_matrix_set(J,21,21,-gsl_sf_cos(a_ec));
    gsl_matrix_set(J,21,22,1);
    gsl_matrix_set(J,21,24,r_ec*gsl_sf_sin(a_ec));
    gsl_matrix_set(J,21,30,-r_ge*gsl_sf_sin(a_ge+phi_ge));
    gsl_matrix_set(J,21,31,gsl_sf_cos(a_ge+phi_ge));
    gsl_matrix_set(J,21,36,-1);

    gsl_matrix_set(J,22,21,gsl_sf_sin(a_ec));
    gsl_matrix_set(J,22,23,1);
    gsl_matrix_set(J,22,24,r_ec*gsl_sf_cos(a_ec));
    gsl_matrix_set(J,22,30,-r_ge*gsl_sf_cos(a_ge+phi_ge));
    gsl_matrix_set(J,22,31,-gsl_sf_sin(a_ge+phi_ge));
    gsl_matrix_set(J,22,32,-1);


    // Point F continuity
    gsl_matrix_set(J,23,26,-gsl_sf_cos(a_fe));
    gsl_matrix_set(J,23,27,1);
    gsl_matrix_set(J,23,29,r_fe*gsl_sf_sin(a_fe));
    gsl_matrix_set(J,23,15,-r_df*gsl_sf_sin(a_df+phi_df));
    gsl_matrix_set(J,23,16,gsl_sf_cos(a_df+phi_df));
    gsl_matrix_set(J,23,17,-1);
    gsl_matrix_set(J,23,19,-r_df*gsl_sf_sin(a_df+phi_df));

    gsl_matrix_set(J,24,26,gsl_sf_sin(a_fe));
    gsl_matrix_set(J,24,28,1);
    gsl_matrix_set(J,24,29,r_fe*gsl_sf_cos(a_fe));
    gsl_matrix_set(J,24,15,-r_df*gsl_sf_cos(a_df+phi_df));
    gsl_matrix_set(J,24,16,-gsl_sf_sin(a_df+phi_df));
    gsl_matrix_set(J,24,18,-1);
    gsl_matrix_set(J,24,19,-r_df*gsl_sf_cos(a_df+phi_df));

    gsl_matrix_set(J,25,26,-gsl_sf_cos(a_fe));
    gsl_matrix_set(J,25,27,1);
    gsl_matrix_set(J,25,29,r_fe*gsl_sf_sin(a_fe));
    gsl_matrix_set(J,25,33,r_gf*gsl_sf_sin(a_gf-phi_gf));
    gsl_matrix_set(J,25,34,gsl_sf_cos(a_gf-phi_gf));
    gsl_matrix_set(J,25,36,-1);

    gsl_matrix_set(J,26,26,gsl_sf_sin(a_fe));
    gsl_matrix_set(J,26,28,1);
    gsl_matrix_set(J,26,29,r_fe*gsl_sf_cos(a_fe));
    gsl_matrix_set(J,26,33,r_gf*gsl_sf_cos(a_gf-phi_gf));
    gsl_matrix_set(J,26,34,-gsl_sf_sin(a_gf-phi_gf));
    gsl_matrix_set(J,26,35,-1);

    // Point G continuity
    gsl_matrix_set(J,27,31,-1);
    gsl_matrix_set(J,27,32,1);
    gsl_matrix_set(J,27,34,1);
    gsl_matrix_set(J,27,35,-1);


    // Point D steadiness
    gsl_matrix_set(J,28,0,-p_top*r_ad*gsl_sf_cos(a_ad+phi_ad));
    gsl_matrix_set(J,28,1,-p_top*gsl_sf_sin(a_ad+phi_ad));
    gsl_matrix_set(J,28,4,-p_top*r_ad*gsl_sf_cos(a_ad+phi_ad));
    gsl_matrix_set(J,28,11,(p_top-p_mid)*gsl_sf_sin(a_dc));
    gsl_matrix_set(J,28,14,(p_top-p_mid)*r_dc*gsl_sf_cos(a_dc));
    gsl_matrix_set(J,28,16,p_mid*gsl_sf_sin(a_df));
    gsl_matrix_set(J,28,19,p_mid*r_df*gsl_sf_cos(a_df));
    gsl_matrix_set(J,28,37,-r_ad*gsl_sf_sin(a_ad+phi_ad) + r_dc*gsl_sf_sin(a_dc));
    gsl_matrix_set(J,28,38,-r_dc*gsl_sf_sin(a_dc) + r_df*gsl_sf_sin(a_df));

    gsl_matrix_set(J,29,0,p_top*r_ad*gsl_sf_sin(a_ad+phi_ad));
    gsl_matrix_set(J,29,1,-p_top*gsl_sf_cos(a_ad+phi_ad));
    gsl_matrix_set(J,29,4,p_top*r_ad*gsl_sf_sin(a_ad+phi_ad));
    gsl_matrix_set(J,29,11,(p_top-p_mid)*gsl_sf_cos(a_dc));
    gsl_matrix_set(J,29,14,-(p_top-p_mid)*r_dc*gsl_sf_sin(a_dc));
    gsl_matrix_set(J,29,16,p_mid*gsl_sf_cos(a_df));
    gsl_matrix_set(J,29,19,-p_mid*r_df*gsl_sf_sin(a_df));
    gsl_matrix_set(J,29,37,-r_ad*gsl_sf_cos(a_ad+phi_ad) + r_dc*gsl_sf_cos(a_dc));
    gsl_matrix_set(J,29,38,-r_dc*gsl_sf_cos(a_dc) + r_df*gsl_sf_cos(a_df));


    // Point C steadiness
    gsl_matrix_set(J,30,6,(p_top-params->p_ac)*gsl_sf_sin(a_cb));
    gsl_matrix_set(J,30,9,(p_top-params->p_ac)*r_cb*gsl_sf_cos(a_cb));
    gsl_matrix_set(J,30,10,-(p_top-p_mid)*r_dc*gsl_sf_cos(a_dc+phi_dc));
    gsl_matrix_set(J,30,11,-(p_top-p_mid)*gsl_sf_sin(a_dc+phi_dc));
    gsl_matrix_set(J,30,14,-(p_top-p_mid)*r_dc*gsl_sf_cos(a_dc+phi_dc));
    gsl_matrix_set(J,30,20,-(p_mid-params->p_ac)*r_ec*gsl_sf_cos(a_ec+phi_ec));
    gsl_matrix_set(J,30,21,-(p_mid-params->p_ac)*gsl_sf_sin(a_ec+phi_ec));
    gsl_matrix_set(J,30,24,-(p_mid-params->p_ac)*r_ec*gsl_sf_cos(a_ec+phi_ec));
    gsl_matrix_set(J,30,37,r_cb*gsl_sf_sin(a_cb) - r_dc*gsl_sf_sin(a_dc+phi_dc));
    gsl_matrix_set(J,30,38,r_dc*gsl_sf_sin(a_dc+phi_dc) - r_ec*gsl_sf_sin(a_ec+phi_ec));

    gsl_matrix_set(J,31,6,(p_top-params->p_ac)*gsl_sf_cos(a_cb));
    gsl_matrix_set(J,31,9,-(p_top-params->p_ac)*r_cb*gsl_sf_sin(a_cb));
    gsl_matrix_set(J,31,10,(p_top-p_mid)*r_dc*gsl_sf_sin(a_dc+phi_dc));
    gsl_matrix_set(J,31,11,-(p_top-p_mid)*gsl_sf_cos(a_dc+phi_dc));
    gsl_matrix_set(J,31,14,(p_top-p_mid)*r_dc*gsl_sf_sin(a_dc+phi_dc));
    gsl_matrix_set(J,31,20,(p_mid-params->p_ac)*r_ec*gsl_sf_sin(a_ec+phi_ec));
    gsl_matrix_set(J,31,21,-(p_mid-params->p_ac)*gsl_sf_cos(a_ec+phi_ec));
    gsl_matrix_set(J,31,24,(p_mid-params->p_ac)*r_ec*gsl_sf_sin(a_ec+phi_ec));
    gsl_matrix_set(J,31,37,r_cb*gsl_sf_cos(a_cb) - r_dc*gsl_sf_cos(a_dc+phi_dc));
    gsl_matrix_set(J,31,38,r_dc*gsl_sf_cos(a_dc+phi_dc) - r_ec*gsl_sf_cos(a_ec+phi_ec));


    // Point E steadiness
    gsl_matrix_set(J,32,21,(p_mid-params->p_ac)*gsl_sf_sin(a_ec));
    gsl_matrix_set(J,32,24,(p_mid-params->p_ac)*r_ec*gsl_sf_cos(a_ec));
    gsl_matrix_set(J,32,25,-(p_mid-p_bot)*r_fe*gsl_sf_cos(a_fe+phi_fe));
    gsl_matrix_set(J,32,26,-(p_mid-p_bot)*gsl_sf_sin(a_fe+phi_fe));
    gsl_matrix_set(J,32,29,-(p_mid-p_bot)*r_fe*gsl_sf_cos(a_fe+phi_fe));
    gsl_matrix_set(J,32,30,-(p_bot-params->p_ac)*r_ge*gsl_sf_cos(a_ge+phi_ge));
    gsl_matrix_set(J,32,31,-(p_bot-params->p_ac)*gsl_sf_sin(a_ge+phi_ge));
    gsl_matrix_set(J,32,38,r_ec*gsl_sf_sin(a_ec) - r_fe*gsl_sf_sin(a_fe+phi_fe));
    gsl_matrix_set(J,32,39,r_fe*gsl_sf_sin(a_fe+phi_fe) - r_ge*gsl_sf_sin(a_ge+phi_ge));

    gsl_matrix_set(J,33,21,(p_mid-params->p_ac)*gsl_sf_cos(a_ec));
    gsl_matrix_set(J,33,24,-(p_mid-params->p_ac)*r_ec*gsl_sf_sin(a_ec));
    gsl_matrix_set(J,33,25,(p_mid-p_bot)*r_fe*gsl_sf_sin(a_fe+phi_fe));
    gsl_matrix_set(J,33,26,-(p_mid-p_bot)*gsl_sf_cos(a_fe+phi_fe));
    gsl_matrix_set(J,33,29,(p_mid-p_bot)*r_fe*gsl_sf_sin(a_fe+phi_fe));
    gsl_matrix_set(J,33,30,(p_bot-params->p_ac)*r_ge*gsl_sf_sin(a_ge+phi_ge));
    gsl_matrix_set(J,33,31,-(p_bot-params->p_ac)*gsl_sf_cos(a_ge+phi_ge));
    gsl_matrix_set(J,33,38,r_ec*gsl_sf_cos(a_ec) - r_fe*gsl_sf_cos(a_fe+phi_fe));
    gsl_matrix_set(J,33,39,r_fe*gsl_sf_cos(a_fe+phi_fe) - r_ge*gsl_sf_cos(a_ge+phi_ge));


    // Point F steadiness
    gsl_matrix_set(J,34,15,-p_mid*r_df*gsl_sf_cos(a_df+phi_df));
    gsl_matrix_set(J,34,16,-p_mid*gsl_sf_sin(a_df+phi_df));
    gsl_matrix_set(J,34,19,-p_mid*r_df*gsl_sf_cos(a_df+phi_df));
    gsl_matrix_set(J,34,26,(p_mid-p_bot)*gsl_sf_sin(a_fe));
    gsl_matrix_set(J,34,29,(p_mid-p_bot)*r_fe*gsl_sf_cos(a_fe));
    gsl_matrix_set(J,34,33,-p_bot*r_gf*gsl_sf_cos(a_gf-phi_gf));
    gsl_matrix_set(J,34,34,p_bot*gsl_sf_sin(a_gf-phi_gf));
    gsl_matrix_set(J,34,38,-r_df*gsl_sf_sin(a_df+phi_df) + r_fe*gsl_sf_sin(a_fe));
    gsl_matrix_set(J,34,39,-r_fe*gsl_sf_sin(a_fe) + r_gf*gsl_sf_sin(a_gf-phi_gf));

    gsl_matrix_set(J,35,15,p_mid*r_df*gsl_sf_sin(a_df+phi_df));
    gsl_matrix_set(J,35,16,-p_mid*gsl_sf_cos(a_df+phi_df));
    gsl_matrix_set(J,35,19,p_mid*r_df*gsl_sf_sin(a_df+phi_df));
    gsl_matrix_set(J,35,26,(p_mid-p_bot)*gsl_sf_cos(a_fe));
    gsl_matrix_set(J,35,29,-(p_mid-p_bot)*r_fe*gsl_sf_sin(a_fe));
    gsl_matrix_set(J,35,33,p_bot*r_gf*gsl_sf_sin(a_gf-phi_gf));
    gsl_matrix_set(J,35,34,p_bot*gsl_sf_cos(a_gf-phi_gf));
    gsl_matrix_set(J,35,38,-r_df*gsl_sf_cos(a_df+phi_df) + r_fe*gsl_sf_cos(a_fe));
    gsl_matrix_set(J,35,39,-r_fe*gsl_sf_cos(a_fe) + r_gf*gsl_sf_cos(a_gf-phi_gf));


    // Point E steadiness
    gsl_matrix_set(J,36,34,-p_bot);
    gsl_matrix_set(J,36,31,(p_bot-params->p_ac));
    gsl_matrix_set(J,36,39,r_ge-r_gf);

    
    // Balloons pressures
    gsl_matrix_set(J,37,37,1);
    gsl_matrix_set(J,38,38,1);
    gsl_matrix_set(J,39,39,1);

    return GSL_SUCCESS;
}


int system_3_levels_fdf(const gsl_vector *x, void *p, gsl_vector *f, gsl_matrix *J)
{
    system_3_levels_f(x,p,f);
    system_3_levels_df(x,p,J);

    return GSL_SUCCESS;
}



static int __compute_init_levels_positions(const gsl_vector *A, const gsl_vector *B, double r_top, double r_mid, double r_bot, double phi_ad, double phi_dc, double phi_df, double phi_fe, struct system_3_levels_params *params, gsl_vector *x0)
{
    gsl_vector *center_top = gsl_vector_alloc(2);
    if(center_from_points_and_radius(A,B,r_top,center_top) != GSL_SUCCESS)
    {
        fprintf(stderr,"Failed to compute top level center point\n");
        exit(1);
    }

    gsl_vector *uX = gsl_vector_calloc(2);
    gsl_vector_set(uX,0,-1);
    gsl_vector *cA = vector_centred(A,center_top);
    double a_ad;
    if(vectors_ang_clockwise(uX,cA,&a_ad) != GSL_SUCCESS)
    {
        fprintf(stderr,"Failed to compute alpha AD\n");
        exit(1);
    }

    double a_dc = a_ad+phi_ad;
    double a_cb = a_dc+phi_dc;

    gsl_vector *D = gsl_vector_alloc(2);
    if(point_from_alpha(a_dc,r_top,center_top,D) != GSL_SUCCESS)
    {
        fprintf(stderr,"Failed to compute point D\n");
        exit(1);
    }

    gsl_vector *C = gsl_vector_alloc(2);
    if(point_from_alpha(a_cb,r_top,center_top,C) != GSL_SUCCESS)
    {
        fprintf(stderr,"Failed to compute point C\n");
        exit(1);
    }

    gsl_vector *cB = vector_centred(B,center_top);
    gsl_vector *cC_top = vector_centred(C,center_top);
    double phi_cb;
    if(vectors_ang_clockwise(cC_top,cB,&phi_cb) != GSL_SUCCESS)
    {
        fprintf(stderr,"Failed to compute phi CB\n");
        exit(1);
    }

    gsl_vector *center_mid = gsl_vector_alloc(2);
    if(center_from_points_and_radius(C,D,r_mid,center_mid) != GSL_SUCCESS)
    {
        fprintf(stderr,"Failed to compute middle level center point\n");
        exit(1);
    }

    gsl_vector *cD_mid = vector_centred(D,center_mid);
    double a_df;
    if(vectors_ang_clockwise(uX,cD_mid,&a_df) != GSL_SUCCESS)
    {
        fprintf(stderr,"Failed to compute alpha DF");
        exit(1);
    }

    double a_fe = a_df+phi_df;
    double a_ec = a_fe+phi_fe;

    gsl_vector *F = gsl_vector_alloc(2);
    if(point_from_alpha(a_fe,r_mid,center_mid,F) != GSL_SUCCESS)
    {
        fprintf(stderr,"Failed to compute point F");
        exit(1);
    }

    gsl_vector *E = gsl_vector_alloc(2);
    if(point_from_alpha(a_ec,r_mid,center_mid,E) != GSL_SUCCESS)
    {
        fprintf(stderr,"Failed to compute point E");
        exit(1);
    }

    gsl_vector *cE_mid = vector_centred(E,center_mid);
    gsl_vector *cC_mid = vector_centred(C,center_mid);
    double phi_ec;
    if(vectors_ang_clockwise(cE_mid,cC_mid,&phi_ec) != GSL_SUCCESS)
    {
        fprintf(stderr,"Failed to compute phi EC");
        exit(1);
    }

    gsl_vector *center_bot = gsl_vector_alloc(2);
    if(center_from_points_and_radius(E,F,r_bot,center_bot) != GSL_SUCCESS)
    {
        fprintf(stderr,"Failed to compute bottom level center point\n");
        exit(1);
    }

    gsl_vector *G = gsl_vector_alloc(2);
    gsl_vector_set(G,0,gsl_vector_get(center_bot,0));
    gsl_vector_set(G,1,gsl_vector_get(center_bot,1)-r_bot);
    // gsl_vector_fprintf(stderr,E,"%g");

    gsl_vector *cG = vector_centred(G,center_bot);
    gsl_vector *cF_bot = vector_centred(F,center_bot);
    gsl_vector *cE_bot = vector_centred(E,center_bot);

    double phi_gf;
    if(vectors_ang_clockwise(cF_bot,cG,&phi_gf) != GSL_SUCCESS)
    {
        fprintf(stderr,"Failed to compute phi GF\n");
        exit(1);
    }

    double phi_ge;
    if(vectors_ang_clockwise(cG,cE_bot,&phi_ge) != GSL_SUCCESS)
    {
        fprintf(stderr,"Failed to compute phi GE\n");
        exit(1);
    }

    params->Ax = gsl_vector_get(A,0);
    params->Ay = gsl_vector_get(A,1);
    params->Bx = gsl_vector_get(B,0);
    params->By = gsl_vector_get(B,1);
    params->phi_ad_0 = phi_ad;
    params->phi_cb_0 = phi_cb;
    params->phi_dc_0 = phi_dc;
    params->phi_ec_0 = phi_ec;
    params->phi_fe_0 = phi_fe;
    params->phi_df_0 = phi_df;
    params->phi_ge_0 = phi_ge;
    params->phi_gf_0 = phi_gf;
    params->r_bot_0 = r_bot;
    params->r_mid_0 = r_mid;
    params->r_top_0 = r_top;

    double x_center_top = gsl_vector_get(center_top,0);
    double y_center_top = gsl_vector_get(center_top,1);
    double x_center_mid = gsl_vector_get(center_mid,0);
    double y_center_mid = gsl_vector_get(center_mid,1);
    double x_center_bot = gsl_vector_get(center_bot,0);
    double y_center_bot = gsl_vector_get(center_bot,1);
    gsl_vector_set(x0,0,phi_ad);
    gsl_vector_set(x0,1,r_top);
    gsl_vector_set(x0,2,x_center_top);
    gsl_vector_set(x0,3,y_center_top);
    gsl_vector_set(x0,4,a_ad);
    gsl_vector_set(x0,5,phi_cb);
    gsl_vector_set(x0,6,r_top);
    gsl_vector_set(x0,7,x_center_top);
    gsl_vector_set(x0,8,y_center_top);
    gsl_vector_set(x0,9,a_cb);
    gsl_vector_set(x0,10,phi_dc);
    gsl_vector_set(x0,11,r_top);
    gsl_vector_set(x0,12,x_center_top);
    gsl_vector_set(x0,13,y_center_top);
    gsl_vector_set(x0,14,a_dc);
    gsl_vector_set(x0,15,phi_df);
    gsl_vector_set(x0,16,r_mid);
    gsl_vector_set(x0,17,x_center_mid);
    gsl_vector_set(x0,18,y_center_mid);
    gsl_vector_set(x0,19,a_df);
    gsl_vector_set(x0,20,phi_ec);
    gsl_vector_set(x0,21,r_mid);
    gsl_vector_set(x0,22,x_center_mid);
    gsl_vector_set(x0,23,y_center_mid);
    gsl_vector_set(x0,24,a_ec);
    gsl_vector_set(x0,25,phi_fe);
    gsl_vector_set(x0,26,r_mid);
    gsl_vector_set(x0,27,x_center_mid);
    gsl_vector_set(x0,28,y_center_mid);
    gsl_vector_set(x0,29,a_fe);
    gsl_vector_set(x0,30,phi_ge);
    gsl_vector_set(x0,31,r_bot);
    gsl_vector_set(x0,32,y_center_bot);
    gsl_vector_set(x0,33,phi_gf);
    gsl_vector_set(x0,34,r_bot);
    gsl_vector_set(x0,35,y_center_bot);
    gsl_vector_set(x0,36,x_center_bot);

    gsl_vector_free(cA);
    gsl_vector_free(cB);
    gsl_vector_free(cC_mid);
    gsl_vector_free(cC_top);
    gsl_vector_free(C);
    gsl_vector_free(cD_mid);
    gsl_vector_free(D);
    gsl_vector_free(E);
    gsl_vector_free(cE_mid);
    gsl_vector_free(cE_bot);
    gsl_vector_free(F);
    gsl_vector_free(cF_bot);
    gsl_vector_free(G);
    gsl_vector_free(cG);
    gsl_vector_free(uX);
    gsl_vector_free(center_bot);
    gsl_vector_free(center_mid);
    gsl_vector_free(center_top);

    return GSL_SUCCESS;
}


int system_3_levels_compute_init_config(double Ax, double Ay, double Bx, double By, double phi_ad, double phi_dc, double phi_df, double phi_fe, double r_top, double r_mid, double r_bot, double k, double p_top, double p_mid, double p_bot, double p, double p_ac, struct system_3_levels_params *params, gsl_vector *x0)
{
    gsl_vector *A = gsl_vector_alloc(2);
    gsl_vector_set(A,0,Ax);
    gsl_vector_set(A,1,Ay);
    gsl_vector *B = gsl_vector_alloc(2);
    gsl_vector_set(B,0,Bx);
    gsl_vector_set(B,1,By);
    if(__compute_init_levels_positions(A,B,r_top,r_mid,r_bot,phi_ad,phi_dc,phi_df,phi_fe,params,x0) != GSL_SUCCESS)
    {
        fprintf(stderr,"Failed to compute initial levels position\n");
        exit(1);
    }
    params->p_top_0 = p_top;
    params->p_mid_0 = p_mid;
    params->p_bot_0 = p_bot;
    params->p_ac = p_ac;

    gsl_vector_set(x0,37,p_top);
    gsl_vector_set(x0,38,p_mid);
    gsl_vector_set(x0,39,p_bot);

    gsl_vector_free(A);
    gsl_vector_free(B);

    return GSL_SUCCESS;
}


int system_3_levels_eval()
{
    struct system_3_levels_params params;
    gsl_vector *x0 = gsl_vector_alloc(N_eq);
    system_3_levels_compute_init_config(0.482,1.4,0.28,0.85,3.129,1.162,1.8,1.0,0.5,0.35,0.2,1,20000,6500,3000,101325,1000,&params,x0);

    const gsl_multiroot_fdfsolver_type *T = gsl_multiroot_fdfsolver_newton;
    gsl_multiroot_fdfsolver *s = gsl_multiroot_fdfsolver_alloc(T,N_eq);

    gsl_multiroot_function_fdf fdf;
    fdf.f = system_3_levels_f;
    fdf.df = system_3_levels_df;
    fdf.fdf = system_3_levels_fdf;
    fdf.n = N_eq;
    fdf.params = &params;

    FILE *f_diff = fopen("diff.txt","w");
    print_J_diff(f_diff,x0,&params,fdf.f,fdf.df);
    fclose(f_diff);

    FILE *fx0 = fopen("3_levels_init.txt","w");
    gsl_vector_fprintf(fx0,x0,"%g");
    fclose(fx0);

    gsl_multiroot_fdfsolver_set(s,&fdf,x0);

    size_t max_iters = 100;
    size_t iter = 0;
    double eps = 1e-7;
    int status;
    do
    {
        status = gsl_multiroot_fdfsolver_iterate(s);
        if(status)
            break;
        
        status = gsl_multiroot_test_residual(s->f,eps);
        ++iter;
    } while(status == GSL_CONTINUE && iter < max_iters);

    FILE *fx = fopen("3_levels.txt","w");
    gsl_vector_fprintf(fx,s->x,"%g");
    fclose(fx);

    gsl_vector_free(x0);
    gsl_multiroot_fdfsolver_free(s);

    return GSL_SUCCESS;
}