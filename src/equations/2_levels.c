#include <equations/2_levels.h>
#include <equations/utils.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>
#include <gsl/gsl_sf_trig.h>
#include <gsl/gsl_blas.h>

static const size_t N_eq = 24;


void __system_2_levels_f_general(const gsl_vector *x, const struct system_2_levels_params *params, gsl_vector *f)
{
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
    const double phi_ed = gsl_vector_get(x,15);
    const double r_ed = gsl_vector_get(x,16);
    const double y_ed = gsl_vector_get(x,17);
    const double phi_ec = gsl_vector_get(x,18);
    const double r_ec = gsl_vector_get(x,19);
    const double y_ec = gsl_vector_get(x,20);
    const double x_bot = gsl_vector_get(x,21);
    const double p_top = gsl_vector_get(x,22);
    const double p_bot = gsl_vector_get(x,23);

    const double a_ec = 3*M_PI_2;
    const double a_ed = 3*M_PI_2;

    // Preservation of length
    double eq1 = r_ad*phi_ad - params->r_top_0*params->phi_ad_0;
    gsl_vector_set(f,0,eq1);
    double eq2 = r_cb*phi_cb - params->r_top_0*params->phi_cb_0;
    gsl_vector_set(f,1,eq2);
    double eq3 = r_dc*phi_dc - params->r_top_0*params->phi_dc_0;
    gsl_vector_set(f,2,eq3);
    double eq4 = r_ec*phi_ec + r_ed*phi_ed - params->r_bot_0*params->phi_ec_0 - params->r_bot_0*params->phi_ed_0;
    gsl_vector_set(f,3,eq4);

    // Point A continuity
    double eq5 = -r_ad*gsl_sf_cos(a_ad) + x_ad - params->Ax;
    gsl_vector_set(f,4,eq5);
    double eq6 = r_ad*gsl_sf_sin(a_ad) + y_ad - params->Ay;
    gsl_vector_set(f,5,eq6);

    // Point B continuity
    double eq7 = -r_cb*gsl_sf_cos(a_cb+phi_cb) + x_cb - params->Bx;
    gsl_vector_set(f,6,eq7);
    double eq8 = r_cb*gsl_sf_sin(a_cb+phi_cb) + y_cb - params->By;
    gsl_vector_set(f,7,eq8);

    // Point C continuity
    double eq9 = -r_dc*gsl_sf_cos(a_dc+phi_dc) + x_dc - (-r_cb*gsl_sf_cos(a_cb) + x_cb);
    gsl_vector_set(f,8,eq9);
    double eq10 = r_dc*gsl_sf_sin(a_dc+phi_dc) + y_dc - (r_cb*gsl_sf_sin(a_cb) + y_cb);
    gsl_vector_set(f,9,eq10);
    double eq11 = -r_dc*gsl_sf_cos(a_dc+phi_dc) + x_dc - (-r_ec*gsl_sf_cos(a_ec+phi_ec) + x_bot);
    gsl_vector_set(f,10,eq11);
    double eq12 = r_dc*gsl_sf_sin(a_dc+phi_dc) + y_dc - (r_ec*gsl_sf_sin(a_ec+phi_ec) + y_ec);
    gsl_vector_set(f,11,eq12);

    // Point D continuity
    double eq13 = -r_ad*gsl_sf_cos(a_ad+phi_ad) + x_ad - (-r_dc*gsl_sf_cos(a_dc) + x_dc);
    gsl_vector_set(f,12,eq13);
    double eq14 = r_ad*gsl_sf_sin(a_ad+phi_ad) + y_ad - (r_dc*gsl_sf_sin(a_dc) + y_dc);
    gsl_vector_set(f,13,eq14);
    double eq15 = -r_ad*gsl_sf_cos(a_ad+phi_ad) + x_ad - (-r_ed*gsl_sf_cos(a_ed-phi_ed) + x_bot);
    gsl_vector_set(f,14,eq15);
    double eq16 = r_ad*gsl_sf_sin(a_ad+phi_ad) + y_ad - (r_ed*gsl_sf_sin(a_ed-phi_ed) + y_ed);
    gsl_vector_set(f,15,eq16);

    // Point E continuity
    double eq17 = (-r_ed + y_ed) - (-r_ec + y_ec);
    gsl_vector_set(f,16,eq17);

    // Point D steadiness
    double eq18 = -p_top*r_ad*gsl_sf_sin(a_ad+phi_ad) + (p_top-p_bot)*r_dc*gsl_sf_sin(a_dc) + p_bot*r_ed*gsl_sf_sin(a_ed-phi_ed);
    gsl_vector_set(f,17,eq18);
    double eq19 = -p_top*r_ad*gsl_sf_cos(a_ad+phi_ad) + (p_top-p_bot)*r_dc*gsl_sf_cos(a_dc) + p_bot*r_ed*gsl_sf_cos(a_ed-phi_ed);
    gsl_vector_set(f,18,eq19);

    // Point C steadiness
    double eq20 = (p_top-params->p_ac)*r_cb*gsl_sf_sin(a_cb) - (p_top-p_bot)*r_dc*gsl_sf_sin(a_dc+phi_dc) - (p_bot-params->p_ac)*r_ec*gsl_sf_sin(a_ec+phi_ec);
    gsl_vector_set(f,19,eq20);
    double eq21 = (p_top-params->p_ac)*r_cb*gsl_sf_cos(a_cb) - (p_top-p_bot)*r_dc*gsl_sf_cos(a_dc+phi_dc) - (p_bot-params->p_ac)*r_ec*gsl_sf_cos(a_ec+phi_ec);
    gsl_vector_set(f,20,eq21);

    // Point E steadiness
    double eq22 = (p_bot-params->p_ac)*r_ec - p_bot*r_ed;
    gsl_vector_set(f,21,eq22);
}


void __system_2_levels_df_general(const gsl_vector *x, const struct system_2_levels_params *params, gsl_matrix *J)
{
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
    const double phi_ed = gsl_vector_get(x,15);
    const double r_ed = gsl_vector_get(x,16);
    //const double y_ed = gsl_vector_get(x,17);
    const double phi_ec = gsl_vector_get(x,18);
    const double r_ec = gsl_vector_get(x,19);
    //const double y_ec = gsl_vector_get(x,20);
    //const double x_bot = gsl_vector_get(x,21);
    const double p_top = gsl_vector_get(x,22);
    const double p_bot = gsl_vector_get(x,23);

    const double a_ec = 3*M_PI_2;
    const double a_ed = 3*M_PI_2;

    gsl_matrix_set_all(J,0);

    // Preservation of length
    gsl_matrix_set(J,0,0,r_ad);
    gsl_matrix_set(J,0,1,phi_ad);
    
    gsl_matrix_set(J,1,5,r_cb);
    gsl_matrix_set(J,1,6,phi_cb);

    gsl_matrix_set(J,2,10,r_dc);
    gsl_matrix_set(J,2,11,phi_dc);
    
    gsl_matrix_set(J,3,15,r_ed);
    gsl_matrix_set(J,3,16,phi_ed);
    gsl_matrix_set(J,3,18,r_ec);
    gsl_matrix_set(J,3,19,phi_ec);


    // Point A continuity
    gsl_matrix_set(J,4,1,-gsl_sf_cos(a_ad));
    gsl_matrix_set(J,4,2,1);
    gsl_matrix_set(J,4,4,r_ad*gsl_sf_sin(a_ad));
    
    gsl_matrix_set(J,5,1,gsl_sf_sin(a_ad));
    gsl_matrix_set(J,5,3,1);
    gsl_matrix_set(J,5,4,r_ad*gsl_sf_cos(a_ad));
    

    // Point B continuity
    gsl_matrix_set(J,6,5,r_cb*gsl_sf_sin(phi_cb+a_cb));
    gsl_matrix_set(J,6,6,-gsl_sf_cos(phi_cb+a_cb));
    gsl_matrix_set(J,6,7,1);
    gsl_matrix_set(J,6,9,r_cb*gsl_sf_sin(phi_cb+a_cb));
    
    gsl_matrix_set(J,7,5,r_cb*gsl_sf_cos(phi_cb+a_cb));
    gsl_matrix_set(J,7,6,gsl_sf_sin(phi_cb+a_cb));
    gsl_matrix_set(J,7,8,1);
    gsl_matrix_set(J,7,9,r_cb*gsl_sf_cos(phi_cb+a_cb));
    

    // Point C continuity
    gsl_matrix_set(J,8,6,gsl_sf_cos(a_cb));
    gsl_matrix_set(J,8,7,-1);
    gsl_matrix_set(J,8,9,-r_cb*gsl_sf_sin(a_cb));
    gsl_matrix_set(J,8,10,r_dc*gsl_sf_sin(phi_dc+a_dc));
    gsl_matrix_set(J,8,11,-gsl_sf_cos(phi_dc+a_dc));
    gsl_matrix_set(J,8,12,1);
    gsl_matrix_set(J,8,14,r_dc*gsl_sf_sin(phi_dc+a_dc));

    gsl_matrix_set(J,9,6,-gsl_sf_sin(a_cb));
    gsl_matrix_set(J,9,8,-1);
    gsl_matrix_set(J,9,9,-r_cb*gsl_sf_cos(a_cb));
    gsl_matrix_set(J,9,10,r_dc*gsl_sf_cos(phi_dc+a_dc));
    gsl_matrix_set(J,9,11,gsl_sf_sin(phi_dc+a_dc));
    gsl_matrix_set(J,9,13,1);
    gsl_matrix_set(J,9,14,r_dc*gsl_sf_cos(phi_dc+a_dc));
    
    gsl_matrix_set(J,10,10,r_dc*gsl_sf_sin(phi_dc+a_dc));
    gsl_matrix_set(J,10,11,-gsl_sf_cos(phi_dc+a_dc));
    gsl_matrix_set(J,10,12,1);
    gsl_matrix_set(J,10,14,r_dc*gsl_sf_sin(phi_dc+a_dc));
    gsl_matrix_set(J,10,18,-r_ec*gsl_sf_sin(a_ec+phi_ec));
    gsl_matrix_set(J,10,19,gsl_sf_cos(a_ec+phi_ec));
    gsl_matrix_set(J,10,21,-1);

    gsl_matrix_set(J,11,10,r_dc*gsl_sf_cos(phi_dc+a_dc));
    gsl_matrix_set(J,11,11,gsl_sf_sin(phi_dc+a_dc));
    gsl_matrix_set(J,11,13,1);
    gsl_matrix_set(J,11,14,r_dc*gsl_sf_cos(phi_dc+a_dc));
    gsl_matrix_set(J,11,18,-r_ec*gsl_sf_cos(a_ec+phi_ec));
    gsl_matrix_set(J,11,19,-gsl_sf_sin(a_ec+phi_ec));
    gsl_matrix_set(J,11,20,-1);


    // Point D continuity
    gsl_matrix_set(J,12,0,r_ad*gsl_sf_sin(phi_ad+a_ad));
    gsl_matrix_set(J,12,1,-gsl_sf_cos(phi_ad+a_ad));
    gsl_matrix_set(J,12,2,1);
    gsl_matrix_set(J,12,4,r_ad*gsl_sf_sin(phi_ad+a_ad));
    gsl_matrix_set(J,12,11,gsl_sf_cos(a_dc));
    gsl_matrix_set(J,12,12,-1);
    gsl_matrix_set(J,12,14,-r_dc*gsl_sf_sin(a_dc));

    gsl_matrix_set(J,13,0,r_ad*gsl_sf_cos(phi_ad+a_ad));
    gsl_matrix_set(J,13,1,gsl_sf_sin(phi_ad+a_ad));
    gsl_matrix_set(J,13,3,1);
    gsl_matrix_set(J,13,4,r_ad*gsl_sf_cos(phi_ad+a_ad));
    gsl_matrix_set(J,13,11,-gsl_sf_sin(a_dc));
    gsl_matrix_set(J,13,13,-1);
    gsl_matrix_set(J,13,14,-r_dc*gsl_sf_cos(a_dc));

    gsl_matrix_set(J,14,0,r_ad*gsl_sf_sin(phi_ad+a_ad));
    gsl_matrix_set(J,14,1,-gsl_sf_cos(phi_ad+a_ad));
    gsl_matrix_set(J,14,2,1);
    gsl_matrix_set(J,14,4,r_ad*gsl_sf_sin(phi_ad+a_ad));
    gsl_matrix_set(J,14,15,r_ed*gsl_sf_sin(a_ed-phi_ed));
    gsl_matrix_set(J,14,16,gsl_sf_cos(a_ed-phi_ed));
    gsl_matrix_set(J,14,21,-1);

    gsl_matrix_set(J,15,0,r_ad*gsl_sf_cos(phi_ad+a_ad));
    gsl_matrix_set(J,15,1,gsl_sf_sin(phi_ad+a_ad));
    gsl_matrix_set(J,15,3,1);
    gsl_matrix_set(J,15,4,r_ad*gsl_sf_cos(phi_ad+a_ad));
    gsl_matrix_set(J,15,15,r_ed*gsl_sf_cos(a_ed-phi_ed));
    gsl_matrix_set(J,15,16,-gsl_sf_sin(a_ed-phi_ed));
    gsl_matrix_set(J,15,17,-1);

    // Point E continuity
    gsl_matrix_set(J,16,16,-1);
    gsl_matrix_set(J,16,17,1);
    gsl_matrix_set(J,16,19,1);
    gsl_matrix_set(J,16,20,-1);


    // Point D steadiness
    gsl_matrix_set(J,17,0,-p_top*r_ad*gsl_sf_cos(a_ad+phi_ad));
    gsl_matrix_set(J,17,1,-p_top*gsl_sf_sin(a_ad+phi_ad));
    gsl_matrix_set(J,17,4,-p_top*r_ad*gsl_sf_cos(a_ad+phi_ad));
    gsl_matrix_set(J,17,11,(p_top-p_bot)*gsl_sf_sin(a_dc));
    gsl_matrix_set(J,17,14,(p_top-p_bot)*r_dc*gsl_sf_cos(a_dc));
    gsl_matrix_set(J,17,15,-p_bot*r_ed*gsl_sf_cos(a_ed-phi_ed));
    gsl_matrix_set(J,17,16,p_bot*gsl_sf_sin(a_ed-phi_ed));
    gsl_matrix_set(J,17,22,-r_ad*gsl_sf_sin(a_ad+phi_ad) + r_dc*gsl_sf_sin(a_dc));
    gsl_matrix_set(J,17,23,-r_dc*gsl_sf_sin(a_dc) + r_ed*gsl_sf_sin(a_ed-phi_ed));

    gsl_matrix_set(J,18,0,p_top*r_ad*gsl_sf_sin(a_ad+phi_ad));
    gsl_matrix_set(J,18,1,-p_top*gsl_sf_cos(a_ad+phi_ad));
    gsl_matrix_set(J,18,4,p_top*r_ad*gsl_sf_sin(a_ad+phi_ad));
    gsl_matrix_set(J,18,11,(p_top-p_bot)*gsl_sf_cos(a_dc));
    gsl_matrix_set(J,18,14,-(p_top-p_bot)*r_dc*gsl_sf_sin(a_dc));
    gsl_matrix_set(J,18,15,p_bot*r_ed*gsl_sf_sin(a_ed-phi_ed));
    gsl_matrix_set(J,18,16,p_bot*gsl_sf_cos(a_ed-phi_ed));
    gsl_matrix_set(J,18,22,-r_ad*gsl_sf_cos(a_ad+phi_ad) + r_dc*gsl_sf_cos(a_dc));
    gsl_matrix_set(J,18,23,-r_dc*gsl_sf_cos(a_dc) + r_ed*gsl_sf_cos(a_ed-phi_ed));


    // Point C steadiness
    gsl_matrix_set(J,19,6,(p_top-params->p_ac)*gsl_sf_sin(a_cb));
    gsl_matrix_set(J,19,9,(p_top-params->p_ac)*r_cb*gsl_sf_cos(a_cb));
    gsl_matrix_set(J,19,10,-(p_top-p_bot)*r_dc*gsl_sf_cos(a_dc+phi_dc));
    gsl_matrix_set(J,19,11,-(p_top-p_bot)*gsl_sf_sin(a_dc+phi_dc));
    gsl_matrix_set(J,19,14,-(p_top-p_bot)*r_dc*gsl_sf_cos(a_dc+phi_dc));
    gsl_matrix_set(J,19,18,-(p_bot-params->p_ac)*r_ec*gsl_sf_cos(a_ec+phi_ec));
    gsl_matrix_set(J,19,19,-(p_bot-params->p_ac)*gsl_sf_sin(a_ec+phi_ec));
    gsl_matrix_set(J,19,22,r_cb*gsl_sf_sin(a_cb) - r_dc*gsl_sf_sin(a_dc+phi_dc));
    gsl_matrix_set(J,19,23,r_dc*gsl_sf_sin(a_dc+phi_dc) - r_ec*gsl_sf_sin(a_ec+phi_ec));

    gsl_matrix_set(J,20,6,(p_top-params->p_ac)*gsl_sf_cos(a_cb));
    gsl_matrix_set(J,20,9,-(p_top-params->p_ac)*r_cb*gsl_sf_sin(a_cb));
    gsl_matrix_set(J,20,10,(p_top-p_bot)*r_dc*gsl_sf_sin(a_dc+phi_dc));
    gsl_matrix_set(J,20,11,-(p_top-p_bot)*gsl_sf_cos(a_dc+phi_dc));
    gsl_matrix_set(J,20,14,(p_top-p_bot)*r_dc*gsl_sf_sin(a_dc+phi_dc));
    gsl_matrix_set(J,20,18,(p_bot-params->p_ac)*r_ec*gsl_sf_sin(a_ec+phi_ec));
    gsl_matrix_set(J,20,19,-(p_bot-params->p_ac)*gsl_sf_cos(a_ec+phi_ec));
    gsl_matrix_set(J,20,22,r_cb*gsl_sf_cos(a_cb) - r_dc*gsl_sf_cos(a_dc+phi_dc));
    gsl_matrix_set(J,20,23,r_dc*gsl_sf_cos(a_dc+phi_dc) - r_ec*gsl_sf_cos(a_ec+phi_ec));


    // Point E steadiness
    gsl_matrix_set(J,21,16,-p_bot);
    gsl_matrix_set(J,21,19,(p_bot-params->p_ac));
    gsl_matrix_set(J,21,23,r_ec-r_ed);
}


int system_2_levels_f(const gsl_vector *x, void *p, gsl_vector *f)
{
    struct system_2_levels_params *params = (struct system_2_levels_params*)p;
    
    const double p_top = gsl_vector_get(x,22);
    const double p_bot = gsl_vector_get(x,23);

    __system_2_levels_f_general(x,params,f);

    
    // Balloons pressures
   
    double eq23 = p_top - params->p_top_0;
    gsl_vector_set(f,22,eq23);
    double eq24 = p_bot - params->p_bot_0;
    gsl_vector_set(f,23,eq24);

    return GSL_SUCCESS;
}


int system_2_levels_df(const gsl_vector *x, void *p, gsl_matrix *J)
{
    //const size_t N_eq = 24;
    struct system_2_levels_params *params = (struct system_2_levels_params*)p;

    __system_2_levels_df_general(x,params,J);
    
    // Balloons pressures
    gsl_matrix_set(J,22,22,1);
    gsl_matrix_set(J,23,23,1);

    return GSL_SUCCESS;
}


int system_2_levels_fdf(const gsl_vector *x, void *p, gsl_vector *f, gsl_matrix *J)
{
    system_2_levels_f(x,p,f);
    system_2_levels_df(x,p,J);

    return GSL_SUCCESS;
}


const int N_top_PpA = 5;    // Points per Arc
const int N_top = N_top_PpA*3;
const int N_bot_PpA = 5;    // Points per Arc
const int N_bot = N_bot_PpA*2;


int system_2_levels_adiabatic_f(const gsl_vector *x, void *p, gsl_vector *f)
{
    struct system_2_levels_params *params = (struct system_2_levels_params*)p;
    
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
    const double phi_ed = gsl_vector_get(x,15);
    const double r_ed = gsl_vector_get(x,16);
    const double y_ed = gsl_vector_get(x,17);
    const double phi_ec = gsl_vector_get(x,18);
    const double r_ec = gsl_vector_get(x,19);
    const double y_ec = gsl_vector_get(x,20);
    const double x_bot = gsl_vector_get(x,21);
    const double p_top = gsl_vector_get(x,22);
    const double p_bot = gsl_vector_get(x,23);

    const double a_ec = 3*M_PI_2;
    const double a_ed = 3*M_PI_2;

    __system_2_levels_f_general(x,params,f);

    double xp_top[N_top+1], yp_top[N_top+1];
    double xp_bot[N_bot+1], yp_bot[N_bot+1];
    for(int pt_i = 0; pt_i < N_top_PpA; ++pt_i)
    {
        double ang_step = (double)pt_i/(double)N_top_PpA;
        double ang_ad = a_ad + ang_step*phi_ad;
        xp_top[pt_i] = x_ad - r_ad*gsl_sf_cos(ang_ad);
        yp_top[pt_i] = y_ad + r_ad*gsl_sf_sin(ang_ad);

        double ang_dc = a_dc + ang_step*phi_dc;
        xp_top[pt_i+N_top_PpA] = x_dc - r_dc*gsl_sf_cos(ang_dc);
        yp_top[pt_i+N_top_PpA] = y_dc + r_dc*gsl_sf_sin(ang_dc);
    
        double ang_cb = a_cb + ang_step*phi_cb;
        xp_top[pt_i+N_top_PpA*2] = x_cb - r_cb*gsl_sf_cos(ang_cb);
        yp_top[pt_i+N_top_PpA*2] = y_cb + r_cb*gsl_sf_sin(ang_cb);
    }
    xp_top[N_top] = xp_top[0];
    yp_top[N_top] = yp_top[0];

    for(int pb_i = 0; pb_i < N_bot_PpA; ++pb_i)
    {
        double ang_step = (double)pb_i/(double)N_bot_PpA;
        double ang_ec = a_ec + ang_step*phi_ec;
        xp_bot[pb_i] = x_bot - r_ec*gsl_sf_cos(ang_ec);
        yp_bot[pb_i] = y_ec + r_ec*gsl_sf_sin(ang_ec);

        double ang_ed = a_ed + (ang_step-1)*phi_ed;
        xp_bot[pb_i+N_bot_PpA] = x_bot - r_ed*gsl_sf_cos(ang_ed);
        yp_bot[pb_i+N_bot_PpA] = y_ed + r_ed*gsl_sf_sin(ang_ed);
    }
    xp_bot[N_bot] = xp_bot[0];
    yp_bot[N_bot] = yp_bot[0];

    double S_top = 0.0;
    for(int pt_i = 0; pt_i < N_top; ++pt_i)
    {
        double v1_x = xp_top[pt_i];
        double v1_y = yp_top[pt_i];
        double v2_x = xp_top[pt_i+1];
        double v2_y = yp_top[pt_i+1];
        S_top += v1_x*v2_y - v1_y*v2_x;
    }
    S_top /= 2;
    S_top = fabs(S_top);

    double S_bot = 0.0;
    for(int pt_i = 0; pt_i < N_bot; ++pt_i)
    {
        double v1_x = xp_bot[pt_i];
        double v1_y = yp_bot[pt_i];
        double v2_x = xp_bot[pt_i+1];
        double v2_y = yp_bot[pt_i+1];
        S_bot += v1_x*v2_y - v1_y*v2_x;
    }
    S_bot /= 2;
    S_bot = fabs(S_bot);

    //printf("S_top: %f; S_bot: %f\n",S_top,S_bot);
    //printf("S_top_0: %f; S_bot_0: %f\n\n",params->S_top_0,params->S_bot_0);

    double eq23 = (params->p_top_0)*pow((params->S_top_0),params->k) - (p_top)*pow(S_top,params->k);
    gsl_vector_set(f,22,eq23);
    double eq24 = (params->p_bot_0)*pow((params->S_bot_0),params->k) - (p_bot)*pow(S_bot,params->k);
    gsl_vector_set(f,23,eq24);

    return GSL_SUCCESS;
}


int system_2_levels_adiabatic_df(const gsl_vector *x, void *p, gsl_matrix *J)
{
    struct system_2_levels_params *params = (struct system_2_levels_params*)p;
    
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
    const double phi_ed = gsl_vector_get(x,15);
    const double r_ed = gsl_vector_get(x,16);
    const double y_ed = gsl_vector_get(x,17);
    const double phi_ec = gsl_vector_get(x,18);
    const double r_ec = gsl_vector_get(x,19);
    const double y_ec = gsl_vector_get(x,20);
    const double x_bot = gsl_vector_get(x,21);
    const double p_top = gsl_vector_get(x,22);
    const double p_bot = gsl_vector_get(x,23);

    const double a_ec = 3*M_PI_2;
    const double a_ed = 3*M_PI_2;

    __system_2_levels_df_general(x,params,J);

    double xp_top[N_top+1], yp_top[N_top+1];
    double xp_bot[N_bot+1], yp_bot[N_bot+1];
    gsl_vector *grad_xp_top[N_top+1], *grad_yp_top[N_top+1];
    gsl_vector *grad_xp_bot[N_bot+1], *grad_yp_bot[N_bot+1];
    for(int i = 0; i < N_top+1; ++i)
    {
        grad_xp_top[i] = gsl_vector_calloc(N_eq);
        grad_yp_top[i] = gsl_vector_calloc(N_eq);
    }
    for(int i = 0; i < N_bot+1; ++i)
    {
        grad_xp_bot[i] = gsl_vector_calloc(N_eq);
        grad_yp_bot[i] = gsl_vector_calloc(N_eq);
    }
    for(int pt_i = 0; pt_i < N_top_PpA; ++pt_i)
    {
        double ang_step = (double)pt_i/(double)N_top_PpA;
        double ang_ad = a_ad + ang_step*phi_ad;
        xp_top[pt_i] = x_ad - r_ad*gsl_sf_cos(ang_ad);
        yp_top[pt_i] = y_ad + r_ad*gsl_sf_sin(ang_ad);
        gsl_vector_set(grad_xp_top[pt_i],0,r_ad*ang_step*gsl_sf_sin(ang_ad));
        gsl_vector_set(grad_xp_top[pt_i],1,-gsl_sf_cos(ang_ad));
        gsl_vector_set(grad_xp_top[pt_i],2,1);
        gsl_vector_set(grad_xp_top[pt_i],4,r_ad*gsl_sf_sin(ang_ad));
        gsl_vector_set(grad_yp_top[pt_i],0,r_ad*ang_step*gsl_sf_cos(ang_ad));
        gsl_vector_set(grad_yp_top[pt_i],1,gsl_sf_sin(ang_ad));
        gsl_vector_set(grad_yp_top[pt_i],3,1);
        gsl_vector_set(grad_yp_top[pt_i],4,r_ad*gsl_sf_cos(ang_ad));

        double ang_dc = a_dc + ang_step*phi_dc;
        xp_top[pt_i+N_top_PpA] = x_dc - r_dc*gsl_sf_cos(ang_dc);
        yp_top[pt_i+N_top_PpA] = y_dc + r_dc*gsl_sf_sin(ang_dc);
        gsl_vector_set(grad_xp_top[pt_i+N_top_PpA],10,r_dc*ang_step*gsl_sf_sin(ang_dc));
        gsl_vector_set(grad_xp_top[pt_i+N_top_PpA],11,-gsl_sf_cos(ang_dc));
        gsl_vector_set(grad_xp_top[pt_i+N_top_PpA],12,1);
        gsl_vector_set(grad_xp_top[pt_i+N_top_PpA],14,r_dc*gsl_sf_sin(ang_dc));
        gsl_vector_set(grad_yp_top[pt_i+N_top_PpA],10,r_dc*ang_step*gsl_sf_cos(ang_dc));
        gsl_vector_set(grad_yp_top[pt_i+N_top_PpA],11,gsl_sf_sin(ang_dc));
        gsl_vector_set(grad_yp_top[pt_i+N_top_PpA],13,1);
        gsl_vector_set(grad_yp_top[pt_i+N_top_PpA],14,r_dc*gsl_sf_cos(ang_dc));
    
        double ang_cb = a_cb + ang_step*phi_cb;
        xp_top[pt_i+N_top_PpA*2] = x_cb - r_cb*gsl_sf_cos(ang_cb);
        yp_top[pt_i+N_top_PpA*2] = y_cb + r_cb*gsl_sf_sin(ang_cb);
        gsl_vector_set(grad_xp_top[pt_i+N_top_PpA*2],5,r_cb*ang_step*gsl_sf_sin(ang_cb));
        gsl_vector_set(grad_xp_top[pt_i+N_top_PpA*2],6,-gsl_sf_cos(ang_cb));
        gsl_vector_set(grad_xp_top[pt_i+N_top_PpA*2],7,1);
        gsl_vector_set(grad_xp_top[pt_i+N_top_PpA*2],9,r_cb*gsl_sf_sin(ang_cb));
        gsl_vector_set(grad_yp_top[pt_i+N_top_PpA*2],5,r_cb*ang_step*gsl_sf_cos(ang_cb));
        gsl_vector_set(grad_yp_top[pt_i+N_top_PpA*2],6,gsl_sf_sin(ang_cb));
        gsl_vector_set(grad_yp_top[pt_i+N_top_PpA*2],8,1);
        gsl_vector_set(grad_yp_top[pt_i+N_top_PpA*2],9,r_cb*gsl_sf_cos(ang_cb));
    }
    xp_top[N_top] = xp_top[0];
    yp_top[N_top] = yp_top[0];
    gsl_vector_memcpy(grad_xp_top[N_top],grad_xp_top[0]);
    gsl_vector_memcpy(grad_yp_top[N_top],grad_yp_top[0]);

    for(int pb_i = 0; pb_i < N_bot_PpA; ++pb_i)
    {
        double ang_step = (double)pb_i/(double)N_bot_PpA;
        double ang_ec = a_ec + ang_step*phi_ec;
        xp_bot[pb_i] = x_bot - r_ec*gsl_sf_cos(ang_ec);
        yp_bot[pb_i] = y_ec + r_ec*gsl_sf_sin(ang_ec);
        gsl_vector_set(grad_xp_bot[pb_i],18,r_ec*ang_step*gsl_sf_sin(ang_ec));
        gsl_vector_set(grad_xp_bot[pb_i],19,-gsl_sf_cos(ang_ec));
        gsl_vector_set(grad_xp_bot[pb_i],21,1);
        gsl_vector_set(grad_yp_bot[pb_i],18,r_ec*ang_step*gsl_sf_cos(ang_ec));
        gsl_vector_set(grad_yp_bot[pb_i],19,gsl_sf_sin(ang_ec));
        gsl_vector_set(grad_yp_bot[pb_i],20,1);

        double ang_ed = a_ed + (ang_step-1)*phi_ed;
        xp_bot[pb_i+N_bot_PpA] = x_bot - r_ed*gsl_sf_cos(ang_ed);
        yp_bot[pb_i+N_bot_PpA] = y_ed + r_ed*gsl_sf_sin(ang_ed);
        gsl_vector_set(grad_xp_bot[pb_i+N_bot_PpA],15,r_ed*(ang_step-1)*gsl_sf_sin(ang_ed));
        gsl_vector_set(grad_xp_bot[pb_i+N_bot_PpA],16,-gsl_sf_cos(ang_ed));
        gsl_vector_set(grad_xp_bot[pb_i+N_bot_PpA],21,1);
        gsl_vector_set(grad_yp_bot[pb_i+N_bot_PpA],15,r_ed*(ang_step-1)*gsl_sf_cos(ang_ed));
        gsl_vector_set(grad_yp_bot[pb_i+N_bot_PpA],16,gsl_sf_sin(ang_ed));
        gsl_vector_set(grad_yp_bot[pb_i+N_bot_PpA],17,1);
    }
    xp_bot[N_bot] = xp_bot[0];
    yp_bot[N_bot] = yp_bot[0];
    gsl_vector_memcpy(grad_xp_bot[N_bot],grad_xp_bot[0]);
    gsl_vector_memcpy(grad_yp_bot[N_bot],grad_yp_bot[0]);

    gsl_vector J_V_top = gsl_matrix_row(J,22).vector;
    gsl_vector *grad_v = gsl_vector_calloc(N_eq);
    double S_top = 0.0;
    for(int pt_i = 0; pt_i < N_top; ++pt_i)
    {
        double v1_x = xp_top[pt_i];
        double v1_y = yp_top[pt_i];
        double v2_x = xp_top[pt_i+1];
        double v2_y = yp_top[pt_i+1];
        S_top += v1_x*v2_y - v1_y*v2_x;
        gsl_vector_memcpy(grad_v,grad_xp_top[pt_i]);
        gsl_vector_scale(grad_v,v2_y);
        gsl_vector_add(&J_V_top,grad_v);
        gsl_vector_memcpy(grad_v,grad_yp_top[pt_i]);
        gsl_vector_scale(grad_v,-v2_x);
        gsl_vector_add(&J_V_top,grad_v);
        gsl_vector_memcpy(grad_v,grad_xp_top[pt_i+1]);
        gsl_vector_scale(grad_v,-v1_y);
        gsl_vector_add(&J_V_top,grad_v);
        gsl_vector_memcpy(grad_v,grad_yp_top[pt_i+1]);
        gsl_vector_scale(grad_v,v1_x);
        gsl_vector_add(&J_V_top,grad_v);
    }
    if(S_top >= 0)
        gsl_vector_scale(&J_V_top,0.5);
    else
        gsl_vector_scale(&J_V_top,-0.5);
    S_top /= 2;
    S_top = fabs(S_top);
    double S_top_pow_k_1 = pow(S_top,params->k-1);
    gsl_vector_scale(&J_V_top,-(p_top)*params->k*S_top_pow_k_1);
    gsl_vector_set(&J_V_top,22,S_top_pow_k_1*S_top);

    gsl_vector J_V_bot = gsl_matrix_row(J,23).vector;
    double S_bot = 0.0;
    for(int pb_i = 0; pb_i < N_bot; ++pb_i)
    {
        double v1_x = xp_bot[pb_i];
        double v1_y = yp_bot[pb_i];
        double v2_x = xp_bot[pb_i+1];
        double v2_y = yp_bot[pb_i+1];
        S_bot += v1_x*v2_y - v1_y*v2_x;
        gsl_vector_memcpy(grad_v,grad_xp_bot[pb_i]);
        gsl_vector_scale(grad_v,v2_y);
        gsl_vector_add(&J_V_bot,grad_v);
        gsl_vector_memcpy(grad_v,grad_yp_bot[pb_i]);
        gsl_vector_scale(grad_v,-v2_x);
        gsl_vector_add(&J_V_bot,grad_v);
        gsl_vector_memcpy(grad_v,grad_xp_bot[pb_i+1]);
        gsl_vector_scale(grad_v,-v1_y);
        gsl_vector_add(&J_V_bot,grad_v);
        gsl_vector_memcpy(grad_v,grad_yp_bot[pb_i+1]);
        gsl_vector_scale(grad_v,v1_x);
        gsl_vector_add(&J_V_bot,grad_v);
    }
    if(S_bot >= 0)
        gsl_vector_scale(&J_V_bot,0.5);
    else
        gsl_vector_scale(&J_V_bot,-0.5);
    S_bot /= 2;
    S_bot = fabs(S_bot);
    double S_bot_pow_k_1 = pow(S_bot,params->k-1);
    gsl_vector_scale(&J_V_bot,-(p_bot)*params->k*S_bot_pow_k_1);
    gsl_vector_set(&J_V_bot,23,S_bot_pow_k_1*S_bot);
    
    gsl_vector_free(grad_v);

    return GSL_SUCCESS;
}


int system_2_levels_adiabatic_fdf(const gsl_vector *x, void *p, gsl_vector *f, gsl_matrix *J)
{
    system_2_levels_adiabatic_f(x,p,f);
    system_2_levels_adiabatic_df(x,p,J);

    return GSL_SUCCESS;
}


int system_2_levels_compute_init_config(const struct system_2_levels_user_params *user_params, gsl_vector *x0, struct system_2_levels_params *params)
{
    params->Ax = user_params->Ax;
    params->Ay = user_params->Ay;
    params->Bx = user_params->Bx;
    params->By = user_params->By;
    params->p_ac = user_params->p_ac;
    params->p_atm = user_params->p_atm;
    params->p_bot_0 = user_params->p_bot_0;
    params->p_top_0 = user_params->p_top_0;
    params->phi_ad_0 = user_params->phi_ad_0;
    params->phi_dc_0 = user_params->phi_dc_0;
    params->r_bot_0 = user_params->r_bot_0;
    params->r_top_0 = user_params->r_top_0;
    params->k = user_params->k;

    gsl_vector *A = gsl_vector_alloc(2);
    gsl_vector_set(A,0,params->Ax);
    gsl_vector_set(A,1,params->Ay);
    gsl_vector *B = gsl_vector_alloc(2);
    gsl_vector_set(B,0,params->Bx);
    gsl_vector_set(B,1,params->By);

    gsl_vector *center_top = gsl_vector_alloc(2);
    if(center_from_points_and_radius(A,B,params->r_top_0,center_top) != GSL_SUCCESS)
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

    double a_dc = a_ad+params->phi_ad_0;
    double a_cb = a_dc+params->phi_dc_0;

    gsl_vector *D = gsl_vector_alloc(2);
    if(point_from_alpha(a_dc,params->r_top_0,center_top,D) != GSL_SUCCESS)
    {
        fprintf(stderr,"Failed to compute point D\n");
        exit(1);
    }

    gsl_vector *C = gsl_vector_alloc(2);
    if(point_from_alpha(a_cb,params->r_top_0,center_top,C) != GSL_SUCCESS)
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

    gsl_vector *center_bot = gsl_vector_alloc(2);
    if(center_from_points_and_radius(C,D,params->r_bot_0,center_bot) != GSL_SUCCESS)
    {
        fprintf(stderr,"Failed to compute bottom level center point\n");
        exit(1);
    }

    gsl_vector *E = gsl_vector_alloc(2);
    gsl_vector_set(E,0,gsl_vector_get(center_bot,0));
    gsl_vector_set(E,1,gsl_vector_get(center_bot,1)-params->r_bot_0);
    // gsl_vector_fprintf(stderr,E,"%g");

    gsl_vector *cE = vector_centred(E,center_bot);
    gsl_vector *cD_bot = vector_centred(D,center_bot);
    gsl_vector *cC_bot = vector_centred(C,center_bot);

    double phi_ed;
    if(vectors_ang_clockwise(cD_bot,cE,&phi_ed) != GSL_SUCCESS)
    {
        fprintf(stderr,"Failed to compute phi ED\n");
        exit(1);
    }

    double phi_ec;
    if(vectors_ang_clockwise(cE,cC_bot,&phi_ec) != GSL_SUCCESS)
    {
        fprintf(stderr,"Failed to compute phi EC\n");
        exit(1);
    }

    gsl_vector *cD_top = vector_centred(D,center_top);
    double S_top = 0, S_bot = 0;
    S_bot += area_between_vectors_triangle(cC_bot,cD_bot);
    double phi_dc_bot = add_angs(phi_ec,phi_ed);
    S_bot += area_segment(phi_dc_bot,params->r_bot_0);
    S_top += area_between_vectors_triangle(cC_top,cD_top);
    S_top += area_between_vectors_triangle(cA,cB);
    S_top += area_segment(phi_cb,params->r_top_0);
    S_top += area_segment(params->phi_ad_0,params->r_top_0);

    params->phi_cb_0 = phi_cb;
    params->phi_ec_0 = phi_ec;
    params->phi_ed_0 = phi_ed;
    params->S_top_0 = S_top;
    params->S_bot_0 = S_bot;

    double x_center_top = gsl_vector_get(center_top,0);
    double y_center_top = gsl_vector_get(center_top,1);
    double x_center_bot = gsl_vector_get(center_bot,0);
    double y_center_bot = gsl_vector_get(center_bot,1);
    gsl_vector_set(x0,0,params->phi_ad_0);
    gsl_vector_set(x0,1,params->r_top_0);
    gsl_vector_set(x0,2,x_center_top);
    gsl_vector_set(x0,3,y_center_top);
    gsl_vector_set(x0,4,a_ad);
    gsl_vector_set(x0,5,phi_cb);
    gsl_vector_set(x0,6,params->r_top_0);
    gsl_vector_set(x0,7,x_center_top);
    gsl_vector_set(x0,8,y_center_top);
    gsl_vector_set(x0,9,a_cb);
    gsl_vector_set(x0,10,params->phi_dc_0);
    gsl_vector_set(x0,11,params->r_top_0);
    gsl_vector_set(x0,12,x_center_top);
    gsl_vector_set(x0,13,y_center_top);
    gsl_vector_set(x0,14,a_dc);
    gsl_vector_set(x0,15,phi_ed);
    gsl_vector_set(x0,16,params->r_bot_0);
    gsl_vector_set(x0,17,y_center_bot);
    gsl_vector_set(x0,18,phi_ec);
    gsl_vector_set(x0,19,params->r_bot_0);
    gsl_vector_set(x0,20,y_center_bot);
    gsl_vector_set(x0,21,x_center_bot);
    gsl_vector_set(x0,22,params->p_top_0);
    gsl_vector_set(x0,23,params->p_bot_0);

    gsl_vector_free(cA);
    gsl_vector_free(cB);
    gsl_vector_free(cC_bot);
    gsl_vector_free(cC_top);
    gsl_vector_free(C);
    gsl_vector_free(cD_bot);
    gsl_vector_free(D);
    gsl_vector_free(cE);
    gsl_vector_free(E);
    gsl_vector_free(uX);
    gsl_vector_free(center_bot);
    gsl_vector_free(center_top);
    gsl_vector_free(A);
    gsl_vector_free(B);

    return GSL_SUCCESS;
}


int system_2_levels_eval_f()
{
    struct system_2_levels_user_params user_params;
    struct system_2_levels_params params;
    gsl_vector *x0 = gsl_vector_alloc(N_eq);
    user_params.Ax = 0.482;
    user_params.Ay = 1.4;
    user_params.Bx = 0.28;
    user_params.By = 0.85;
    user_params.phi_ad_0 = 3.129;
    user_params.phi_dc_0 = 1.162;
    user_params.r_top_0 = 0.5;
    user_params.r_bot_0 = 0.35;
    user_params.p_top_0 = 20000;
    user_params.p_bot_0 = 6500;
    user_params.p_atm = 101325;
    user_params.p_ac = 1500;
    system_2_levels_compute_init_config(&user_params,x0,&params);

    const gsl_multiroot_fdfsolver_type *T = gsl_multiroot_fdfsolver_hybridsj;
    gsl_multiroot_fdfsolver *s = gsl_multiroot_fdfsolver_alloc(T,N_eq);

    gsl_multiroot_function_fdf fdf;
    fdf.f = system_2_levels_f;
    fdf.df = system_2_levels_df;
    fdf.fdf = system_2_levels_fdf;
    fdf.n = N_eq;
    fdf.params = &params;

    FILE *f_diff = fopen("diff.txt","w");
    print_J_diff(f_diff,x0,&params,fdf.f,fdf.df);
    fclose(f_diff);

    FILE *fx0 = fopen("2_levels_init.txt","w");
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

    FILE *fx = fopen("2_levels.txt","w");
    gsl_vector_fprintf(fx,s->x,"%g");
    fclose(fx);

    gsl_multiroot_fdfsolver_free(s);
    gsl_vector_free(x0);

    return GSL_SUCCESS;
}


void system_2_levels_x_to_res(const gsl_vector *x, struct system_2_levels_result *result)
{
    result->phi_ad = gsl_vector_get(x,0);
    result->r_ad = gsl_vector_get(x,1);
    result->x_ad = gsl_vector_get(x,2);
    result->y_ad = gsl_vector_get(x,3);
    result->a_ad = gsl_vector_get(x,4);
    result->phi_cb = gsl_vector_get(x,5);
    result->r_cb = gsl_vector_get(x,6);
    result->x_cb = gsl_vector_get(x,7);
    result->y_cb = gsl_vector_get(x,8);
    result->a_cb = gsl_vector_get(x,9);
    result->phi_dc = gsl_vector_get(x,10);
    result->r_dc = gsl_vector_get(x,11);
    result->x_dc = gsl_vector_get(x,12);
    result->y_dc = gsl_vector_get(x,13);
    result->a_dc = gsl_vector_get(x,14);
    result->phi_ed = gsl_vector_get(x,15);
    result->r_ed = gsl_vector_get(x,16);
    result->y_ed = gsl_vector_get(x,17);
    result->phi_ec = gsl_vector_get(x,18);
    result->r_ec = gsl_vector_get(x,19);
    result->y_ec = gsl_vector_get(x,20);
    result->x_bot = gsl_vector_get(x,21);
    result->p_top = gsl_vector_get(x,22);
    result->p_bot = gsl_vector_get(x,23);
}


int __system_2_levels_eval_general(const struct system_2_levels_user_params *user_params, struct system_2_levels_result *result, gsl_solver_f_t f_ptr, gsl_solver_df_t df_ptr, gsl_solver_fdf_t fdf_ptr)
{
    struct system_2_levels_params params;
    gsl_vector *x0 = gsl_vector_alloc(N_eq);
    system_2_levels_compute_init_config(user_params,x0,&params);

    const gsl_multiroot_fdfsolver_type *T = gsl_multiroot_fdfsolver_hybridsj;
    gsl_multiroot_fdfsolver *s = gsl_multiroot_fdfsolver_alloc(T,N_eq);

    gsl_multiroot_function_fdf fdf;
    fdf.f = f_ptr;
    fdf.df = df_ptr;
    fdf.fdf = fdf_ptr;
    fdf.n = N_eq;
    fdf.params = &params;

    gsl_multiroot_fdfsolver_set(s,&fdf,x0);

    size_t max_iters = 1000;
    size_t iter = 0;
    double eps = 1e-7;
    int status;
    do
    {
        gsl_multiroot_fdfsolver_iterate(s);
        status = gsl_multiroot_test_residual(s->f,eps);
        ++iter;
    } while(status == GSL_CONTINUE && iter < max_iters);

    system_2_levels_x_to_res(s->x,result);

    gsl_multiroot_fdfsolver_free(s);
    gsl_vector_free(x0);

    return GSL_SUCCESS;
}


int system_2_levels_eval(const struct system_2_levels_user_params *user_params, struct system_2_levels_result *result)
{
    if(!user_params || !result)
        return -1;

    return __system_2_levels_eval_general(user_params,result,system_2_levels_f,system_2_levels_df,system_2_levels_fdf);
}


int system_2_levels_adiabatic_eval(const struct system_2_levels_user_params *user_params, struct system_2_levels_result *result)
{
    if(!user_params || !result)
        return -1;

    return __system_2_levels_eval_general(user_params,result,system_2_levels_adiabatic_f,system_2_levels_adiabatic_df,system_2_levels_adiabatic_fdf);
}