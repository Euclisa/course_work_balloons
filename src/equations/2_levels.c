#include <equations/2_levels.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>
#include <gsl/gsl_sf_trig.h>
#include <gsl/gsl_blas.h>

const size_t N_eq = 24;

const int N_top_PpA = 5;    // Points per Arc
const int N_top = N_top_PpA*3;
const int N_bot_PpA = 5;    // Points per Arc
const int N_bot = N_bot_PpA*2;


int system_2_levels_f(const gsl_vector *x, void *p, gsl_vector *f)
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
    // p_top*r_ad*gsl_sf_sin(a_ad+phi_ad) - (p_top-p_bot)*r_dc*gsl_sf_sin(a_dc) - p_bot*r_ed*gsl_sf_sin(a_ed-phi_ed);

    gsl_matrix_set(J,17,0,-p_top*r_ad*gsl_sf_cos(a_ad+phi_ad));
    gsl_matrix_set(J,17,1,-p_top*gsl_sf_sin(a_ad+phi_ad));
    gsl_matrix_set(J,17,4,-p_top*r_ad*gsl_sf_cos(a_ad+phi_ad));
    gsl_matrix_set(J,17,11,(p_top-p_bot)*gsl_sf_sin(a_dc));
    gsl_matrix_set(J,17,14,(p_top-p_bot)*r_dc*gsl_sf_cos(a_dc));
    gsl_matrix_set(J,17,15,-p_bot*r_ed*gsl_sf_cos(a_ed-phi_ed));
    gsl_matrix_set(J,17,16,p_bot*gsl_sf_sin(a_ed-phi_ed));
    gsl_matrix_set(J,17,22,-r_ad*gsl_sf_sin(a_ad+phi_ad) + r_dc*gsl_sf_sin(a_dc));
    gsl_matrix_set(J,17,23,-r_dc*gsl_sf_sin(a_dc) + r_ed*gsl_sf_sin(a_ed-phi_ed));

    // -p_top*r_ad*gsl_sf_cos(a_ad+phi_ad) + (p_top-p_bot)*r_dc*gsl_sf_cos(a_dc) + p_bot*r_ed*gsl_sf_cos(a_ed-phi_ed)
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

    //p_top*r_cb*gsl_sf_sin(a_cb) - (p_top-p_bot)*r_dc*gsl_sf_sin(a_dc+phi_dc) - p_bot*r_ec*gsl_sf_sin(a_ec+phi_ec)
    gsl_matrix_set(J,19,6,(p_top-params->p_ac)*gsl_sf_sin(a_cb));
    gsl_matrix_set(J,19,9,(p_top-params->p_ac)*r_cb*gsl_sf_cos(a_cb));
    gsl_matrix_set(J,19,10,-(p_top-p_bot)*r_dc*gsl_sf_cos(a_dc+phi_dc));
    gsl_matrix_set(J,19,11,-(p_top-p_bot)*gsl_sf_sin(a_dc+phi_dc));
    gsl_matrix_set(J,19,14,-(p_top-p_bot)*r_dc*gsl_sf_cos(a_dc+phi_dc));
    gsl_matrix_set(J,19,18,-(p_bot-params->p_ac)*r_ec*gsl_sf_cos(a_ec+phi_ec));
    gsl_matrix_set(J,19,19,-(p_bot-params->p_ac)*gsl_sf_sin(a_ec+phi_ec));
    gsl_matrix_set(J,19,22,r_cb*gsl_sf_sin(a_cb) - r_dc*gsl_sf_sin(a_dc+phi_dc));
    gsl_matrix_set(J,19,23,r_dc*gsl_sf_sin(a_dc+phi_dc) - r_ec*gsl_sf_sin(a_ec+phi_ec));

    //-p_top*r_cb*gsl_sf_cos(a_cb) + (p_top-p_bot)*r_dc*gsl_sf_cos(a_dc+phi_dc) + p_bot*r_ec*gsl_sf_cos(a_ec+phi_ec)
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


static inline gsl_vector *vector_centred(const gsl_vector *v, const gsl_vector *center)
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


double area_between_vectors_triangle(const gsl_vector *a, const gsl_vector *b)
{
    double a1 = gsl_vector_get(a,0);
    double a2 = gsl_vector_get(a,1);
    double b1 = gsl_vector_get(b,0);
    double b2 = gsl_vector_get(b,1);
    double area = (a1*b2 - a2*b1)/2;
    area = fabs(area);

    return area;
}


static inline double area_segment(double ang, double r)
{
    return ang/2 * gsl_pow_2(r); 
}


int point_from_alpha(double alpha, double norm, gsl_vector *center, gsl_vector *v)
{
    gsl_vector_set(v,0,-gsl_sf_cos(alpha));
    gsl_vector_set(v,1,gsl_sf_sin(alpha));
    gsl_vector_scale(v,norm);
    gsl_vector_add(v,center);

    return GSL_SUCCESS;
}


static inline double add_angs(double ang1, double ang2)
{
    double res = ang1 + ang2;
    double w = res/M_PI/2;
    double w0 = (int)w;
    w -= w0;
    res = w * 2*M_PI;

    return res;
}


static inline double sub_angs(double ang1, double ang2)
{
    double res = ang1 - ang2;
    double w = res/M_PI/2;
    double w0 = (int)w;
    w -= w0;
    res = w * 2*M_PI;
    if(res < 0)
        res += 2*M_PI;

    return res;
}


int compute_init_levels_positions(const gsl_vector *A, const gsl_vector *B, double r_top, double r_bot, double phi_ad, double phi_dc, struct system_2_levels_params *params, gsl_vector *x0)
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

    double a_dc = add_angs(a_ad,phi_ad);
    double a_cb = add_angs(a_dc,phi_dc);

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

    gsl_vector *center_bot = gsl_vector_alloc(2);
    if(center_from_points_and_radius(C,D,r_bot,center_bot) != GSL_SUCCESS)
    {
        fprintf(stderr,"Failed to compute bottom level center point\n");
        exit(1);
    }

    gsl_vector *E = gsl_vector_alloc(2);
    gsl_vector_set(E,0,gsl_vector_get(center_bot,0));
    gsl_vector_set(E,1,gsl_vector_get(center_bot,1)-r_bot);
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
    S_bot += area_segment(phi_dc_bot,r_bot);
    S_top += area_between_vectors_triangle(cC_top,cD_top);
    S_top += area_between_vectors_triangle(cA,cB);
    S_top += area_segment(phi_cb,r_top);
    S_top += area_segment(phi_ad,r_top);

    params->Ax = gsl_vector_get(A,0);
    params->Ay = gsl_vector_get(A,1);
    params->Bx = gsl_vector_get(B,0);
    params->By = gsl_vector_get(B,1);
    params->phi_ad_0 = phi_ad;
    params->phi_cb_0 = phi_cb;
    params->phi_dc_0 = phi_dc;
    params->phi_ec_0 = phi_ec;
    params->phi_ed_0 = phi_ed;
    params->r_bot_0 = r_bot;
    params->r_top_0 = r_top;
    params->S_bot_0 = S_bot;
    params->S_top_0 = S_top;

    double x_center_top = gsl_vector_get(center_top,0);
    double y_center_top = gsl_vector_get(center_top,1);
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
    gsl_vector_set(x0,15,phi_ed);
    gsl_vector_set(x0,16,r_bot);
    gsl_vector_set(x0,17,y_center_bot);
    gsl_vector_set(x0,18,phi_ec);
    gsl_vector_set(x0,19,r_bot);
    gsl_vector_set(x0,20,y_center_bot);
    gsl_vector_set(x0,21,x_center_bot);

    gsl_vector_free(cA);
    gsl_vector_free(cB);
    gsl_vector_free(cC_bot);
    gsl_vector_free(cC_top);
    gsl_vector_free(C);
    gsl_vector_free(cD_bot);
    gsl_vector_free(cD_top);
    gsl_vector_free(D);
    gsl_vector_free(cE);
    gsl_vector_free(E);
    gsl_vector_free(uX);
    gsl_vector_free(center_bot);
    gsl_vector_free(center_top);

    return GSL_SUCCESS;
}


int compute_init_config(double Ax, double Ay, double Bx, double By, double phi_ad, double phi_dc, double r_top, double r_bot, double k, double p_top, double p_bot, double p, double p_ac, struct system_2_levels_params *params, gsl_vector *x0)
{
    gsl_vector *A = gsl_vector_alloc(2);
    gsl_vector_set(A,0,Ax);
    gsl_vector_set(A,1,Ay);
    gsl_vector *B = gsl_vector_alloc(2);
    gsl_vector_set(B,0,Bx);
    gsl_vector_set(B,1,By);
    if(compute_init_levels_positions(A,B,r_top,r_bot,phi_ad,phi_dc,params,x0) != GSL_SUCCESS)
    {
        fprintf(stderr,"Failed to compute initial levels position\n");
        exit(1);
    }
    params->k = k;
    params->p = p;
    params->p_top_0 = p_top;
    params->p_bot_0 = p_bot;
    params->p_ac = p_ac;

    gsl_vector_set(x0,22,p_top);
    gsl_vector_set(x0,23,p_bot);

    gsl_vector_free(A);
    gsl_vector_free(B);

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


void system_2_levels_J_estimate(const gsl_vector *x0, struct system_2_levels_params *params, gsl_matrix *J_est)
{
    gsl_vector *f = gsl_vector_alloc(N_eq);
    gsl_matrix *J = gsl_matrix_alloc(N_eq,N_eq);
    system_2_levels_fdf(x0,params,f,J);
    gsl_vector *x0_new = gsl_vector_alloc(N_eq);
    gsl_vector_memcpy(x0_new,x0);
    const double eps = 1e-12;
    for(int i = 0; i < N_eq; ++i)
    {
        gsl_vector df_dxi = gsl_matrix_column(J_est,i).vector;
        double old_xi_v = gsl_vector_get(x0,i);
        gsl_vector_set(x0_new,i,old_xi_v+eps);
        system_2_levels_f(x0_new,params,&df_dxi);
        gsl_vector_sub(&df_dxi,f);
        gsl_vector_scale(&df_dxi,1.0/eps);
        gsl_vector_set(x0_new,i,old_xi_v);
    }
}


static inline void system_2_levels_print_J_diff(FILE *stream, const gsl_vector *x, struct system_2_levels_params *params)
{
    gsl_matrix *J_est = gsl_matrix_alloc(N_eq,N_eq);
    gsl_matrix *J_diff = gsl_matrix_alloc(N_eq,N_eq);
    gsl_matrix *J = gsl_matrix_alloc(N_eq,N_eq);
    system_2_levels_df(x,params,J);
    system_2_levels_J_estimate(x,params,J_est);
    gsl_matrix_memcpy(J_diff,J);
    gsl_matrix_sub(J_diff,J_est);
    fprintf(stream,"Jacobian true:\n");
    print_matrix(stream,J);
    fprintf(stream,"Jacobian est:\n");
    print_matrix(stream,J_est);
    fprintf(stream,"Jacobian diff:\n");
    print_matrix(stream,J_diff);
}


int system_2_levels_eval()
{
    printf("Enter\n");
    struct system_2_levels_params params;
    gsl_vector *x0 = gsl_vector_alloc(N_eq);
    compute_init_config(0.482,1.4,0.28,0.85,3.129,1.162,0.5,0.35,1,20000,6500,101325,1500,&params,x0);
    printf("Enter1\n");

    const gsl_multiroot_fdfsolver_type *T = gsl_multiroot_fdfsolver_newton;
    gsl_multiroot_fdfsolver *s = gsl_multiroot_fdfsolver_alloc(T,N_eq);

    gsl_multiroot_function_fdf fdf;
    fdf.f = system_2_levels_f;
    fdf.df = system_2_levels_df;
    fdf.fdf = system_2_levels_fdf;
    fdf.n = N_eq;
    fdf.params = &params;

    FILE *fx0 = fopen("2_levels_init.txt","w");
    printf("Enter2\n");
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

    return GSL_SUCCESS;
}
