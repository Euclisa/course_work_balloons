#include <equations/2_levels.h>
#include <equations/3_levels.h>
#include <gsl/gsl_multiroots.h>
#include <gtk-3.0/gtk/gtk.h>
#include <pthread.h>
#include <stdbool.h>
#include <stdio.h>


struct adiabatic_mode_widgets
{
    GtkLabel *adiabatic_constant_label;
    GtkSpinButton *adiabatic_constant_spin;
};

struct app_level_context_l2
{
    pthread_mutex_t params_lock;
    pthread_cond_t params_cond;
    struct system_2_levels_user_params user_params;
    bool params_dirty;
    bool adiabatic;
    pthread_t drawing_thread;
    pthread_mutex_t result_lock;
    struct system_2_levels_result result;
    struct adiabatic_mode_widgets adia_widgets;
};

struct app_level_context_l3
{
    pthread_mutex_t params_lock;
    pthread_cond_t params_cond;
    struct system_3_levels_user_params user_params;
    bool params_dirty;
    bool adiabatic;
    pthread_t drawing_thread;
    pthread_mutex_t result_lock;
    struct system_3_levels_result result;
    struct adiabatic_mode_widgets adia_widgets;
};

struct app_level_context_l2 l2_context;
struct app_level_context_l3 l3_context;


static void draw_function_l2(GtkDrawingArea *area, cairo_t *cr, gpointer data)
{
    pthread_mutex_lock(&l2_context.result_lock);
    
    const int width = gtk_widget_get_allocated_width(GTK_WIDGET(area));
    const int height = gtk_widget_get_allocated_height(GTK_WIDGET(area));

    const int pixels_per_meter = MIN(width,height)/2;

    double r_ad = l2_context.result.r_ad*pixels_per_meter;
    double r_cb = l2_context.result.r_cb*pixels_per_meter;
    double r_dc = l2_context.result.r_dc*pixels_per_meter;
    double r_ec = l2_context.result.r_ec*pixels_per_meter;
    double r_ed = l2_context.result.r_ed*pixels_per_meter;
    double x_ad = l2_context.result.x_ad*pixels_per_meter;
    double y_ad = l2_context.result.y_ad*pixels_per_meter;
    double x_cb = l2_context.result.x_cb*pixels_per_meter;
    double y_cb = l2_context.result.y_cb*pixels_per_meter;
    double x_dc = l2_context.result.x_dc*pixels_per_meter;
    double y_dc = l2_context.result.y_dc*pixels_per_meter;
    double y_ec = l2_context.result.y_ec*pixels_per_meter;
    double y_ed = l2_context.result.y_ed*pixels_per_meter;
    double x_bot = l2_context.result.x_bot*pixels_per_meter;

    cairo_translate(cr,0,height);
    cairo_scale(cr,1,-1);

    double a_ad = G_PI-l2_context.result.a_ad;
    cairo_arc(cr,x_ad,y_ad,r_ad,a_ad-l2_context.result.phi_ad,a_ad);
    cairo_stroke(cr);

    double a_dc = G_PI-l2_context.result.a_dc;
    cairo_arc(cr,x_dc,y_dc,r_dc,a_dc-l2_context.result.phi_dc,a_dc);
    cairo_stroke(cr);
    
    double a_cb = G_PI-l2_context.result.a_cb;
    cairo_arc(cr,x_cb,y_cb,r_cb,a_cb-l2_context.result.phi_cb,a_cb);
    cairo_stroke(cr);
    
    double a_ec = -G_PI_2;
    cairo_arc(cr,x_bot,y_ec,r_ec,a_ec-l2_context.result.phi_ec,a_ec);
    cairo_stroke(cr);

    double a_ed = -G_PI_2;
    cairo_arc(cr,x_bot,y_ed,r_ed,a_ed,a_ed+l2_context.result.phi_ed);
    cairo_stroke(cr);
    
    
    cairo_select_font_face(cr,"monospace",CAIRO_FONT_SLANT_NORMAL,CAIRO_FONT_WEIGHT_NORMAL);
    cairo_set_font_size(cr,20);
    cairo_scale(cr,1,-1);

    double letter_offset = height / 20;

    double x_a = x_ad + r_ad*cos(a_ad) - letter_offset;
    double y_a = -y_ad - r_ad*sin(a_ad);
    cairo_move_to(cr,x_a,y_a);
    cairo_show_text(cr,"A");

    double x_b = x_cb + r_cb*cos(a_cb-l2_context.result.phi_cb) - letter_offset;
    double y_b = -y_cb - r_cb*sin(a_cb-l2_context.result.phi_cb);
    cairo_move_to(cr,x_b,y_b);
    cairo_show_text(cr,"B");

    double x_c = x_cb + r_cb*cos(a_cb) - letter_offset;
    double y_c = -y_cb - r_cb*sin(a_cb) + letter_offset/2;
    cairo_move_to(cr,x_c,y_c);
    cairo_show_text(cr,"C");

    double x_d = x_ad + r_ad*cos(a_ad+l2_context.result.phi_ad) + letter_offset/4;
    double y_d = -y_ad - r_ad*sin(a_ad+l2_context.result.phi_ad) + letter_offset;
    cairo_move_to(cr,x_d,y_d);
    cairo_show_text(cr,"D");

    double x_e = x_bot + r_ec*cos(a_ec);
    double y_e = -y_ec - r_ec*sin(a_ec) + letter_offset;
    cairo_move_to(cr,x_e,y_e);
    cairo_show_text(cr,"E");

    pthread_mutex_unlock(&l2_context.result_lock);
}


static void draw_function_l3(GtkDrawingArea *area, cairo_t *cr, gpointer data)
{
    pthread_mutex_lock(&l3_context.result_lock);
    
    const int width = gtk_widget_get_allocated_width(GTK_WIDGET(area));
    const int height = gtk_widget_get_allocated_height(GTK_WIDGET(area));

    const int pixels_per_meter = MIN(width,height)/3;

    double r_ad = l3_context.result.r_ad*pixels_per_meter;
    double r_cb = l3_context.result.r_cb*pixels_per_meter;
    double r_dc = l3_context.result.r_dc*pixels_per_meter;
    double r_df = l3_context.result.r_df*pixels_per_meter;
    double r_ec = l3_context.result.r_ec*pixels_per_meter;
    double r_fe = l3_context.result.r_fe*pixels_per_meter;
    double r_ge = l3_context.result.r_ge*pixels_per_meter;
    double r_gf = l3_context.result.r_gf*pixels_per_meter;
    double x_ad = l3_context.result.x_ad*pixels_per_meter;
    double y_ad = l3_context.result.y_ad*pixels_per_meter;
    double x_cb = l3_context.result.x_cb*pixels_per_meter;
    double y_cb = l3_context.result.y_cb*pixels_per_meter;
    double x_dc = l3_context.result.x_dc*pixels_per_meter;
    double y_dc = l3_context.result.y_dc*pixels_per_meter;
    double x_df = l3_context.result.x_df*pixels_per_meter;
    double y_df = l3_context.result.y_df*pixels_per_meter;
    double x_ec = l3_context.result.x_ec*pixels_per_meter;
    double y_ec = l3_context.result.y_ec*pixels_per_meter;
    double x_fe = l3_context.result.x_fe*pixels_per_meter;
    double y_fe = l3_context.result.y_fe*pixels_per_meter;
    double y_ge = l3_context.result.y_ge*pixels_per_meter;
    double y_gf = l3_context.result.y_gf*pixels_per_meter;
    double x_bot = l3_context.result.x_bot*pixels_per_meter;

    cairo_translate(cr,0,height);
    cairo_scale(cr,1,-1);

    double a_ad = G_PI-l3_context.result.a_ad;
    cairo_arc(cr,x_ad,y_ad,r_ad,a_ad-l3_context.result.phi_ad,a_ad);
    cairo_stroke(cr);

    double a_dc = G_PI-l3_context.result.a_dc;
    cairo_arc(cr,x_dc,y_dc,r_dc,a_dc-l3_context.result.phi_dc,a_dc);
    cairo_stroke(cr);
    
    double a_cb = G_PI-l3_context.result.a_cb;
    cairo_arc(cr,x_cb,y_cb,r_cb,a_cb-l3_context.result.phi_cb,a_cb);
    cairo_stroke(cr);

    double a_df = G_PI-l3_context.result.a_df;
    cairo_arc(cr,x_df,y_df,r_df,a_df-l3_context.result.phi_df,a_df);
    cairo_stroke(cr);

    double a_ec = G_PI-l3_context.result.a_ec;
    cairo_arc(cr,x_ec,y_ec,r_ec,a_ec-l3_context.result.phi_ec,a_ec);
    cairo_stroke(cr);

    double a_fe = G_PI-l3_context.result.a_fe;
    cairo_arc(cr,x_fe,y_fe,r_fe,a_fe-l3_context.result.phi_fe,a_fe);
    cairo_stroke(cr);
    
    double a_ge = -G_PI_2;
    cairo_arc(cr,x_bot,y_ge,r_ge,a_ge-l3_context.result.phi_ge,a_ge);
    cairo_stroke(cr);

    double a_gf = -G_PI_2;
    cairo_arc(cr,x_bot,y_gf,r_gf,a_gf,a_gf+l3_context.result.phi_gf);
    cairo_stroke(cr);
    
    
    cairo_select_font_face(cr,"monospace",CAIRO_FONT_SLANT_NORMAL,CAIRO_FONT_WEIGHT_NORMAL);
    cairo_set_font_size(cr,20);
    cairo_scale(cr,1,-1);

    double letter_offset = height / 20;

    double x_a = x_ad + r_ad*cos(a_ad) - letter_offset;
    double y_a = -y_ad - r_ad*sin(a_ad);
    cairo_move_to(cr,x_a,y_a);
    cairo_show_text(cr,"A");

    double x_b = x_cb + r_cb*cos(a_cb-l2_context.result.phi_cb) - letter_offset;
    double y_b = -y_cb - r_cb*sin(a_cb-l2_context.result.phi_cb);
    cairo_move_to(cr,x_b,y_b);
    cairo_show_text(cr,"B");

    double x_c = x_cb + r_cb*cos(a_cb) - letter_offset;
    double y_c = -y_cb - r_cb*sin(a_cb) + letter_offset/2;
    cairo_move_to(cr,x_c,y_c);
    cairo_show_text(cr,"C");

    double x_d = x_ad + r_ad*cos(a_ad+l3_context.result.phi_ad) + letter_offset/4;
    double y_d = -y_ad - r_ad*sin(a_ad+l3_context.result.phi_ad) + letter_offset;
    cairo_move_to(cr,x_d,y_d);
    cairo_show_text(cr,"D");

    double x_e = x_ec + r_ec*cos(a_ec) - letter_offset;
    double y_e = -y_ec - r_ec*sin(a_ec) + letter_offset;
    cairo_move_to(cr,x_e,y_e);
    cairo_show_text(cr,"E");

    double x_f = x_bot + r_gf*cos(a_gf+l3_context.result.phi_gf) + letter_offset;
    double y_f = -y_gf - r_gf*sin(a_gf+l3_context.result.phi_gf) - letter_offset/4;
    cairo_move_to(cr,x_f,y_f);
    cairo_show_text(cr,"F");

    double x_g = x_bot + r_gf*cos(a_gf);
    double y_g = -y_gf - r_gf*sin(a_gf) + letter_offset;
    cairo_move_to(cr,x_g,y_g);
    cairo_show_text(cr,"G");


    pthread_mutex_unlock(&l3_context.result_lock);
}


static void queue_update_picture_l2(GtkDrawingArea *area, const struct system_2_levels_user_params *params_extracted, bool adiabatic_extracted)
{
    struct system_2_levels_result result_local;

    if(adiabatic_extracted)
        system_2_levels_adiabatic_eval(params_extracted,&result_local);
    else
        system_2_levels_eval(params_extracted,&result_local);

    pthread_mutex_lock(&l2_context.result_lock);
    memcpy(&l2_context.result,&result_local,sizeof(struct system_2_levels_result));
    gtk_widget_queue_draw(GTK_WIDGET(area));
    pthread_mutex_unlock(&l2_context.result_lock);
}

static void queue_update_picture_l3(GtkDrawingArea *area, const struct system_3_levels_user_params *params_extracted)
{
    struct system_3_levels_result result_local;
    system_3_levels_eval(params_extracted,&result_local);

    pthread_mutex_lock(&l3_context.result_lock);
    memcpy(&l3_context.result,&result_local,sizeof(struct system_3_levels_result));
    gtk_widget_queue_draw(GTK_WIDGET(area));
    pthread_mutex_unlock(&l3_context.result_lock);
}


static void *update_picture_l2(GtkDrawingArea *area)
{
    while(true)
    {
        pthread_mutex_lock(&l2_context.params_lock);
        if(!l2_context.params_dirty)
            pthread_cond_wait(&l2_context.params_cond,&l2_context.params_lock);
        
        struct system_2_levels_user_params params_extracted;
        memcpy(&params_extracted,&l2_context.user_params,sizeof(struct system_2_levels_user_params));
        bool adiabatic_extracted = l2_context.adiabatic;
        l2_context.params_dirty = false;

        pthread_mutex_unlock(&l2_context.params_lock);

        queue_update_picture_l2(area,&params_extracted,adiabatic_extracted);
    }

    return NULL;
}


static void *update_picture_l3(GtkDrawingArea *area)
{
    while(true)
    {
        pthread_mutex_lock(&l3_context.params_lock);
        if(!l3_context.params_dirty)
            pthread_cond_wait(&l3_context.params_cond,&l3_context.params_lock);
        
        struct system_3_levels_user_params params_extracted;
        memcpy(&params_extracted,&l3_context.user_params,sizeof(struct system_3_levels_user_params));
        l3_context.params_dirty = false;

        pthread_mutex_unlock(&l3_context.params_lock);

        queue_update_picture_l3(area,&params_extracted);
    }

    return NULL;
}


static void spin_button_value_changed_cb_l2(GtkSpinButton *spin_button, double *param)
{
    pthread_mutex_lock(&l2_context.params_lock);
    double value = (double)gtk_spin_button_get_value(spin_button);
    *param = value;
    l2_context.params_dirty = true;
    pthread_cond_signal(&l2_context.params_cond);
    pthread_mutex_unlock(&l2_context.params_lock);
}

static void spin_button_value_changed_cb_l3(GtkSpinButton *spin_button, double *param)
{
    pthread_mutex_lock(&l3_context.params_lock);
    double value = (double)gtk_spin_button_get_value(spin_button);
    *param = value;
    l3_context.params_dirty = true;
    pthread_cond_signal(&l3_context.params_cond);
    pthread_mutex_unlock(&l3_context.params_lock);
}

static void adiabatic_mode_changed_cb_l2(GtkSwitch *adiabatic_sw, gpointer *data)
{
    pthread_mutex_lock(&l2_context.params_lock);
    l2_context.adiabatic = (bool)gtk_switch_get_active(adiabatic_sw);
    gtk_widget_set_sensitive(GTK_WIDGET(l2_context.adia_widgets.adiabatic_constant_label),l2_context.adiabatic);
    gtk_widget_set_sensitive(GTK_WIDGET(l2_context.adia_widgets.adiabatic_constant_spin),l2_context.adiabatic);

    l2_context.params_dirty = true;
    pthread_cond_signal(&l2_context.params_cond);
    pthread_mutex_unlock(&l2_context.params_lock);
}


static void init_widgets_l2(GtkBuilder *builder)
{
    GObject *ax_spin = gtk_builder_get_object(builder,"ax_l2");
    g_signal_connect(ax_spin,"value-changed",G_CALLBACK(spin_button_value_changed_cb_l2),&l2_context.user_params.Ax);
    spin_button_value_changed_cb_l2(GTK_SPIN_BUTTON(ax_spin),&l2_context.user_params.Ax);

    GObject *ay_spin = gtk_builder_get_object(builder,"ay_l2");
    g_signal_connect(ay_spin,"value-changed",G_CALLBACK(spin_button_value_changed_cb_l2),&l2_context.user_params.Ay);
    spin_button_value_changed_cb_l2(GTK_SPIN_BUTTON(ay_spin),&l2_context.user_params.Ay);

    GObject *bx_spin = gtk_builder_get_object(builder,"bx_l2");
    g_signal_connect(bx_spin,"value-changed",G_CALLBACK(spin_button_value_changed_cb_l2),&l2_context.user_params.Bx);
    spin_button_value_changed_cb_l2(GTK_SPIN_BUTTON(bx_spin),&l2_context.user_params.Bx);

    GObject *by_spin = gtk_builder_get_object(builder,"by_l2");
    g_signal_connect(by_spin,"value-changed",G_CALLBACK(spin_button_value_changed_cb_l2),&l2_context.user_params.By);
    spin_button_value_changed_cb_l2(GTK_SPIN_BUTTON(by_spin),&l2_context.user_params.By);

    GObject *phi_ad_spin = gtk_builder_get_object(builder,"phi_ad_l2");
    g_signal_connect(phi_ad_spin,"value-changed",G_CALLBACK(spin_button_value_changed_cb_l2),&l2_context.user_params.phi_ad_0);
    spin_button_value_changed_cb_l2(GTK_SPIN_BUTTON(phi_ad_spin),&l2_context.user_params.phi_ad_0);

    GObject *phi_dc_spin = gtk_builder_get_object(builder,"phi_dc_l2");
    g_signal_connect(phi_dc_spin,"value-changed",G_CALLBACK(spin_button_value_changed_cb_l2),&l2_context.user_params.phi_dc_0);
    spin_button_value_changed_cb_l2(GTK_SPIN_BUTTON(phi_dc_spin),&l2_context.user_params.phi_dc_0);

    GObject *r_top_spin = gtk_builder_get_object(builder,"r_top_l2");
    g_signal_connect(r_top_spin,"value-changed",G_CALLBACK(spin_button_value_changed_cb_l2),&l2_context.user_params.r_top_0);
    spin_button_value_changed_cb_l2(GTK_SPIN_BUTTON(r_top_spin),&l2_context.user_params.r_top_0);

    GObject *r_bot_spin = gtk_builder_get_object(builder,"r_bot_l2");
    g_signal_connect(r_bot_spin,"value-changed",G_CALLBACK(spin_button_value_changed_cb_l2),&l2_context.user_params.r_bot_0);
    spin_button_value_changed_cb_l2(GTK_SPIN_BUTTON(r_bot_spin),&l2_context.user_params.r_bot_0);

    GObject *p_top_spin = gtk_builder_get_object(builder,"p_top_l2");
    g_signal_connect(p_top_spin,"value-changed",G_CALLBACK(spin_button_value_changed_cb_l2),&l2_context.user_params.p_top_0);
    spin_button_value_changed_cb_l2(GTK_SPIN_BUTTON(p_top_spin),&l2_context.user_params.p_top_0);

    GObject *p_bot_spin = gtk_builder_get_object(builder,"p_bot_l2");
    g_signal_connect(p_bot_spin,"value-changed",G_CALLBACK(spin_button_value_changed_cb_l2),&l2_context.user_params.p_bot_0);
    spin_button_value_changed_cb_l2(GTK_SPIN_BUTTON(p_bot_spin),&l2_context.user_params.p_bot_0);

    GObject *p_ship_spin = gtk_builder_get_object(builder,"p_ship_l2");
    g_signal_connect(p_ship_spin,"value-changed",G_CALLBACK(spin_button_value_changed_cb_l2),&l2_context.user_params.p_ac);
    spin_button_value_changed_cb_l2(GTK_SPIN_BUTTON(p_ship_spin),&l2_context.user_params.p_ac);

    GObject *p_atm_spin = gtk_builder_get_object(builder,"p_atm_l2");
    g_signal_connect(p_atm_spin,"value-changed",G_CALLBACK(spin_button_value_changed_cb_l2),&l2_context.user_params.p_atm);
    spin_button_value_changed_cb_l2(GTK_SPIN_BUTTON(p_atm_spin),&l2_context.user_params.p_atm);


    l2_context.adia_widgets.adiabatic_constant_label = GTK_LABEL(gtk_builder_get_object(builder,"adiabatic_constant_label"));
    l2_context.adia_widgets.adiabatic_constant_spin = GTK_SPIN_BUTTON(gtk_builder_get_object(builder,"adiabatic_constant_spin"));
    g_signal_connect(l2_context.adia_widgets.adiabatic_constant_spin,"value-changed",G_CALLBACK(spin_button_value_changed_cb_l2),&l2_context.user_params.k);
    spin_button_value_changed_cb_l2(GTK_SPIN_BUTTON(l2_context.adia_widgets.adiabatic_constant_spin),&l2_context.user_params.k);

    GObject *adiabatic_sw = gtk_builder_get_object(builder,"adiabatic_l2");
    g_signal_connect(adiabatic_sw,"state-set",G_CALLBACK(adiabatic_mode_changed_cb_l2),NULL);
    adiabatic_mode_changed_cb_l2(GTK_SWITCH(adiabatic_sw),NULL);

    GtkDrawingArea *area = GTK_DRAWING_AREA(gtk_builder_get_object(builder,"drawing_area_l2"));
    g_signal_connect(G_OBJECT(area),"draw",G_CALLBACK(draw_function_l2),NULL);

    l2_context.params_dirty = true;

    pthread_create(&l2_context.drawing_thread,NULL,(void *(*)(void*))update_picture_l2,area);
}


static void init_widgets_l3(GtkBuilder *builder)
{
    GObject *ax_spin = gtk_builder_get_object(builder,"ax_l3");
    g_signal_connect(ax_spin,"value-changed",G_CALLBACK(spin_button_value_changed_cb_l3),&l3_context.user_params.Ax);
    spin_button_value_changed_cb_l3(GTK_SPIN_BUTTON(ax_spin),&l3_context.user_params.Ax);

    GObject *ay_spin = gtk_builder_get_object(builder,"ay_l3");
    g_signal_connect(ay_spin,"value-changed",G_CALLBACK(spin_button_value_changed_cb_l3),&l3_context.user_params.Ay);
    spin_button_value_changed_cb_l3(GTK_SPIN_BUTTON(ay_spin),&l3_context.user_params.Ay);

    GObject *bx_spin = gtk_builder_get_object(builder,"bx_l3");
    g_signal_connect(bx_spin,"value-changed",G_CALLBACK(spin_button_value_changed_cb_l3),&l3_context.user_params.Bx);
    spin_button_value_changed_cb_l3(GTK_SPIN_BUTTON(bx_spin),&l3_context.user_params.Bx);

    GObject *by_spin = gtk_builder_get_object(builder,"by_l3");
    g_signal_connect(by_spin,"value-changed",G_CALLBACK(spin_button_value_changed_cb_l3),&l3_context.user_params.By);
    spin_button_value_changed_cb_l3(GTK_SPIN_BUTTON(by_spin),&l3_context.user_params.By);

    GObject *phi_ad_spin = gtk_builder_get_object(builder,"phi_ad_l3");
    g_signal_connect(phi_ad_spin,"value-changed",G_CALLBACK(spin_button_value_changed_cb_l3),&l3_context.user_params.phi_ad_0);
    spin_button_value_changed_cb_l3(GTK_SPIN_BUTTON(phi_ad_spin),&l3_context.user_params.phi_ad_0);

    GObject *phi_dc_spin = gtk_builder_get_object(builder,"phi_dc_l3");
    g_signal_connect(phi_dc_spin,"value-changed",G_CALLBACK(spin_button_value_changed_cb_l3),&l3_context.user_params.phi_dc_0);
    spin_button_value_changed_cb_l3(GTK_SPIN_BUTTON(phi_dc_spin),&l3_context.user_params.phi_dc_0);

    GObject *phi_df_spin = gtk_builder_get_object(builder,"phi_df_l3");
    g_signal_connect(phi_df_spin,"value-changed",G_CALLBACK(spin_button_value_changed_cb_l3),&l3_context.user_params.phi_df_0);
    spin_button_value_changed_cb_l3(GTK_SPIN_BUTTON(phi_df_spin),&l3_context.user_params.phi_df_0);

    GObject *phi_fe_spin = gtk_builder_get_object(builder,"phi_fe_l3");
    g_signal_connect(phi_fe_spin,"value-changed",G_CALLBACK(spin_button_value_changed_cb_l3),&l3_context.user_params.phi_fe_0);
    spin_button_value_changed_cb_l3(GTK_SPIN_BUTTON(phi_fe_spin),&l3_context.user_params.phi_fe_0);

    GObject *r_top_spin = gtk_builder_get_object(builder,"r_top_l3");
    g_signal_connect(r_top_spin,"value-changed",G_CALLBACK(spin_button_value_changed_cb_l3),&l3_context.user_params.r_top_0);
    spin_button_value_changed_cb_l3(GTK_SPIN_BUTTON(r_top_spin),&l3_context.user_params.r_top_0);

    GObject *r_mid_spin = gtk_builder_get_object(builder,"r_mid_l3");
    g_signal_connect(r_mid_spin,"value-changed",G_CALLBACK(spin_button_value_changed_cb_l3),&l3_context.user_params.r_mid_0);
    spin_button_value_changed_cb_l3(GTK_SPIN_BUTTON(r_mid_spin),&l3_context.user_params.r_mid_0);

    GObject *r_bot_spin = gtk_builder_get_object(builder,"r_bot_l3");
    g_signal_connect(r_bot_spin,"value-changed",G_CALLBACK(spin_button_value_changed_cb_l3),&l3_context.user_params.r_bot_0);
    spin_button_value_changed_cb_l3(GTK_SPIN_BUTTON(r_bot_spin),&l3_context.user_params.r_bot_0);

    GObject *p_top_spin = gtk_builder_get_object(builder,"p_top_l3");
    g_signal_connect(p_top_spin,"value-changed",G_CALLBACK(spin_button_value_changed_cb_l3),&l3_context.user_params.p_top_0);
    spin_button_value_changed_cb_l3(GTK_SPIN_BUTTON(p_top_spin),&l3_context.user_params.p_top_0);

    GObject *p_mid_spin = gtk_builder_get_object(builder,"p_mid_l3");
    g_signal_connect(p_mid_spin,"value-changed",G_CALLBACK(spin_button_value_changed_cb_l3),&l3_context.user_params.p_mid_0);
    spin_button_value_changed_cb_l3(GTK_SPIN_BUTTON(p_mid_spin),&l3_context.user_params.p_mid_0);

    GObject *p_bot_spin = gtk_builder_get_object(builder,"p_bot_l3");
    g_signal_connect(p_bot_spin,"value-changed",G_CALLBACK(spin_button_value_changed_cb_l3),&l3_context.user_params.p_bot_0);
    spin_button_value_changed_cb_l3(GTK_SPIN_BUTTON(p_bot_spin),&l3_context.user_params.p_bot_0);

    GObject *p_ship_spin = gtk_builder_get_object(builder,"p_ship_l3");
    g_signal_connect(p_ship_spin,"value-changed",G_CALLBACK(spin_button_value_changed_cb_l3),&l3_context.user_params.p_ac);
    spin_button_value_changed_cb_l3(GTK_SPIN_BUTTON(p_ship_spin),&l3_context.user_params.p_ac);

    GObject *p_atm_spin = gtk_builder_get_object(builder,"p_atm_l3");
    g_signal_connect(p_atm_spin,"value-changed",G_CALLBACK(spin_button_value_changed_cb_l3),&l3_context.user_params.p_atm);
    spin_button_value_changed_cb_l3(GTK_SPIN_BUTTON(p_atm_spin),&l3_context.user_params.p_atm);

    GtkDrawingArea *area = GTK_DRAWING_AREA(gtk_builder_get_object(builder,"drawing_area_l3"));
    g_signal_connect(G_OBJECT(area),"draw",G_CALLBACK(draw_function_l3),NULL);

    l3_context.params_dirty = true;

    pthread_create(&l3_context.drawing_thread,NULL,(void *(*)(void*))update_picture_l3,area);
}


static void activate(GtkApplication *app, gpointer data)
{
    GtkBuilder *builder = gtk_builder_new_from_file(CW_PROJECT_DIR "/res/course_work_gtk.glade");

    GObject *window = gtk_builder_get_object(builder,"window");
    g_signal_connect(window,"destroy",G_CALLBACK(gtk_main_quit),NULL);

    gtk_window_set_application(GTK_WINDOW(window),GTK_APPLICATION(app));
    gtk_widget_show(GTK_WIDGET(window));

    init_widgets_l2(builder);
    init_widgets_l3(builder);
}


int main(int argc, char *argv[])
{
    pthread_mutex_init(&l2_context.params_lock,0);
    pthread_cond_init(&l2_context.params_cond,NULL);
    pthread_mutex_init(&l2_context.result_lock,0);
    l2_context.adiabatic = false;
    l2_context.params_dirty = false;

    pthread_mutex_init(&l3_context.params_lock,0);
    pthread_cond_init(&l3_context.params_cond,NULL);
    pthread_mutex_init(&l3_context.result_lock,0);
    l3_context.adiabatic = false;
    l3_context.params_dirty = false;

    GtkApplication *app = gtk_application_new("org.cw.ui",G_APPLICATION_DEFAULT_FLAGS);

    g_signal_connect(app,"activate",G_CALLBACK(activate),NULL);
    int status = g_application_run(G_APPLICATION(app),argc,argv);


    pthread_cancel(l2_context.drawing_thread);
    pthread_cond_destroy(&l2_context.params_cond);
    pthread_mutex_destroy(&l2_context.params_lock);
    pthread_mutex_destroy(&l2_context.result_lock);

    pthread_cancel(l3_context.drawing_thread);
    pthread_cond_destroy(&l3_context.params_cond);
    pthread_mutex_destroy(&l3_context.params_lock);
    pthread_mutex_destroy(&l3_context.result_lock);

    g_object_unref(app);

    return status;
}