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

pthread_mutex_t params_lock;
pthread_cond_t params_cond;
struct system_2_levels_user_params user_params;
bool params_dirty = false;
bool adiabatic = false;
pthread_t drawing_thread;
pthread_mutex_t result_lock;
struct system_2_levels_result result;
struct adiabatic_mode_widgets adia_widgets;


static void draw_function(GtkDrawingArea *area, cairo_t *cr, gpointer data)
{
    pthread_mutex_lock(&result_lock);
    
    const int width = gtk_widget_get_allocated_width(GTK_WIDGET(area));
    const int height = gtk_widget_get_allocated_height(GTK_WIDGET(area));

    const int pixels_per_meter = MIN(width,height)/2;

    double r_ad = result.r_ad*pixels_per_meter;
    double r_cb = result.r_cb*pixels_per_meter;
    double r_dc = result.r_dc*pixels_per_meter;
    double r_ec = result.r_ec*pixels_per_meter;
    double r_ed = result.r_ed*pixels_per_meter;
    double x_ad = result.x_ad*pixels_per_meter;
    double y_ad = result.y_ad*pixels_per_meter;
    double x_cb = result.x_cb*pixels_per_meter;
    double y_cb = result.y_cb*pixels_per_meter;
    double x_dc = result.x_dc*pixels_per_meter;
    double y_dc = result.y_dc*pixels_per_meter;
    double y_ec = result.y_ec*pixels_per_meter;
    double y_ed = result.y_ed*pixels_per_meter;
    double x_bot = result.x_bot*pixels_per_meter;

    cairo_translate(cr,0,height);
    cairo_scale(cr,1,-1);

    double a_ad = G_PI-result.a_ad;
    cairo_arc(cr,x_ad,y_ad,r_ad,a_ad-result.phi_ad,a_ad);
    cairo_stroke(cr);

    double a_dc = G_PI-result.a_dc;
    cairo_arc(cr,x_dc,y_dc,r_dc,a_dc-result.phi_dc,a_dc);
    cairo_stroke(cr);
    
    double a_cb = G_PI-result.a_cb;
    cairo_arc(cr,x_cb,y_cb,r_cb,a_cb-result.phi_cb,a_cb);
    cairo_stroke(cr);
    
    double a_ec = -G_PI_2;
    cairo_arc(cr,x_bot,y_ec,r_ec,a_ec-result.phi_ec,a_ec);
    cairo_stroke(cr);

    double a_ed = -G_PI_2;
    cairo_arc(cr,x_bot,y_ed,r_ed,a_ed,a_ed+result.phi_ed);
    cairo_stroke(cr);
    
    
    cairo_select_font_face(cr,"monospace",CAIRO_FONT_SLANT_NORMAL,CAIRO_FONT_WEIGHT_NORMAL);
    cairo_set_font_size(cr,20);
    cairo_scale(cr,1,-1);

    double letter_offset = height / 20;

    double x_a = x_ad + r_ad*cos(a_ad) - letter_offset;
    double y_a = -y_ad - r_ad*sin(a_ad);
    cairo_move_to(cr,x_a,y_a);
    cairo_show_text(cr,"A");

    double x_b = x_cb + r_cb*cos(a_cb-result.phi_cb) - letter_offset;
    double y_b = -y_cb - r_cb*sin(a_cb-result.phi_cb);
    cairo_move_to(cr,x_b,y_b);
    cairo_show_text(cr,"B");

    double x_c = x_cb + r_cb*cos(a_cb) - letter_offset;
    double y_c = -y_cb - r_cb*sin(a_cb) + letter_offset/2;
    cairo_move_to(cr,x_c,y_c);
    cairo_show_text(cr,"C");

    double x_d = x_ad + r_ad*cos(a_ad+result.phi_ad) + letter_offset/4;
    double y_d = -y_ad - r_ad*sin(a_ad+result.phi_ad) + letter_offset;
    cairo_move_to(cr,x_d,y_d);
    cairo_show_text(cr,"D");

    double x_e = x_bot + r_ec*cos(a_ec);
    double y_e = -y_ec - r_ec*sin(a_ec) + letter_offset;
    cairo_move_to(cr,x_e,y_e);
    cairo_show_text(cr,"E");

    pthread_mutex_unlock(&result_lock);
}


static void queue_update_picture(GtkDrawingArea *area, const struct system_2_levels_user_params *params_extracted, bool adiabatic_extracted)
{
    struct system_2_levels_result result_local;

    if(adiabatic_extracted)
        system_2_levels_adiabatic_eval(params_extracted,&result_local);
    else
        system_2_levels_eval(params_extracted,&result_local);


    pthread_mutex_lock(&result_lock);

    memcpy(&result,&result_local,sizeof(struct system_2_levels_result));
    gtk_widget_queue_draw(GTK_WIDGET(area));

    pthread_mutex_unlock(&result_lock);
}


static void *update_picture(GtkDrawingArea *area)
{
    while(true)
    {
        pthread_mutex_lock(&params_lock);
        if(!params_dirty)
            pthread_cond_wait(&params_cond,&params_lock);
        
        struct system_2_levels_user_params params_extracted;
        memcpy(&params_extracted,&user_params,sizeof(struct system_2_levels_user_params));
        bool adiabatic_extracted = adiabatic;
        params_dirty = false;

        pthread_mutex_unlock(&params_lock);

        queue_update_picture(area,&params_extracted,adiabatic_extracted);
    }

    return NULL;
}


static void spin_button_value_changed_cb(GtkSpinButton *spin_button, double *param)
{
    pthread_mutex_lock(&params_lock);
    double value = (double)gtk_spin_button_get_value(spin_button);
    *param = value;
    params_dirty = true;
    pthread_cond_signal(&params_cond);
    pthread_mutex_unlock(&params_lock);
}

static void adiabatic_mode_changed_cb(GtkSwitch *adiabatic_sw, gpointer *data)
{
    pthread_mutex_lock(&params_lock);
    adiabatic = (bool)gtk_switch_get_active(adiabatic_sw);
    gtk_widget_set_sensitive(GTK_WIDGET(adia_widgets.adiabatic_constant_label),adiabatic);
    gtk_widget_set_sensitive(GTK_WIDGET(adia_widgets.adiabatic_constant_spin),adiabatic);

    params_dirty = true;
    pthread_cond_signal(&params_cond);
    pthread_mutex_unlock(&params_lock);
}


static void activate(GtkApplication *app, gpointer data)
{
    GtkBuilder *builder = gtk_builder_new_from_file(CW_PROJECT_DIR "/res/course_work_gtk.glade");

    GObject *window = gtk_builder_get_object(builder,"window");

    gtk_window_set_application(GTK_WINDOW(window),GTK_APPLICATION(app));
    gtk_widget_show(GTK_WIDGET(window));

    GObject *ax_spin = gtk_builder_get_object(builder,"ax");
    g_signal_connect(ax_spin,"value-changed",G_CALLBACK(spin_button_value_changed_cb),&user_params.Ax);
    spin_button_value_changed_cb(GTK_SPIN_BUTTON(ax_spin),&user_params.Ax);

    GObject *ay_spin = gtk_builder_get_object(builder,"ay");
    g_signal_connect(ay_spin,"value-changed",G_CALLBACK(spin_button_value_changed_cb),&user_params.Ay);
    spin_button_value_changed_cb(GTK_SPIN_BUTTON(ay_spin),&user_params.Ay);

    GObject *bx_spin = gtk_builder_get_object(builder,"bx");
    g_signal_connect(bx_spin,"value-changed",G_CALLBACK(spin_button_value_changed_cb),&user_params.Bx);
    spin_button_value_changed_cb(GTK_SPIN_BUTTON(bx_spin),&user_params.Bx);

    GObject *by_spin = gtk_builder_get_object(builder,"by");
    g_signal_connect(by_spin,"value-changed",G_CALLBACK(spin_button_value_changed_cb),&user_params.By);
    spin_button_value_changed_cb(GTK_SPIN_BUTTON(by_spin),&user_params.By);

    GObject *phi_ad_spin = gtk_builder_get_object(builder,"phi_ad");
    g_signal_connect(phi_ad_spin,"value-changed",G_CALLBACK(spin_button_value_changed_cb),&user_params.phi_ad_0);
    spin_button_value_changed_cb(GTK_SPIN_BUTTON(phi_ad_spin),&user_params.phi_ad_0);

    GObject *phi_dc_spin = gtk_builder_get_object(builder,"phi_dc");
    g_signal_connect(phi_dc_spin,"value-changed",G_CALLBACK(spin_button_value_changed_cb),&user_params.phi_dc_0);
    spin_button_value_changed_cb(GTK_SPIN_BUTTON(phi_dc_spin),&user_params.phi_dc_0);

    GObject *r_top_spin = gtk_builder_get_object(builder,"r_top");
    g_signal_connect(r_top_spin,"value-changed",G_CALLBACK(spin_button_value_changed_cb),&user_params.r_top_0);
    spin_button_value_changed_cb(GTK_SPIN_BUTTON(r_top_spin),&user_params.r_top_0);

    GObject *r_bot_spin = gtk_builder_get_object(builder,"r_bot");
    g_signal_connect(r_bot_spin,"value-changed",G_CALLBACK(spin_button_value_changed_cb),&user_params.r_bot_0);
    spin_button_value_changed_cb(GTK_SPIN_BUTTON(r_bot_spin),&user_params.r_bot_0);

    GObject *p_top_spin = gtk_builder_get_object(builder,"p_top");
    g_signal_connect(p_top_spin,"value-changed",G_CALLBACK(spin_button_value_changed_cb),&user_params.p_top_0);
    spin_button_value_changed_cb(GTK_SPIN_BUTTON(p_top_spin),&user_params.p_top_0);

    GObject *p_bot_spin = gtk_builder_get_object(builder,"p_bot");
    g_signal_connect(p_bot_spin,"value-changed",G_CALLBACK(spin_button_value_changed_cb),&user_params.p_bot_0);
    spin_button_value_changed_cb(GTK_SPIN_BUTTON(p_bot_spin),&user_params.p_bot_0);

    GObject *p_ship_spin = gtk_builder_get_object(builder,"p_ship");
    g_signal_connect(p_ship_spin,"value-changed",G_CALLBACK(spin_button_value_changed_cb),&user_params.p_ac);
    spin_button_value_changed_cb(GTK_SPIN_BUTTON(p_ship_spin),&user_params.p_ac);

    GObject *p_atm_spin = gtk_builder_get_object(builder,"p_atm");
    g_signal_connect(p_atm_spin,"value-changed",G_CALLBACK(spin_button_value_changed_cb),&user_params.p_atm);
    spin_button_value_changed_cb(GTK_SPIN_BUTTON(p_atm_spin),&user_params.p_atm);


    adia_widgets.adiabatic_constant_label = GTK_LABEL(gtk_builder_get_object(builder,"adiabatic_constant_label"));
    adia_widgets.adiabatic_constant_spin = GTK_SPIN_BUTTON(gtk_builder_get_object(builder,"adiabatic_constant_spin"));
    g_signal_connect(adia_widgets.adiabatic_constant_spin,"value-changed",G_CALLBACK(spin_button_value_changed_cb),&user_params.k);
    spin_button_value_changed_cb(GTK_SPIN_BUTTON(adia_widgets.adiabatic_constant_spin),&user_params.k);

    GObject *adiabatic_sw = gtk_builder_get_object(builder,"adiabatic");
    g_signal_connect(adiabatic_sw,"state-set",G_CALLBACK(adiabatic_mode_changed_cb),NULL);
    adiabatic_mode_changed_cb(GTK_SWITCH(adiabatic_sw),NULL);

    GtkDrawingArea *area = GTK_DRAWING_AREA(gtk_builder_get_object(builder,"drawing_area"));
    g_signal_connect(G_OBJECT(area),"draw",G_CALLBACK(draw_function),NULL);

    params_dirty = true;

    pthread_create(&drawing_thread,NULL,(void *(*)(void*))update_picture,area);
}


int main(int argc, char *argv[])
{
    system_2_levels_eval_f();
    pthread_mutex_init(&params_lock,0);
    pthread_cond_init(&params_cond,NULL);
    pthread_mutex_init(&result_lock,0);

    struct system_2_levels_params params;

    GtkApplication *app = gtk_application_new("org.cw.ui",G_APPLICATION_DEFAULT_FLAGS);

    g_signal_connect(app,"activate",G_CALLBACK(activate),&params);
    int status = g_application_run(G_APPLICATION(app),argc,argv);

    pthread_cancel(drawing_thread);
    g_object_unref(app);

    return status;
}