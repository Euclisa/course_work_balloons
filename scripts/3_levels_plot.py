import matplotlib.pyplot as plt
import numpy as np
from math import pi


def plot_data(fn: str, linestyle='solid', linecolor='black'):
    with open(fn,'r') as f:
        data_raw = f.read().split('\n')
        if data_raw[-1] == '':
            data_raw = data_raw[:-1]
        data = list(map(lambda x: float(x), data_raw))

    phi_ad = data[0]
    r_ad = data[1]
    x_ad = data[2]
    y_ad = data[3]
    a_ad = pi-data[4]
    phi_cb = data[5]
    r_cb = data[6]
    x_cb = data[7]
    y_cb = data[8]
    a_cb = pi-data[9]
    phi_dc = data[10]
    r_dc = data[11]
    x_dc = data[12]
    y_dc = data[13]
    a_dc = pi-data[14]
    phi_df = data[15]
    r_df = data[16]
    x_df = data[17]
    y_df = data[18]
    a_df = pi-data[19]
    phi_ec = data[20]
    r_ec = data[21]
    x_ec = data[22]
    y_ec = data[23]
    a_ec = pi-data[24]
    phi_fe = data[25]
    r_fe = data[26]
    x_fe = data[27]
    y_fe = data[28]
    a_fe = pi-data[29]
    phi_ge = data[30]
    r_ge = data[31]
    y_ge = data[32]
    a_ge = -pi/2
    phi_gf = data[33]
    r_gf = data[34]
    y_gf = data[35]
    a_gf = -pi/2
    x_bot = data[36]

    grain = 100

    plt.axis('equal')

    t_ad = np.linspace(a_ad-phi_ad,a_ad,grain)
    xs_ad = x_ad + r_ad*np.cos(t_ad)
    ys_ad = y_ad + r_ad*np.sin(t_ad)
    plt.plot(xs_ad,ys_ad,linestyle=linestyle,color=linecolor)

    t_cb = np.linspace(a_cb-phi_cb,a_cb,grain)
    xs_cb = x_cb + r_cb*np.cos(t_cb)
    ys_cb = y_cb + r_cb*np.sin(t_cb)
    plt.plot(xs_cb,ys_cb,linestyle=linestyle,color=linecolor)

    t_dc = np.linspace(a_dc-phi_dc,a_dc,grain)
    xs_dc = x_dc + r_dc*np.cos(t_dc)
    ys_dc = y_dc + r_dc*np.sin(t_dc)
    plt.plot(xs_dc,ys_dc,linestyle=linestyle,color=linecolor)

    t_ec = np.linspace(a_ec-phi_ec,a_ec,grain)
    xs_ec = x_ec + r_ec*np.cos(t_ec)
    ys_ec = y_ec + r_ec*np.sin(t_ec)
    plt.plot(xs_ec,ys_ec,linestyle=linestyle,color=linecolor)

    t_fe = np.linspace(a_fe-phi_fe,a_fe,grain)
    xs_fe = x_fe + r_fe*np.cos(t_fe)
    ys_fe = y_fe + r_fe*np.sin(t_fe)
    plt.plot(xs_fe,ys_fe,linestyle=linestyle,color=linecolor)

    t_df = np.linspace(a_df-phi_df,a_df,grain)
    xs_df = x_df + r_df*np.cos(t_df)
    ys_df = y_df + r_df*np.sin(t_df)
    plt.plot(xs_df,ys_df,linestyle=linestyle,color=linecolor)

    t_ge = np.linspace(a_ge-phi_ge,a_ge,grain)
    xs_ge = x_bot + r_ge*np.cos(t_ge)
    ys_ge = y_ge + r_ge*np.sin(t_ge)
    plt.plot(xs_ge,ys_ge,linestyle=linestyle,color=linecolor)

    t_gf = np.linspace(a_gf+phi_gf,a_gf,grain)
    xs_gf = x_bot + r_gf*np.cos(t_gf)
    ys_gf = y_gf + r_gf*np.sin(t_gf)
    plt.plot(xs_gf,ys_gf,linestyle=linestyle,color=linecolor)


plot_data('3_levels.txt')
plot_data('3_levels_init.txt',linestyle='dotted')
plt.show()