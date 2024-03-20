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
    phi_ed = data[15]
    r_ed = data[16]
    y_ed = data[17]
    a_ed = 3*pi/2
    phi_ec = data[18]
    r_ec = data[19]
    y_ec = data[20]
    a_ec = 3*pi/2
    x_bot = data[21]


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

    t_ed = np.linspace(a_ed,a_ed+phi_ed,grain)
    xs_ed = x_bot + r_ed*np.cos(t_ed)
    ys_ed = y_ed + r_ed*np.sin(t_ed)
    plt.plot(xs_ed,ys_ed,linestyle=linestyle,color=linecolor)

    t_ec = np.linspace(a_ec-phi_ec,a_ec,grain)
    xs_ec = x_bot + r_ec*np.cos(t_ec)
    ys_ec = y_ec + r_ec*np.sin(t_ec)
    plt.plot(xs_ec,ys_ec,linestyle=linestyle,color=linecolor)


plot_data('2_levels.txt')
plot_data('2_levels_init.txt',linestyle='dotted')
plt.show()