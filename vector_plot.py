from matplotlib.pyplot import cm
import warnings
warnings.filterwarnings("ignore")
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import fsolve
import sympy as sm
from matplotlib import style as style
from matplotlib.backends.backend_pdf import PdfPages
import seaborn.apionly as sns
plt.style.use('ggplot')
plt.rc('text', usetex=True)
plt.rc('font', family='sans-serif')
plt.rcParams['font.size'] = 10
plt.rcParams['axes.labelsize'] = 10
#plt.rcParams['axes.labelweight'] = 'bold'
plt.rcParams['xtick.labelsize'] = 8
plt.rcParams['ytick.labelsize'] = 8
plt.rcParams['legend.fontsize'] = 10
plt.rcParams['figure.titlesize'] = 12

#https://www.getdatajoy.com/examples/python-plots/vector-fields


def find_roots():
    def equations(p):
        return (((1 / (1 + (y1 / 150))) + 0.1 * (1 - (1 / (1 + (y1 / 150))))) * ((1 / (1 + (x1 / 100) ** 3)) +10 * (1 - (1 / (1 + (x1 / 100) ** 3)))) * 4 - 0.1 * x1,
                ((1 / (1 + (x1 / 150))) + 0.1 * (1 - (1 / (1 + (x1 / 150))))) * ((1 / (1 + (y1 / 100) ** 3)) + 10 * (1 - (1 / (1 + (y1 / 100) ** 3)))) * 4 - 0.1 * y1)
        #return (((1 / (1 + (y1 / 200) ** 3)) + 0.1 * (1 - (1 / (1 + (y1 / 200) ** 3)))) * 40 - 0.1 * x1,
        #        ((1 / (1 + (x1 / 200) ** 3)) + 0.1 * (1 - (1 / (1 + (x1 / 200) ** 3)))) * 40 - 0.1 * y1)

    x1, y1 = fsolve(equations, (0, 0))
    x2, y2 = fsolve(equations, (10, 0))
    x3, y3 = fsolve(equations, (176, 176))
    x4, y4 = fsolve(equations, (235, 105))
    x5, y5 = fsolve(equations, (360, 15))
    roots_x_st = [x1, x3, x5]
    roots_y_st = [y1, y3, y5]

    roots_x_ust = [x2, x4]
    roots_y_ust = [y2, y4]

    return roots_x_st, roots_y_st, roots_x_ust, roots_y_ust

def plot_strm(x, y, UN, VN, roots_x_st, roots_y_st, roots_x_ust, roots_y_ust, speed):
    pp = PdfPages("stream_plot_DP.pdf")
    plot3 = plt.figure()
    fft_axes = plot3.add_subplot(111)
    fft_axes.set_xlim([0,10])
    fft_axes.set_ylim([0,10])
    plt.streamplot(x, y, UN, VN,          # data
               color=speed,         # array that determines the colour
               cmap=cm.autumn,        # colour map
               linewidth=1.5,         # line thickness
               arrowstyle='->',     # arrow style
               arrowsize=1)       # arrow size

    plt.colorbar(shrink=.5)
    plt.scatter(roots_x_st, roots_y_st, color='#27408B', s=50)
    plt.scatter(roots_x_ust, roots_y_ust, color='grey',s=50)
    plt.xlabel('XX')
    plt.ylabel('YY')
    pp.savefig()
    plt.close()
    pp.close()
    return

y, x = np.mgrid[0:10:18j, 0:10:18j]
#DP
U = ((1 / (1 + (y / 150))) + 0.1 * (1 - (1 / (1 + (y / 150))))) * ((1 / (1 + (x / 100) ** 3)) +10 * (1 - (1 / (1 + (x / 100) ** 3)))) * 4 - 0.1 * x
V = ((1 / (1 + (x / 150))) + 0.1 * (1 - (1 / (1 + (x / 150))))) * ((1 / (1 + (y / 100) ** 3)) + 10 * (1 - (1 / (1 + (y / 100) ** 3)))) * 4 - 0.1 * y
#CS
U_CS = ((1 / (1 + (y / 200) ** 3)) +0.1 * (1 - (1 / (1 + (y / 200) ** 3)))) * 40 -0.1 * x
V_CS = ((1 / (1 + (x / 200) ** 3)) + 0.1 * (1 - (1 / (1 + (x / 200) ** 3)))) * 40 -0.1 * y


speed = np.sqrt(U**2 + V**2)
UN = U/speed
VN = V/speed

roots_x_st, roots_y_st, roots_x_ust, roots_y_ust = find_roots()
plot_strm(x, y, UN, VN, roots_x_st, roots_y_st, roots_x_ust, roots_y_ust, speed)