#
import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse, Polygon
#
from ..simulators import TimeDomainSimulator

CLASS_DIR       = os.path.dirname(__file__)
# Meshes associated with the geometry
MESH_DIR        = os.path.join(CLASS_DIR, 'meshes')
COARSE_MESH     = os.path.join(MESH_DIR, 'coarse_mesh.msh')
FINE_MESH       = os.path.join(MESH_DIR, 'fine_mesh.msh')
# Templates for the .edp files
TEMPLATE_DIR    = os.path.join(CLASS_DIR, 'edp_template')
EDP_HEAD        = os.path.join(TEMPLATE_DIR, 'head')
EDP_BASEFLOW    = os.path.join(TEMPLATE_DIR, 'baseflow')
EDP_OPENLOOP    = os.path.join(TEMPLATE_DIR, 'openloop')
# Geometry of the step
ymax            =  1.0
ymin            = -1.0
xmax            =  12.0
xmin            = -2.0

class Step(TimeDomainSimulator):
    def __init__(self):
        self.DATABASE_FILE = os.path.join(CLASS_DIR,'STEP.db')
        # Initialise database and default parameters for the benchmark
        self.main_init()

    def init_default_parameters(self):
        self.name       = 'default'
        # Physical parameters ans setting
        self.Re         = 100
        # Control input(s)
        self.actuators  = []
        act             = GaussianIO(-0.05, 0.01, 0.01, 0.1)
        self.add_to_list(self.actuators, act)
        # Noise input(s)
        self.noise      = []
        noi             = GaussianNoiseIO(-0.5, 0.25, 0.1, 0.01)
        self.add_to_list(self.noise, noi)
        # Measured output(s)
        self.sensors    = []
        up_sens         = BorderIntegralIO([-0.35, -0.25], 'bottom')
        self.add_to_list(self.sensors, up_sens)
        # Performance output(s)
        self.perf       = []
        down_sens       = BorderIntegralIO([10.3, 10.7],'bottom')
        self.add_to_list(self.perf, down_sens)
        # Simulation parameter
        self.mesh       = COARSE_MESH   # mesh to be used
        self.N          = 250000        # number of iterations
        self.dt         = 0.002         # integration time step
        self.NL         = False         # activation of non-linear term

    def add_to_list(self, li, e):
        if e.outside(xmin, xmax, ymin, ymax):
            print('Element outside the box of the step, x=[%.0f,%.0f] y=[%.0f,%.0f]'%(xmin, xmax, ymin,ymax))
            return
        li.append(e)

    def plot_config(self):
        [fig, ax]   = plt.subplots()
        # Numerical setup
        txt = '\n'.join((\
        'Re : %.2f' %(self.Re, ),\
        'dt  : %.3f' %(self.dt, ),\
        'N   : %d' %(self.N,),\
        'tf   : %.2f' %(self.tf,)\
        ))
        props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
        ax.text(3, 5, txt, fontsize=14,
                verticalalignment='top', bbox=props)
        # The box
        plt.plot(np.array([xmin, xmin, xmax, xmax, xmin]),\
                 np.array([ymin, ymax, ymax, ymin, ymin]), color='black')
        # The step
        plt.fill(np.array([xmin, xmin, 0, 0, xmax,xmin]),\
                 np.array([ymin, 0, 0, ymin, ymin, ymin]), color='grey')
        # The noise(s)
        n_x = xmin -1
        n_y = ymax + 1
        plot_io(plt, ax, 'left', self.noise, 'w', n_x, n_y, 0.5, 'red')
        # The actuator(s)
        u_x = xmin - 1
        u_y = ymin - 1
        plot_io(plt, ax, 'left', self.actuators, 'u', u_x, u_y, -0.5, 'blue')
        # The performane(s)
        p_x = xmax + 1
        p_y = ymax + 1
        plot_io(plt, ax, 'right', self.perf, 'z', p_x, p_y, 0.5, 'orange')
        # The sensor(s)
        s_x = xmax + 1
        s_y = ymin - 1
        plot_io(plt, ax, 'right', self.sensors, 'y', s_x, s_y, -0.5, 'green')
        #
        ax.axis('equal')
        ax.set_xlim(xmin-2,xmax+2)
        ax.set_ylim(-4,5)
        plt.title('Step use-case \"%s\"' %(self.name))
        plt.show()

    @property
    def tf(self):
        return self.N * self.dt

    @tf.setter
    def tf(self, t):
        self.N = int(round(t/self.dt,0))
        print('Number of iterations changed to {}'.format(self.N))

def plot_io(plt, ax, side, li, id, x, y , dy, col):
    if side =='left':
        si = 1
    elif side =='right':
        si = -1
    for i,e in enumerate(li):
        # Name of the input
        txt = '$%s_%d$'%(id, i+1,)
        ax.text(x - si * 0.3, y, txt, transform=ax.transData,fontsize=9,
                horizontalalignment=side, verticalalignment='center',color=col)
        # Line to show where it acts
        plt.plot(np.array([x + si* 0.5, e.xp, e.xp]),\
                 np.array([y, y, e.yp]), color=col, linestyle=':')
        # Action
        e.plot(plt, ax, col)
        # Increment for next input
        y = y + dy

class SystemIO:
    def __init__(self, x, y):
        self.x      = x
        self.y      = y
    def outside(self, xmin, xmax, ymin, ymax):
        if self.x < xmin or self.x > xmax or self.y < ymin or self.y > ymax:
            return True
        return False
    @property
    def xp(self):
        return self.x
    @property
    def yp(self):
        return self.y

# Integral action
class BorderIntegralIO(SystemIO):
    def __init__(self, x, side):
        super().__init__(x, [])
        self.side = side

    def outside(self, xmin, xmax, ymin, ymax):
        if self.x[0] < xmin or self.x[1] > xmax:
            return True
        return False
    @property
    def xp(self):
        return np.mean(self.x)

    @property
    def yp(self):
        if self.side=='top':
            return ymax
        else:
            if self.x[1]<0:
                return 0
            else:
                return ymin

    def plot(self, plt, ax, color):
        plt.plot(self.x, self.yp*np.ones(2), color=color, linewidth=3)

# Gaussian action
class GaussianIO(SystemIO):
    def __init__(self, x, y, sigma_x, sigma_y = []):
        super().__init__(x, y)
        self.sigma_x = sigma_x
        if sigma_y:
            self.sigma_y = sigma_y
        else:
            self.sigma_y = sigma_x

    def plot(self, plt, ax, color):
        sigs    = [1, 2, 3]
        al      = 1
        for s in sigs:
            ell = Ellipse((self.x, self.y), s*self.sigma_x, s*self.sigma_y, color=color,zorder=0)
            ax.add_artist(ell)
            ell.set_alpha(al)
            al = al/2

class GaussianNoiseIO(GaussianIO):
    def __init__(self, x, y, sigma_x,  sigma_time, sigma_y = []):
        super().__init__(x, y, sigma_x, sigma_y)
        self.sigma_time = sigma_time
