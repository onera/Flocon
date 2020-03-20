#
import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse, Polygon
#
from ..simulators import TimeDomainSimulator
from ..simulators import simulators as sim


CLASS_DIR       = os.path.dirname(__file__)
WORKING_DIR     = os.path.join(CLASS_DIR, '.work')
# Meshes associated with the geometry
MESH_DIR        = os.path.join(CLASS_DIR, 'meshes')
TEST_MESH       = os.path.join(MESH_DIR, 'test.msh')
COARSE_MESH     = os.path.join(MESH_DIR, 'coarse.msh')
FINE_MESH       = os.path.join(MESH_DIR, 'fine.msh')
# Templates for the .edp files
TEMPLATE_DIR        = os.path.join(CLASS_DIR, 'edp_templates')
TEMPLATE_BASEFLOW   = os.path.join(TEMPLATE_DIR, 'baseflow')
TEMPLATE_OPENLOOP   = os.path.join(TEMPLATE_DIR, 'openloop')
BASEFLOW_EDP        = os.path.join(WORKING_DIR,'baseflow.edp')
BASEFLOW_CONV       = os.path.join(WORKING_DIR, 'baseflow_conv')
BASEFLOW_DATA       = os.path.join(WORKING_DIR, 'baseflow_data')
BASEFLOW_OTHER_DATA = os.path.join(WORKING_DIR, 'baseflow_other_data')

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
        # Initialise database tables
        self.create_table('baseflow',['Re real', 'dt real', 'N real', 'mesh text', 'data text'])
        # Create working temp directory if it does not exists
        if not os.path.exists(WORKING_DIR):
            os.mkdir(WORKING_DIR)

    def init_default_parameters(self):
        self.name       = 'default'
        # Physical parameters ans setting
        self.Re         = 100.0
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
        self.mesh       = 'test'   # mesh to be used
        self.N          = 250000     # number of iterations
        self.dt         = 0.002      # integration time step
        self.NL         = False      # activation of non-linear term

    def add_to_list(self, li, e):
        if e.outside(xmin, xmax, ymin, ymax):
            print('Element outside the box of the step, x=[%.0f,%.0f] y=[%.0f,%.0f]'%(xmin, xmax, ymin,ymax))
            return
        li.append(e)

    @property
    def tf(self):
        return self.N * self.dt

    @tf.setter
    def tf(self, t):
        self.N = int(round(t/self.dt,0))
        print('Number of iterations changed to {}'.format(self.N))

    @property
    def mesh_file(self):
        if self.mesh=='coarse':
            return COARSE_MESH
        elif self.mesh=='fine':
            return FINE_MESH
        elif self.mesh=='test':
            return TEST_MESH


    @property
    def reduced_id(self):
        return '[%.1f,%d,%.3f,%s]'%(self.Re, self.N, self.dt, self.mesh)
    # --------------------------------------------------------------------------
    # BASE FLOW
    # --------------------------------------------------------------------------
    def compute_baseflow(self, force=False, store=True):
        # Attempt to read base-flow data
        data = self.get_baseflow_data()
        if data and not force:
            self.print_msg('Base-flow data already in DB. Returning...')
            return data
        self.print_msg('Base-flow data not in DB. Starting computation...')
        # Otherwise, data must be computed
        target_Re = self.Re
        # Try to find the closest possible baseflow with existing data to
        # initialise the computation
        self.print_msg('Looking for closest existing base-flow...')
        other = self.get_closest_other()
        if other:
            self.print_msg('Closest base-flow found. Restarting from ' + other.reduced_id)
            # if the data exists, write it to the baseflow_data file
            other_data = other.get_baseflow_data()
            sim.write_file(BASEFLOW_OTHER_DATA, other_data)
            is_restart = True
        else:
            is_restart = False
            self.print_msg('No existing base-flow found. No restart...')
        # Create the EDP file for baseflow computation
        self.make_baseflow_edp_file(is_restart)
        # Launch simulation
        self.print_msg('Launching simulation for base-flow computation...')
        self.clean_temp_files()
        sim.launch_edp_file(BASEFLOW_EDP)
        if self.baseflow_has_converged():
            self.print_msg('Base-flow has converged...')
            data = sim.file_to_str(BASEFLOW_DATA)
            if store:
                self.print_msg('Storing base-flow...')
                self.store_baseflow_data(data)
            return data
        # Simulation has not converged
        self.print_msg('Base-flow has not converged...')
        if other:
            # if it was a restart, an intermediate Re is considered
            target_Re       = self.Re
            self.Re         =  0.5*(self.Re + other.Re)
            self.print_msg('Considering intermediate Re=%.1f' %(self.Re))
            self.compute_baseflow()
            self.Re         = target_Re
            return self.compute_baseflow()
        else:
            # it was not a restart, decrease the Re
            target_Re   = self.Re
            self.Re     = self.Re/2
            self.print_msg('Decreasing Re to %.1f' %(self.Re))
            # And attempt to compute baseflow
            self.compute_baseflow()
            # then retry for target_Re
            self.Re     = target_Re
            return self.compute_baseflow()

    def clean_temp_files(self):
        if os.path.exists(BASEFLOW_CONV):
            os.system('rm '+ BASEFLOW_CONV)
        if os.path.exists(BASEFLOW_DATA):
            os.system('rm '+ BASEFLOW_DATA)

    def get_closest_other(self):
        c       = self.db.cursor()
        diff    = np.inf
        c.execute('SELECT * FROM baseflow WHERE dt=? AND N=? AND mesh=?',
                  (self.dt, self.N, self.mesh))
        matches     = c.fetchall()
        if not matches:
            return []
        dist        = []
        for i,m in enumerate(matches):
            dist.append(np.abs(self.Re - m[0]))
        imin = np.argmin(dist)
        return db_data_to_object(matches[imin])

    def baseflow_has_converged(self):
        conv_flag =  np.loadtxt(BASEFLOW_CONV)
        return conv_flag == 1

    def get_baseflow_data(self):
        data    = []
        c       = self.db.cursor()
        c.execute('SELECT data FROM baseflow WHERE Re=? AND dt=? AND N=? AND mesh=?',
                  (self.Re, self.dt, self.N, self.mesh))
        data = c.fetchone()
        if data:
            data = data[0]
        return data

    def store_baseflow_data(self, data_str):
        c = self.db.cursor()
        c.execute('INSERT INTO baseflow VALUES (?,?,?,?,?)',
                  (self.Re, self.dt, self.N, self.mesh, data_str))
        self.db.commit()

    def get_placeholders(self):
        ph = {'MESH':self.mesh_file,
              'BASEFLOW_CONV': BASEFLOW_CONV,
              'BASEFLOW_DATA': BASEFLOW_DATA,
              'BASEFLOW_OTHER_DATA':BASEFLOW_OTHER_DATA}
        return ph

    def make_baseflow_edp_file(self, is_restart = False):
        # Read associated EDP template
        baseflow_temp   = sim.read_template(TEMPLATE_BASEFLOW)
        #
        content = '// Configuration parameters declaration, case \'%s\'\n' %(self.name)
        content = content + sim.assign_freefem_var('Re', self.Re) + '\n'
        content = content + sim.assign_freefem_var('isRestart', is_restart)
        content = content + '\n// End of parameters declaration \n' + baseflow_temp
        content = sim.replace_placeholders(self.get_placeholders(), content)
        #
        sim.write_file(BASEFLOW_EDP, content)
        return content

    def print_msg(self, msg):
        print('[FLOCON]'+self.reduced_id + ' '+ msg)
    # --------------------------------------------------------------------------
    # PLOTING
    # --------------------------------------------------------------------------
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

def db_data_to_object(data):
    s = Step()
    s.Re    = data[0]
    s.dt    = data[1]
    s.N     = data[2]
    s.mesh  = data[3]
    return s
