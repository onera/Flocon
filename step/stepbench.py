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
SIMOUT              = os.path.join(WORKING_DIR, 'simout')
SIMIN               = os.path.join(WORKING_DIR, 'simin')
# Geometry of the step
ymax            =  1.0
ymin            = -1.0
xmax            =  12.0
xmin            = -2.0


# TODO :
#  - check for longer experiment or comensurable ones for baseflow data
#
class Step(TimeDomainSimulator):
    def __init__(self):
        self.DATABASE_FILE = os.path.join(CLASS_DIR,'STEP.db')
        # Initialise database and default parameters for the benchmark
        self.main_init()
        # Initialise database tables
        self.create_table('baseflow',['Re real', 'dt real', 'N real', 'mesh text', 'data text'])
        self.create_table('openloop',['config text', 'ioconfig text', \
                                      'casename text', 'inputs array', 'outputs array'])
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
        self.noises     = []
        noi             = GaussianNoiseIO(-0.5, 0.25, 0.1, 0.01)
        self.add_to_list(self.noises, noi)
        # Measured output(s)
        self.sensors    = []
        up_sens         = BorderIntegralIO([-0.35, -0.25], 'bottom')
        self.add_to_list(self.sensors, up_sens)
        # Performance output(s)
        self.performances = []
        down_sens       = BorderIntegralIO([10.3, 10.7],'bottom')
        self.add_to_list(self.performances , down_sens)
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
        return '%.1f,%d,%.3f,%s,%r'%(self.Re, self.N, self.dt, self.mesh, self.NL)

    @property
    def nu(self):
        return len(self.actuators)
    @property
    def nw(self):
        return len(self.noises)
    @property
    def ny(self):
        return len(self.sensors)
    @property
    def nz(self):
        return len(self.performances)

    @property
    def ioconfig_str(self):
        out = ' U:'
        for a in self.actuators:
            out = out + a.str_id
        out = out + ' W:'
        for n in self.noises:
            out = out + n.str_id
        out = out + ' Y:'
        for s in self.sensors:
            out = out + s.str_id
        out = out + ' Z:'
        for p in self.performances :
            out = out + p.str_id
        return out

    def print_msg(self, msg):
        print('[FLOCON]['+self.reduced_id + '] '+ msg)
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
        if os.path.exists(SIMOUT):
            os.system('rm '+ SIMOUT)

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
              'BASEFLOW_OTHER_DATA':BASEFLOW_OTHER_DATA,
              'SIMOUT':SIMOUT,
              'SIMIN':SIMIN}
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
    # --------------------------------------------------------------------------
    # OPEN LOOP
    # --------------------------------------------------------------------------
    def simulate_openloop(self, name, signals, force=False):
        # Check whether this simulation already exists based on its name
        data = selt.get_openloop_simulation(name)
        if not force and data:
            self.print_msg('Open-loop simulation named %s already found, returning.' %(name))
            return data
        # Try to get baseflow data
        bf_data = self.get_baseflow_data()
        if not bf_data:
            self.print_msg('No existing baseflow for this configuration...Compute it first.')
            return
        # Check dimensions of input signal
        err_msg = self.check_input_signals(signals)
        if err_msg:
            self.print_msg(err_msg)
            return
        # Store the input signal into the file SIMIN
        sim.np_to_freefem_file(SIMIN, signals)
        # At this point, the simulation can be launched
        # Store baseflow data in associated file
        sim.write_file(BASEFLOW_OTHER_DATA, bf_data)
        # Assemble EDP file for open-loop simulation
        self.make_openloop_edp_file()

    def make_openloop_edp_file(self):
        # Read associated EDP template
        openloop_temp   = sim.read_template(TEMPLATE_OPENLOOP)
        #
        content = '// Configuration parameters declaration, case\n'
        content = content + sim.assign_freefem_var('Re', self.Re) + '\n'
        content = content + sim.assign_freefem_var('NL', self.NL) + '\n'

        content = content + '\n// End of parameters declaration \n' + openloop_temp
        content = sim.replace_placeholders(self.get_placeholders(), content)
        #
        sim.write_file(OPENLOOP_EDP, content)
        return content

    def get_openloop_simulation(self, name):
        data    = []
        c       = self.db.cursor()
        c.execute('SELECT inputs AND outputs FROM openloop WHERE config=? AND ioconfig=? AND casename=?',
                  (self.reduced_id, self.ioconfig_str, name))
        data = c.fetchone()
        if data:
            data = data[0]
        return data

    def check_input_signals(self, signals):
        # Number of input signals
        n = size(signals,2)
        if n>self.nu:
            return 'Too many input signals...%d signals,  excepted %d'%(n, self.nu)
        elif n<self.nu:
            return 'Not enough input signals...%d signals,  excepted %d'%(n, self.nu)
        # signal length
        nt = size(signals,1)
        if nt>self.N:
            return 'Signal has too many elements...%d, expected %d'%(nt, self.N)
        elif nt < self.N:
            return 'Signal has not enough elements...%d, expected %d'%(nt, self.N)
        # Otherwise, signal has the right dimensions
        return ''

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
        plot_io(plt, ax, 'left', self.noises, 'w', n_x, n_y, 0.5, 'red')
        # The actuator(s)
        u_x = xmin - 1
        u_y = ymin - 1
        plot_io(plt, ax, 'left', self.actuators, 'u', u_x, u_y, -0.5, 'blue')
        # The performane(s)
        p_x = xmax + 1
        p_y = ymax + 1
        plot_io(plt, ax, 'right', self.performances , 'z', p_x, p_y, 0.5, 'orange')
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

    @property
    def str_id(self):
        return 'BI,%f,%f,%s'%(self.x[0],self.x[1],self.side)

    def to_fem(self, cat, id):
        if self.side=='top':
            side = '(y>=0.1)'
        else:
            side = '(y<=0.1)'
        output_id   = '%sBI%d'%(cat,id)
        decl        = 'Up support_%s=(x>=%.2f)*(x<=%2f)*%s'%(output_id, self.x[0], self.x[1], side)
        init        = 'real %s = int1d(th,2)(dy(u)*support_%s)'%(output_id, output_id)
        assign      = '%s = int1d(th,2)(dy(u)*support_%s)'%(output_id, output_id)
        return [output_id, decl, init, assign]

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
            al  = al/2
    @property
    def str_id(self):
        return 'G,%f,%f,%f,%f'%(self.x, self.y, self.sigma_x, self.sigma_y)

    def to_fem(self, cat, id):
        # ID
        input_id    = '%sG%d'%(cat,id)
        # Declarations
        # Actuator
        x_decl      = sim.assign_freefem_var('x_%s'%(input_id), self.x)
        sx_decl     = sim.assign_freefem_var('sx_%s'%(input_id), self.sigma_x)
        y_decl      = sim.assign_freefem_var('y_%s'%(input_id), self.y)
        sy_decl     = sim.assign_freefem_var('sy_%s'%(input_id), self.sigma_y)
        decl        = '\n'.join([x_decl + sx_decl +y_decl + sy_decl +\
                       'Uvvp [%s_1, %s_2, %s_3];'%(input_id,input_id,input_id),\
                       '{',\
                       'func real in_act_%s(real x, real y)'%(input_id),\
                       'return exp(-(x-x_%s).^2/(2*(sx_%s).^2))*exp(-(y-y_%s).^2/(2*(sy_%s).^2));'%(input_id, input_id, input_id, input_id),\
                       '};',\
                       '[rhs_%s_1, rhs_%s_2, rhs_%s_3]  = [0,in_act_%s(x,y),0]'%(input_id,input_id,input_id,input_id),\
                       '%s_1[]             = MatMass * rhs_%s_1[]'%(input_id,input_id),\
                       '};'])
        # Init: input signal read
        init = ''
        # init = '\n'.join(['// Reading input signal %s from input matrix'%(input_id),\
        #                    'real[int] cmd_%s(N);'%(input_id),\
        #                    'cmd_%s=0.;'%(input_id),\
        #                    '{',\
        #                    'int in_size;',\
        #                    'ifstream file("@SIMIN");',\
        #                    'file >> in_size;',\
        #                    'for (int i=0; i< in_size; i++)',\
        #                    '{',\
        #                    'file >> cmd_%s(i)'%(input_id),\
        #                    '};',\
        #                    '};'])
        # Associated
        assign = 'rhs1[] += input_signals[i,%d]*%s_1[]; // input %s' %(id, input_id, input_id)
        return [input_id, decl, init, assign]


class GaussianNoiseIO(GaussianIO):
    def __init__(self, x, y, sigma_x,  sigma_time, sigma_y = []):
        super().__init__(x, y, sigma_x, sigma_y)
        self.sigma_time = sigma_time
    @property
    def str_id(self):
        str = super().str_id
        str = str.replace('G','GN')
        return str +',%f'%(self.sigma_time)

def db_data_to_object(data):
    s       = Step()
    s.Re    = data[0]
    s.dt    = data[1]
    s.N     = data[2]
    s.mesh  = data[3]
    return s
