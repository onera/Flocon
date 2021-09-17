#
import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse, Polygon
#
from ..simulators import TimeDomainSimulator
from ..simulators import simulators as sim

CLASS_DIR           = os.path.dirname(__file__)
WORKING_DIR         = os.path.join(CLASS_DIR, 'work')
# Meshes associated with the geometry
MESH_DIR            = os.path.join(CLASS_DIR, 'meshes')
BASE_MESH           = os.path.join(MESH_DIR, 'mesh.msh')
# Templates for the .edp files
TEMPLATE_DIR        = os.path.join(CLASS_DIR, 'edp_templates')
TEMPLATE_EIG        = os.path.join(TEMPLATE_DIR, 'eig')

EIG_EDP             = os.path.join(WORKING_DIR, 'eig.edp')


class Ks(TimeDomainSimulator):
    def __init__(self):
        self.DATABASE_FILE = os.path.join(CLASS_DIR,'KS.db')
        # Initialise database and default parameters for the benchmark
        self.main_init()
        # Initialise database tables
        #self.create_table('baseflow',['Re real', 'mesh text', 'data text'])
        # Create working temp directory if it does not exists
        if not os.path.exists(WORKING_DIR):
            os.mkdir(WORKING_DIR)
        # Some preliminary tests
        # FreeFem++ executable file
        out         = sim.find_ex('FreeFem++')
        if out is None:
            raise Exception("FreeFem++ not found. Make sure it is installed.")
        # Trying to determine whether the MUMPS solver is available in FreeFem++ install
        self.test_freefem_solver()

    def init_default_parameters(self):
        self.name       = 'default'
        # Physical parameters and setting
        self.ub         = 1.0
        self.gamma      = 1.0
        self.d          = 1.0
        self.mu0        = 3.95
        #
        self.xmin       = -100.0
        self.xmax       = 100.0
        self.deltax     = 0.05
        self.lambda0    = complex(0.1,0.58)
        # Simulation parameter
        self.mesh       = 'base'     # mesh to be used
        self.N          = 12500     # number of iterations
        self.dt         = 0.02       # integration time step
        self.NL         = True       # activation of non-linear term
    @property
    def tf(self):
        return self.N * self.dt

    @tf.setter
    def tf(self, t):
        self.N = int(round(t/self.dt,0))
        self.print_msg('Number of iterations updated...')

    @property
    def mesh_file(self):
        if self.mesh=='base':
            return BASE_MESH
        else:
            self.print_msg('Unaivalable value for mesh...')

    @property
    def reduced_id(self):
        return '%.2f,%d,%.3f,%s,%r'%(self.mu0, self.N, self.dt, self.mesh, self.NL)

    def print_msg(self, msg):
        print('[FLOCON]['+self.reduced_id + '] '+ msg)

    def get_placeholders(self):
        ph = {'MESH':self.mesh_file,
              'SOLVER':''}
        if self.solver == 'mumps':
            # ph['SOLVER'] = 'load "MUMPS_seq"\ndefaulttoMUMPSseq();'
            ph['SOLVER'] = 'load "MUMPS_seq"'
        return ph

    def get_physical_setting_decl(self):
        content = '// Configuration parameters declaration, case \'%s\'\n' %(self.name)
        content = content + sim.assign_freefem_var('xmin', self.xmin) + '\n'
        content = content + sim.assign_freefem_var('xmax', self.xmax) + '\n'
        content = content + sim.assign_freefem_var('deltax', self.deltax) + '\n'
        content = content + sim.assign_freefem_var('ub', self.ub) + '\n'
        content = content + sim.assign_freefem_var('d', self.d) + '\n'
        content = content + sim.assign_freefem_var('gamma', self.gamma) + '\n'
        content = content + sim.assign_freefem_var('mu0', self.mu0) + '\n'
        content = content + sim.assign_freefem_var('s', self.lambda0)
        content = content + '\n// End of parameters declaration \n'
        return content
    # --------------------------------------------------------------------------
    # EIGENVALUE PROBLEM
    # --------------------------------------------------------------------------
    def compute_eigs(self):
        content = self.make_eig_edp_file()
        sim.write_file(EIG_EDP, content)


    def make_eig_edp_file(self):
        # Read associated EDP template
        eig_temp   = sim.read_template(TEMPLATE_EIG)
        # Physical setting declaration + initial eigenvalue
        content = self.get_physical_setting_decl()
        # Template
        content = content + eig_temp
        # Replace place holders
        content = sim.replace_placeholders(self.get_placeholders(), content)
        #
        sim.write_file(EIG_EDP, content)
        return content
