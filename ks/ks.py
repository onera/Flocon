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
MESH_FILE           = os.path.join(MESH_DIR, 'mesh.msh')
# Templates for the .edp files
TEMPLATE_DIR        = os.path.join(CLASS_DIR, 'edp_templates')
TEMPLATE_EIG        = os.path.join(TEMPLATE_DIR, 'eig')
TEMPLATE_MESH       = os.path.join(TEMPLATE_DIR, 'mesh_generator')
TEMPLATE_LINE       = os.path.join(TEMPLATE_DIR, 'get_line')
# edp files
EIG_EDP             = os.path.join(WORKING_DIR, 'eig.edp')
MESH_GENERATOR_EDP  = os.path.join(WORKING_DIR, 'mesh_generator.edp')
#
EIG_FILE            = os.path.join(WORKING_DIR, 'eig')
EV_FILE             = os.path.join(WORKING_DIR, 'ev')
X0_FILE             = os.path.join(WORKING_DIR, 'x0')

#
LOG_FILE            = os.path.join(WORKING_DIR,'log')


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
        return '%.2f,%d,%.3f,%r'%(self.mu0, self.N, self.dt, self.NL)

    def print_msg(self, msg):
        print('[FLOCON][KS]['+self.reduced_id + '] '+ msg)

    def get_placeholders(self):
        ph = {'MESH':MESH_FILE,
              'EIG_FILE': EIG_FILE,
              'EV_FILE':EV_FILE,
              'SOLVER':''}
        if self.solver == 'mumps':
            # ph['SOLVER'] = 'load "MUMPS_seq"\ndefaulttoMUMPSseq();'
            ph['SOLVER'] = 'load "MUMPS_seq"'
        return ph

    def get_geom_decl(self):
        content = '// Parameters defining the geometry, case \'%s\'\n' %(self.name)
        content = content + sim.assign_freefem_var('xmin', self.xmin) + '\n'
        content = content + sim.assign_freefem_var('xmax', self.xmax) + '\n'
        content = content + sim.assign_freefem_var('deltax', self.deltax) + '\n'
        return content

    def get_physical_setting_decl(self):
        content = self.get_geom_decl()
        content = content + '// Other parameters, case \'%s\'\n' %(self.name)
        content = content + sim.assign_freefem_var('ub', self.ub) + '\n'
        content = content + sim.assign_freefem_var('d', self.d) + '\n'
        content = content + sim.assign_freefem_var('gamma', self.gamma) + '\n'
        content = content + sim.assign_freefem_var('mu0', self.mu0) + '\n'
        content = content + sim.assign_freefem_var('s', self.lambda0)
        content = content + '\n// End of parameters declaration \n'
        return content
    # --------------------------------------------------------------------------
    # GENERATING MESH
    # --------------------------------------------------------------------------
    def make_mesh(self):
        self.mesh_generator_edp()
        sim.launch_edp_file(MESH_GENERATOR_EDP, log=LOG_FILE)
        if os.path.exists(MESH_FILE):
            self.print_msg('Mesh file successfully generated.')

    def mesh_generator_edp(self):
        generator_temp = sim.read_template(TEMPLATE_MESH)
        content = self.get_geom_decl()
        content = content + generator_temp
        content = sim.replace_placeholders(self.get_placeholders(), content)
        sim.write_file(MESH_GENERATOR_EDP, content)
        return content
    # --------------------------------------------------------------------------
    # INITIAL STATE
    # --------------------------------------------------------------------------
    def make_x0(self):
        [eigs, evs] = self.compute_eigs()
        wmax        = np.max(np.abs(np.real(evs[0])))
        x0          = np.array(0.01/ wmax * np.real(evs[0]))
        #
        sim.np_to_freefem_file(X0_FILE, x0)
        return x0
    # --------------------------------------------------------------------------
    # EIGENVALUE PROBLEM
    # --------------------------------------------------------------------------
    def compute_eigs(self, nev=2):
        content = self.make_eig_edp_file(nev)
        sim.write_file(EIG_EDP, content)
        sim.launch_edp_file(EIG_EDP, log=LOG_FILE)
        # Read eigenvalues and eigenvectors
        eigs_ri = sim.freefem_data_file_to_np(EIG_FILE, 0)
        evs     = []
        eigs    = np.zeros((nev,1), dtype=np.complex)
        for i in range(nev):
            eigs[i] = complex(eigs_ri[i][0], eigs_ri[i][1])
            tmp     = sim.freefem_cvec_to_np('%s%d'%(EV_FILE,i))
            evs.append(tmp)
        return [eigs, evs]

    def make_eig_edp_file(self, nev=2):
        # Read associated EDP template
        eig_temp    = sim.read_template(TEMPLATE_EIG)
        # Physical setting declaration + initial eigenvalue
        content     = self.get_physical_setting_decl()
        content     = content + sim.assign_freefem_var('nev', nev) + '\n'
        # Template
        content     = content + eig_temp
        # Replace place holders
        content     = sim.replace_placeholders(self.get_placeholders(), content)
        #
        sim.write_file(EIG_EDP, content)
        return content
