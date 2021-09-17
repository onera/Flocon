import sqlite3
from sqlite3 import Error
import io
import numpy as np
import subprocess
import os


CLASS_DIR           = os.path.dirname(__file__)
WORKING_DIR         = os.path.join(CLASS_DIR, 'work')
# Create working temp directory if it does not exists
if not os.path.exists(WORKING_DIR):
    os.mkdir(WORKING_DIR)
#
TMP_EDP_FILE        = os.path.join(WORKING_DIR,'tmp.edp')
LOG_FILE            = os.path.join(WORKING_DIR,'log')

def adapt_array(arr):
    return arr.tobytes()

def convert_array(text):
    return np.frombuffer(text)

sqlite3.register_adapter(np.ndarray, adapt_array)
sqlite3.register_converter("array", convert_array)


class Simulator:
    def __init__(self):
        self.DATABASE_FILE  = []
        self.db             = []

    def init_database(self):
        try:
            self.db = sqlite3.connect(self.DATABASE_FILE,detect_types=sqlite3.PARSE_DECLTYPES)
        except Error as e:
            print(e)

    def init_default_parameters(self):
        pass

    def main_init(self):
        # Database initialisation
        self.init_database()
        # Parameters for the step benchmark
        self.init_default_parameters()

    def __del__(self):
        self.db.close()

    def table_exists(self, table_name):
        c = self.db.cursor()
        c.execute("SELECT count(name) FROM sqlite_master WHERE type='table' AND name=? ", (table_name,))
        return c.fetchone()[0] == 1

    def create_table(self, table_name, elements):
        if self.table_exists(table_name):
            return
        c = self.db.cursor()
        # Insecure...
        cmd = 'CREATE TABLE %s (' %(table_name)
        cmd = cmd + ','.join(elements) +');'
        #
        c.execute(cmd)
        self.db.commit()

class TimeDomainSimulator(Simulator):
    def __init__(self):
        pass

    def test_freefem_solver(self):
        content = '\n'.join(['try',\
                             '{',\
                             'load "MUMPS_seq"',\
                             '// defaulttoMUMPSseq();',\
                             'cout << "SOLVER FOUND" << endl;',\
                             '}',\
                             'catch(...){',\
                             'cout << "SOLVER NOT FOUND" << endl;',\
                             '}',\
                             ])
        #
        write_file(TMP_EDP_FILE, content)
        launch_edp_file(TMP_EDP_FILE, log=LOG_FILE, opt='ne')
        f   = open(LOG_FILE, 'r')
        str = f.read()
        id  = str.find('SOLVER FOUND')
        if id>=0:
            self.solver = 'mumps'
        else:
            self.print_msg('FreeFem "MUMPS" solver not found, switching to default.')
            self.solver='default'
        f.close()
# ------------------------------------------------------------------------------
# External file
# ------------------------------------------------------------------------------
def write_file(target, content):
    f       = open(target, 'w+')
    f.write(content)
    f.close()

def read_template(target):
    return file_to_str(target)

def file_to_str(target):
    f       = open(target, 'r')
    str     = f.read()
    f.close()
    return str
# ------------------------------------------------------------------------------
# EDP interface
# ------------------------------------------------------------------------------
def launch_edp_file(target, opt='', log=None):
    cmd = find_ex('FreeFem++')
    if cmd is None:
        raise Exception("FreeFem++ not found.")
    if opt:
        opt = '-'+opt
    if log is None:
        log ='log'
    if os.name == 'nt':
        subprocess.run('"%s" %s -f %s > %s'%(cmd, opt, target,log),shell=True)
    else:
        subprocess.run('%s %s %s > %s'%(cmd, opt, target,log),shell=True)


def ABCD_to_freefem_file(target, ABCD):
    np_to_freefem_file(target,ABCD)

def np_to_freefem_file(target, m):
    np.savetxt(target,m, delimiter=' ', fmt='%.16f', header='%d %d'%(m.shape[0], m.shape[1]),comments='')

def freefem_data_file_to_np(file):
    return np.loadtxt(file,skiprows=1)

def assign_freefem_var(name, val):
    if type(val) == bool:
        typ = 'bool'
        val = '{}'.format(val)
        val = val.lower()
    elif type(val) == int:
        typ = 'int'
    elif type(val) == float:
        typ = 'real'
    elif type(val) == complex:
        typ = 'complex'
        if np.imag(val) >=0:
            sign = '+'
        else:
            sign = '-'
        return '{} {} = {} {} {}i;\n'.format(typ, name, np.real(val), sign, abs(np.imag(val)) )
    else:
        return ''
    return '{} {} = {};\n'.format(typ, name, val)

def replace_placeholders(ph_list, str):
    for ph_name, path_str in ph_list.items():
        # path_str = os.path.join('.',*ph_val)
        if os.name=='nt':
            path_str = path_str.replace('\\','\\\\')
        str = str.replace('@'+ph_name, path_str)
    return str


def find_ex(ex):
    path    = os.environ['PATH']
    paths   = path.split(os.pathsep)
    ext     = ''
    if os.name == 'win32' or os.name=='nt':
        ext = '.exe'
    exfile = ex + ext
    for p in paths:
        f = os.path.join(p, exfile)
        if os.path.isfile(f):
            return f
    return None
