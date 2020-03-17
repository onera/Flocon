import sqlite3
from sqlite3 import Error
import io
import numpy as np
import os

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
            self.db = sqlite3.connect(self.DATABASE_FILE)
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
    pass

# ------------------------------------------------------------------------------
# External file
# ------------------------------------------------------------------------------
def write_file(target, content):
    f       = open(target, 'w+')
    f.write(content)
    f.close()

def read_template(target):
    f       = open(target, 'r')
    str     = f.read()
    f.close()
    return str
# ------------------------------------------------------------------------------
# EDP interface
# ------------------------------------------------------------------------------
def launch_edp_file(target):
    os.system('FreeFem++ -nw {} > log'.format(target))

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
