import sqlite3
from sqlite3 import Error

class TimeDomainSimulator:
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
