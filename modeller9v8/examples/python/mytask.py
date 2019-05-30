from modeller import *
from modeller.parallel import task

class MyTask(task):
    """A task to read in a PDB file on the slave, and return the resolution"""
    def run(self, code):
        env = environ()
        env.io.atom_files_directory = ["../atom_files"]
        mdl = model(env, file=code)
        return mdl.resolution
