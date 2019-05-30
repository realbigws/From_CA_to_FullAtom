# Homology modeling with ligand transfer from the template
from modeller import *              # Load standard Modeller classes
from modeller.automodel import *    # Load the automodel class

log.verbose()    # request verbose output
env = environ()  # create a new MODELLER environment to build this model in

# directories for input atom files
env.io.atom_files_directory = ['.', '../atom_files']

# Read in HETATM records from template PDBs
env.io.hetatm = True

a = automodel(env,
              alnfile  = 'align-ligand.ali',  # alignment filename
              knowns   = '5fd1',              # codes of the templates
              sequence = '1fdx')              # code of the target
a.starting_model= 4                 # index of the first model
a.ending_model  = 4                 # index of the last model
                                    # (determines how many models to calculate)
a.make()                            # do the actual homology modeling
