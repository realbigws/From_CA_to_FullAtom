# Very fast homology modeling by the automodel class
from modeller import *
from modeller.automodel import *    # Load the automodel class

log.verbose()
env = environ()
# directories for input atom files
env.io.atom_files_directory = ['.', '../atom_files']

a = automodel(env,
              alnfile='alignment.ali',      # alignment filename
              knowns='5fd1',                # codes of the templates
              sequence='1fdx',              # code of the target
              assess_methods=assess.GA341)  # request GA341 model assessment

a.very_fast()                       # prepare for extremely fast optimization

a.starting_model = 2
a.ending_model = 2
a.final_malign3d = True

a.make()                            # make the homology model
