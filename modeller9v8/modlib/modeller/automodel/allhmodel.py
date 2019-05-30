"""Classes to build all-atom model(s) using template information"""

from automodel import automodel

__docformat__ = "epytext en"

class allhmodel(automodel):
    """Automatically build all-atom model(s) using template information"""

    toplib = '${LIB}/top_allh.lib'

    def __init__(self, env, alnfile, knowns, sequence, deviation=None,
                 library_schedule=None, csrfile=None, inifile=None,
                 assess_methods=None):
        automodel.__init__(self, env, alnfile, knowns, sequence, deviation,
                           library_schedule, csrfile, inifile, assess_methods)
        # Modeling won't work unless we read/write hydrogen atoms!
        self.env.io.hydrogen = True
