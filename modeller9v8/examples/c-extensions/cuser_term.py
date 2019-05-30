from modeller import *
from modeller.scripts import complete_pdb

# Load the C extension module; this needs to be compiled first - see
# cuser_feat.py for suitable commands.
import _cuser_term


env = environ()
log.verbose()
env.io.atom_files_directory = ['../atom_files']
env.libs.topology.read(file='$(LIB)/top_heav.lib')
env.libs.parameters.read(file='$(LIB)/par.lib')

class MyTerm(terms.energy_term):
    """Custom energy term, which tries to force all atoms to one side of
       the x=10.0A plane, implemented as a C extension."""

    _physical_type = physical.absposition

    # Override the __init__ function so that we can pass in a 'strength'
    # parameter
    def __init__(self, strength):
        self.strength = strength
        terms.energy_term.__init__(self)
    def _add_term(self, edat, indx):
        _cuser_term.myterm_create(edat, indx, self._physical_type.get_type(),
                                  self.strength)

t = env.edat.energy_terms
t.append(MyTerm(strength=1.0))

mdl = complete_pdb(env, "1fdn")
sel = selection(mdl)
print sel.energy()
