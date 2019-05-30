from modeller import *
from modeller.scripts import complete_pdb

env = environ()
log.verbose()
env.io.atom_files_directory = ['../atom_files']
env.libs.topology.read(file='$(LIB)/top_heav.lib')
env.libs.parameters.read(file='$(LIB)/par.lib')

class MyTerm(terms.energy_term):
    """Custom energy term, which tries to force all atoms to one side of
       the x=10.0A plane"""

    _physical_type = physical.absposition

    # Override the __init__ function so that we can pass in a 'strength'
    # parameter
    def __init__(self, strength):
        self.strength = strength
        terms.energy_term.__init__(self)

    def eval(self, mdl, deriv, indats):
        atoms = self.indices_to_atoms(mdl, indats)
        e = 0.
        dvx = [0.] * len(indats)
        dvy = [0.] * len(indats)
        dvz = [0.] * len(indats)
        for (num, at) in enumerate(atoms):
            # Enforce a linearly increasing potential in the x direction
            if at.x > 10.0:
                e += (at.x - 10.0) * self.strength
                dvx[num] += self.strength
        if deriv:
            return (e, dvx, dvy, dvz)
        else:
            return e

t = env.edat.energy_terms
t.append(MyTerm(strength=1.0))

mdl = complete_pdb(env, "1fdn")
sel = selection(mdl)
print sel.energy()
