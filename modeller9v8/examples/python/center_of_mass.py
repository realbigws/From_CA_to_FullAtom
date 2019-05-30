from modeller import *
from modeller.scripts import complete_pdb
import math

env = environ()

env.io.atom_files_directory = ['../atom_files']
log.verbose()
env.libs.topology.read(file='$(LIB)/top_heav.lib')
env.libs.parameters.read(file='$(LIB)/par.lib')

mdl = complete_pdb(env, "1fdn")

def center_of_mass(mdl):
    """Calculates the center of mass and total mass of the whole model"""
    totalmass = 0.
    cx = cy = cz = 0.
    for a in mdl.atoms:
        totalmass += a.mass
        cx += a.x * a.mass
        cy += a.y * a.mass
        cz += a.z * a.mass
    cx /= totalmass
    cy /= totalmass
    cz /= totalmass
    return (cx, cy, cz, totalmass)

def radius_gyration(mdl, cx, cy, cz, totalmass):
    """Returns the radius of gyration of the whole molecule, given the
       center of mass and total mass"""
    rgyr = 0.
    for a in mdl.atoms:
        rgyr += ((a.x - cx) ** 2 + (a.y - cy) ** 2 + (a.z - cz) ** 2) * a.mass
    rgyr /= totalmass
    return math.sqrt(rgyr)

(cx, cy, cz, totalmass) = center_of_mass(mdl)
print "COM %f, %f, %f" % (cx, cy, cz)

rgyr = radius_gyration(mdl, cx, cy, cz, totalmass)
print "Rgyr = %f" % rgyr
