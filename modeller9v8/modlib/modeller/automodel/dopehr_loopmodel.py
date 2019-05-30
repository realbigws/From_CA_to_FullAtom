"""Classes for loop modeling using the DOPE-HR potential"""

from dope_loopmodel import dope_loopmodel
from modeller import group_restraints

__docformat__ = "epytext en"

class dopehr_loopmodel(dope_loopmodel):
    """Loop modeling using the DOPE-HR potential"""

    def read_potential(self):
        return group_restraints(self.env, classes='$(LIB)/atmcls-mf.lib',
                                parameters='$(LIB)/dist-mfhr.lib')
