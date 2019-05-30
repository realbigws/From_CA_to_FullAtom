import _modeller
from modeller.error import ModellerError

class action(object):
    """Base class for all periodic optimizer actions"""
    def __init__(self, skip, first, last):
        self.skip = skip
        self.first = first
        self.last = last

class write_structure(action):
    """Write out the current model coordinates"""

    def __init__(self, skip, filepattern, write_all_atoms=True, first=False,
                 last=False, start=0):
        action.__init__(self, skip, first, last)
        self.num = start
        self.filepattern = filepattern
        if not isinstance(filepattern, str):
            raise TypeError("filepattern must be a string")
        self.write_all_atoms = write_all_atoms

    def __call__(self, opt):
        atoms = opt.get_selection()
        try:
            filename = self.filepattern % self.num
        except TypeError:
            filename = self.filepattern + "-%d" % self.num
        if self.write_all_atoms:
            atoms = atoms.get_model()
        atoms.write(file=filename)
        self.num = self.num + 1


class trace(action):
    """Write out optimization energies, etc."""

    def __init__(self, skip, output=None):
        action.__init__(self, skip, True, True)
        if output is None:
            import sys
            output = sys.stdout
        elif not hasattr(output, 'write'):
            output = file(output, 'w')
        self.output = output

    def __call__(self, opt):
        opt.trace(self.output)


class charmm_trajectory(action):
    """Write out a trajectory in CHARMM (DCD) format"""
    __unit = None

    def __init__(self, skip, filename, first=False, last=False):
        action.__init__(self, skip, first, last)
        self.num = 0
        self.filename = filename

    def __del__(self):
        if self.__unit != None:
            _modeller.mod_traj_close(self.__unit)

    def __call__(self, opt):
        atmsel = opt.get_selection()
        (inds, mdl) = atmsel.get_atom_indices()
        if self.__unit is None:
            self.__unit = _modeller.mod_traj_write_header(mdl.modpt,
                                                          self.filename,
                                                          self.skip, inds)
            self.__inds = inds
            self.__mdl = mdl
        else:
            if mdl is not self.__mdl:
                raise ModellerError("all CHARMM trajectory frames must " +
                                    "be for the same model")
            if inds != self.__inds:
                raise ModellerError("all CHARMM trajectory frames must " +
                                    "be for the same set of atoms")
            _modeller.mod_traj_write_set(mdl.modpt, inds, self.__unit)
