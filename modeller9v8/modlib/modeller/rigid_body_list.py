from modeller.util.modlist import VarList
from modeller.rigid_body import rigid_body
import _modeller

class RigidBodyList(VarList):
    def __init__(self, mdl):
        self.__mdl = mdl
        VarList.__init__(self)

    def __len__(self):
        return _modeller.mod_model_nrigid_get(self.__mdl.modpt)
    def _setdimfunc(self, num):
        _modeller.mod_model_nrigid_set(self.__mdl.modpt, num)
    def _getfunc(self, indx):
        mdl = self.__mdl
        (inds, scale) = _modeller.mod_rigid_body_get(mdl.modpt, indx)
        p = rigid_body(*[mdl.atoms[i-1] for i in inds])
        p.scale_factor = scale
        return p
    def _setfunc(self, indx, obj):
        mdl = self.__mdl
        if not isinstance(obj, rigid_body):
            raise TypeError("can only use rigid_body objects here")
        _modeller.mod_rigid_body_set(mdl.modpt, indx, obj._get_base_atoms(mdl),
                                     obj.scale_factor)
