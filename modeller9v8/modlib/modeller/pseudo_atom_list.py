from modeller.util.modlist import AppendList
import pseudo_atom, virtual_atom
import _modeller

class PseudoAtomList(AppendList):
    typemap = None

    def __init__(self, mdl):
        AppendList.__init__(self)
        self.__mdl = mdl
        if self.typemap is None:
            t = {}
            self.__add_typemap(t, 'pseudo_atom')
            self.__add_typemap(t, 'virtual_atom')
            PseudoAtomList.typemap = t

    def __add_typemap(self, t, modname):
        classes = dir(eval(modname))
        for c in classes:
            c = eval("%s.%s" % (modname, c))
            try:
                t[c._builtin_index] = c
            except AttributeError:
                pass

    def __len__(self):
        psd = _modeller.mod_model_psd_get(self.__mdl.modpt)
        return _modeller.mod_pseudo_atoms_npseudo_get(psd)

    def _getfunc(self, indx):
        mdl = self.__mdl
        psdtype, inds = _modeller.mod_pseudo_atom_get(mdl.modpt, indx)
        cls = PseudoAtomList.typemap[psdtype]
        p = cls(*[mdl.atoms[i-1] for i in inds])
        p._set_atom_index(mdl.natm + 1 + indx, mdl)
        return p

    def append(self, obj):
        mdl = self.__mdl
        if not isinstance(obj, pseudo_atom.pseudo_atom):
            raise TypeError("can only use pseudo_atom objects here")
        typ = obj.get_type()
        ind = _modeller.mod_pseudo_atom_append(mdl.modpt, typ,
                                               obj._get_base_atoms(mdl))
        obj._set_atom_index(ind, mdl)
