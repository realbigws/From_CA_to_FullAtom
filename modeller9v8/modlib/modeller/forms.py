import _modeller
import physical

class restraint_form(object):
    """Base class for all restraint forms"""

    _builtin_index = None
    _user_index = None
    _use_array = False
    ndim = None

    def __init__(self, group, features, modality, parameters):
        if '_builtin_index' not in self.__class__.__dict__ \
           and self._user_index is None:
            self.__class__._user_index = \
                _modeller.mod_user_form_new(self.eval, self.vmin, self.vheavy,
                                            self.rvmin, self.rvheavy,
                                            self.min_mean, self.heavy_mean)
        if not isinstance(group, physical.physical_type):
            raise TypeError("group should be a physical_type object")
        self._group = group
        self._features = features
        self._modality = modality
        self._parameters = parameters

    def __get_restraint(self, mdl):
        feats = self._features
        if not isinstance(feats, (tuple, list)):
            feats = (feats,)
        modal = self._modality
        if not isinstance(modal, (tuple, list)):
            modal = (modal,)
        if self.ndim is not None and \
           (self.ndim != len(modal) or self.ndim != len(feats)):
            raise ValueError("Incorrect number of dimensions - did you "
                             "forget to call add_dimension() ?")
        iftyp = []
        atoms = []
        natoms = []
        for ft in feats:
            iftyp.append(ft.get_type())
            (at, ftmdl) = ft.get_atom_indices()
            if ftmdl is None:
                raise ValueError("Feature contains no atoms!")
            elif ftmdl != mdl:
                raise ValueError("Feature is defined on a different model!")
            natoms.append(len(at))
            atoms.extend(at)
        if isinstance(self._use_array, bool):
            prm = self._parameters
            use = self._use_array
        else:
            prm = (self._use_array,)
            use = False
        return (self.get_type(), self._group._num, iftyp, modal, natoms, atoms,
                prm, use)

    def _add_restraint(self, rsr, mdl):
        rsrdata = self.__get_restraint(mdl)
        ret = _modeller.mod_restraints_add(rsr.modpt, *rsrdata)
        if ret > 0:
            return ret

    def get_type(cls):
        """Get the numeric type of this form"""
        if cls._user_index is not None:
            return cls._user_index
        elif cls._builtin_index is not None:
            return cls._builtin_index
        else:
            raise ValueError("Cannot get type - no objects created yet")
    get_type = classmethod(get_type)

    def deltaf(self, feat, mean, iftyp):
        return _modeller.mod_feature_delta(feat, mean, iftyp)

    def eval(self, feats, iftyp, modal, deriv, param):
        raise NotImplementedError("Value not defined for form " + str(self))

    def vmin(self, feats, iftyp, modal, param):
        raise NotImplementedError("Min. violation not defined for form "
                                  + str(self))

    def vheavy(self, feats, iftyp, modal, param):
        raise NotImplementedError("Heavy violation not defined for form "
                                  + str(self))

    def rvmin(self, feats, iftyp, modal, param):
        raise NotImplementedError("Relative min. violation not defined " +
                                  "for form " + str(self))

    def rvheavy(self, feats, iftyp, modal, param):
        raise NotImplementedError("Relative heavy violation not defined " +
                                  "for form " + str(self))

    def min_mean(self, feats, iftyp, modal, param):
        raise NotImplementedError("Min. mean not defined for form " + str(self))

    def heavy_mean(self, feats, iftyp, modal, param):
        raise NotImplementedError("Heavy mean not defined for form "
                                  + str(self))


class lower_bound(restraint_form):
    """Left Gaussian (harmonic lower bound)"""
    _builtin_index = 1
    def __init__(self, group, feature, mean, stdev):
        restraint_form.__init__(self, group, feature, 0, (mean, stdev))


class upper_bound(restraint_form):
    """Right Gaussian (harmonic upper bound)"""
    _builtin_index = 2
    def __init__(self, group, feature, mean, stdev):
        restraint_form.__init__(self, group, feature, 0, (mean, stdev))


class gaussian(restraint_form):
    """Single Gaussian (harmonic potential)"""
    _builtin_index = 3
    def __init__(self, group, feature, mean, stdev):
        restraint_form.__init__(self, group, feature, 0, (mean, stdev))


class multi_gaussian(restraint_form):
    """Multiple Gaussian"""
    _builtin_index = 4
    def __init__(self, group, feature, weights, means, stdevs):
        lv = -1
        for var in (weights, means, stdevs):
            if (lv >= 0 and lv != len(var)) \
               or not isinstance(var, (tuple, list)):
                raise TypeError("weights, means and stdevs should all be "
                                + "sequences of the same length")
            lv = len(var)
        restraint_form.__init__(self, group, feature, len(weights),
                                tuple(weights) + tuple(means) + tuple(stdevs))

class lennard_jones(restraint_form):
    """Lennard-Jones potential"""
    _builtin_index = 5
    def __init__(self, group, feature, A, B):
        restraint_form.__init__(self, group, feature, 0, (A, B))


class coulomb(restraint_form):
    """Coulomb potential"""
    _builtin_index = 6
    def __init__(self, group, feature, q1, q2):
        restraint_form.__init__(self, group, feature, 0, (q1, q2))


class cosine(restraint_form):
    """Cosine potential"""
    _builtin_index = 7
    def __init__(self, group, feature, phase, force, period):
        restraint_form.__init__(self, group, feature, period, (phase, force))


class factor(restraint_form):
    """Simple scaling of feature value"""
    _builtin_index = 8
    def __init__(self, group, feature, factor):
        restraint_form.__init__(self, group, feature, 0, (factor,))


class multi_binormal(restraint_form):
    """Multiple binormal"""
    _builtin_index = 9
    def __init__(self, group, features, weights, means, stdevs, correls):
        lv = -1
        for var in (weights, means, stdevs, correls):
            if (lv >= 0 and lv != len(var)) \
               or not isinstance(var, (tuple, list)):
                raise TypeError("weights, means, stdevs and correls should "
                                + "all be sequences of the same length")
            lv = len(var)
        for var in (means, stdevs):
            for elem in var:
                if not isinstance(elem, (tuple, list)) or len(elem) != 2:
                    raise TypeError("Each element of means and stdevs should "
                                    + "be a list of two elements")
        if not isinstance(features, (tuple, list)) or len(features) != 2:
            raise TypeError("features should be a list of two features")
        params = []
        for m in means:
            params += list(m)
        for m in stdevs:
            params += list(m)
        restraint_form.__init__(self, group, features, [len(weights)]*2,
                                tuple(weights) + tuple(params) + tuple(correls))

class spline(restraint_form):
    """Cubic spline"""
    _builtin_index = 10
    def __init__(self, group, feature, open, low, high, delta, lowderiv,
                 highderiv, values, use_array=False):
        self._use_array = use_array
        if not isinstance(values, (tuple, list)):
            raise TypeError("values must be a sequence of spline points")
        spltyp = {False:-1.0, True: 0.0}
        restraint_form.__init__(self, group, feature, len(values),
                                (spltyp[open], low, high, delta, lowderiv,
                                 highderiv) + tuple(values))

class nd_spline(restraint_form):
    """Multi-dimensional cubic spline"""

    _builtin_index = 10

    def __init__(self, group, values, dimensions=None, use_array=False):
        self._use_array = use_array
        modal = []
        param = []
        if dimensions is not None:
            modal = dimensions[:]
            nparam = 1
            for dim in modal:
                nparam *= dim
            if nparam != len(values):
                raise ValueError("Number of spline parameters is not equal" +
                                 " to the product of supplied dimensions")
            param = list(values)
        else:
            self._parse_values(values, modal, param, 0)
        restraint_form.__init__(self, group, [], modal, param)
        self.ndim = 0

    def add_dimension(self, feature, open, low, high, delta, lowderiv,
                      highderiv):
        """Set up the next dimension of the spline"""
        spltyp = {False:-1.0, True: 0.0}
        self._features.append(feature)
        self._parameters[self.ndim*6:self.ndim*6] = (spltyp[open], low, high,
                                                     delta, lowderiv, highderiv)
        self.ndim = self.ndim + 1

    def _parse_values(self, values, modal, param, idim):
        if len(modal) <= idim:
            modal.append(len(values))
        else:
            if modal[idim] != len(values):
                raise ValueError("Spline parameters should be a " +
                                 "rectangular matrix")
        for elmnt in values:
            if hasattr(elmnt, '__iter__'):
                self._parse_values(elmnt, modal, param, idim + 1)
            else:
                if idim + 1 != len(modal):
                    raise ValueError("Spline parameters should be a " +
                                     "rectangular matrix")
                param.append(elmnt)
