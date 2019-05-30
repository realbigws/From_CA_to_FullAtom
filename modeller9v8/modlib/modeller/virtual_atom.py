import _modeller

from pseudo_atom import pseudo_atom

class ch1(pseudo_atom):
    """Virtual aliphatic proton on a tetrahedral carbon (->CH),
       defined by the central C and the three other substituents"""
    _builtin_index = 2
    numatoms = 4


class ch1a(pseudo_atom):
    """Virtual aromatic proton on a trigonal carbon (=CH),
       defined by the central C and the two C atoms bonded to the central C"""
    _builtin_index = 3
    numatoms = 3


class ch2(pseudo_atom):
    """Virtual aliphatic proton on a tetrahedral carbon (>CH2)
       assigned stereospecifically; defined by the central
       tetrahedral atom and the other two substituents on it"""
    _builtin_index = 5
    numatoms = 3
