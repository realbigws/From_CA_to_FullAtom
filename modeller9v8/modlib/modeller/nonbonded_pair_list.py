import _modeller

class NonbondedPairList(object):
    def __init__(self, mdl):
        self.__mdl = mdl

    def __len__(self):
        return _modeller.mod_model_npairs_get(self.__mdl.modpt)
