from modeller.util.matrix import list_to_matrix

class SuperposeData(object):
    """Data returned from the 'selection.superpose' method"""

    def __init__(self, initial_rms, rms, drms, rotation, translation,
                 num_equiv_pos, num_equiv_dist, num_equiv_cutoff_pos,
                 num_equiv_cutoff_dist, cutoff_rms, cutoff_drms):
        self.initial_rms = initial_rms
        self.rms = rms
        self.drms = drms
        self.rotation = list_to_matrix(rotation)
        self.translation = translation
        self.num_equiv_pos = num_equiv_pos
        self.num_equiv_dist = num_equiv_dist
        self.num_equiv_cutoff_pos = num_equiv_cutoff_pos
        self.num_equiv_cutoff_dist = num_equiv_cutoff_dist
        self.cutoff_rms = cutoff_rms
        self.cutoff_drms = cutoff_drms
