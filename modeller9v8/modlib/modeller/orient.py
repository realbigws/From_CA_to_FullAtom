from modeller.util.matrix import list_to_matrix

class OrientData(object):
    """Data returned from the 'model.orient' method"""

    def __init__(self, rotation, translation):
        self.rotation = list_to_matrix(rotation)
        self.translation = translation
