"""Utility methods for handling matrices"""

def matrix_to_list(matrix):
    """Flatten a 3x3 matrix to a 9-element list"""
    lst = []
    if isinstance(matrix, (list, tuple)) and len(matrix) == 3:
        for row in matrix:
            if isinstance(row, (list, tuple)) and len(row) == 3:
                lst.extend(row)
            else:
                raise ValueError("Matrix must be a 3x3 array")
    else:
        raise ValueError("Matrix must be a 3x3 array")
    return lst

def list_to_matrix(lst):
    """Construct a 3x3 matrix from a 9-element list"""
    if isinstance(lst, (list, tuple)) and len(lst) == 9:
        return [lst[0:3], lst[3:6], lst[6:9]]
    else:
        raise ValueError("Expecting a 9-element list to construct a matrix")
