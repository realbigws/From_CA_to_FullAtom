def handle_seq_indx(seqtype, indx, lookup_func=None, lookup_args=(),
                    require_inrange=True):
    if isinstance(indx, int):
        if indx < 0:
            indx += len(seqtype)
        if indx < 0 or indx >= len(seqtype):
            if require_inrange:
                raise IndexError("list index out of range")
            else:
                indx = max(indx, 0)
                indx = min(indx, len(seqtype))
                return indx
        else:
            return indx
    elif isinstance(indx, slice):
        return range(*indx.indices(len(seqtype)))
    elif lookup_func is not None:
        args = lookup_args + (indx,)
        int_indx = lookup_func(*args)
        if int_indx < 0:
            raise KeyError(indx)
        else:
            return int_indx
    else:
        raise TypeError("expecting an integer index")
