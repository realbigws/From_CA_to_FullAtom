"""Classes for obtaining per-residue energy profiles."""

__docformat__ = "epytext en"

def __smooth_one(vals, n, window):
    """Return the weighted average of C{window} values around C{vals[n]}"""
    i1 = max(0, n - window)
    i2 = min(len(vals) - 1, n + window)
    wtsum = av = 0.
    for j in range(i1, i2 + 1):
        weight = 0.1 * (window - abs(j - n) + 1)
        av += weight * vals[j]
        wtsum += weight
    return av / max(wtsum, 0.1)

def _smooth(vals, window):
    """Return a smoothed copy of the array C{vals}"""
    return [__smooth_one(vals, n, window) for n in range(len(vals))]

class EnergyProfile(object):
    """A per-residue energy profile"""

    def __init__(self, profile, nprofile, min_rms, heavy_rms):
        self.profile = profile
        self.nprofile = nprofile
        self.min_rms = min_rms
        self.heavy_rms = heavy_rms

    def __repr__(self):
        return repr(self.profile)

    def get_normalized(self):
        """Return a new 'normalized' energy profile, in which each residue's
           energy is divided by the number of restraints acting on that
           residue."""
        return EnergyProfile([prof/nprof for prof, nprof in \
                                         zip(self.profile, self.nprofile)],
                             [1]*len(self.profile), self.min_rms,
                             self.heavy_rms)

    def get_smoothed(self, window=1):
        """Return a new energy profile, smoothed by window averaging."""
        return EnergyProfile(_smooth(self.profile, window), self.nprofile,
                             self.min_rms, self.heavy_rms)

    def write_to_file(self, filename):
        """Write the profile to a file"""
        fh = file(filename, "w")
        for (n, val) in enumerate(self.profile):
            print >> fh, "%10d %12.4f" % (n+1, val)
