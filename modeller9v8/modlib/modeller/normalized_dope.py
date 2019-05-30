"""Classes to get normalized DOPE scores."""


class DOPEScorer(object):
    """Use raw DOPE scores to get normalized measures."""

    aa_scores = { 'ALA': -138.432,
                  'CYS': -161.779,
                  'ASP': -10.0430,
                  'GLU': -23.1399,
                  'PHE': -292.925,
                  'GLY': -12.7251,
                  'HIS': -106.559,
                  'ILE': -237.622,
                  'LYS': -46.4323,
                  'LEU': -238.801,
                  'MET': -195.451,
                  'ASN': -4.41923,
                  'PRO': -2.17275,
                  'GLN': -49.1654,
                  'ARG': -52.7026,
                  'SER': -37.4127,
                  'THR': -49.1835,
                  'VAL': -177.472,
                  'TRP': -285.274,
                  'TYR': -236.675
                }

    def __init__(self, mdl):
        self.__mdl = mdl


    def get_z_score(self, value):
        """Given a DOPE score, return a normalized score (z value)"""
        (mean, numstdres) = self._get_mean(self.__mdl)
        return (value - mean) / self._get_stdev(numstdres)


    def get_profile(self, prof):
        """Given a raw DOPE profile, return a normalized profile"""
        from modeller.energy_profile import EnergyProfile
        mdl = self.__mdl
        numstdres = 0
        norm_prof = [0.] * len(mdl.residues)
        for (n, r) in enumerate(mdl.residues):
            try:
#              No need to divide by 2 here, since profile contributions are
#              already divided by the number of atoms (2 in the case of
#              2-atom statistical potentials)
                norm_prof[n] = prof.profile[n] - self.aa_scores[r.name]
                numstdres += 1
            except KeyError:
                pass
        correction = 2194.27 / float(numstdres)
        for n in range(len(mdl.residues)):
            norm_prof[n] -= correction
        return EnergyProfile(norm_prof, prof.nprofile, prof.min_rms,
                             prof.heavy_rms)

    def _get_stdev(self, numstdres):
        """Get the predicted DOPE score standard deviation"""
        return 14.12 * numstdres

    def _get_mean(self, mdl):
        """Get the predicted DOPE score"""
        (score, numstdres) = self._get_aa_scores(mdl)
        return (score + self._get_small_chain_correction(numstdres), numstdres)

    def _get_aa_scores(self, mdl):
        """Get the contribution to the DOPE score from standard amino acids"""
        score = 0.
        numstdres = 0
        for r in mdl.residues:
            try:
                score += self.aa_scores[r.name]
                numstdres += 1
            except KeyError:
                pass
        return (score, numstdres)

    def _get_small_chain_correction(self, numstdres):
        """DOPE score is not linear for very short chains, so calculate a
           correcting term."""
        e = 2.7182818284590451  # we don't want to require the math module for e
        return 3689.07 * (1.0 - e ** (-numstdres / 70.0 ))


class DOPEHRScorer(DOPEScorer):
    """Use raw DOPE-HR scores to get normalized measures, using the new
       training set."""

    aa_scores = { 'ALA': -160.524,
                  'CYS': -30.3241,
                  'ASP':  34.539,
                  'GLU': -45.4734,
                  'PHE': -131.283,
                  'GLY': 14.7563,
                  'HIS': -135.9,
                  'ILE': -160.324,
                  'LYS': -54.1812,
                  'LEU': -145.011,
                  'MET': -3.51822,
                  'ASN': 15.1741,
                  'PRO': 37.3899,
                  'GLN': -109.93,
                  'ARG': -35.056,
                  'SER': 9.42672,
                  'THR': -77.0285,
                  'VAL': -75.3469,
                  'TRP': -44.7595,
                  'TYR': -306.688
                }

    def _get_stdev(self, numstdres):
        """Get the predicted DOPE score standard deviation"""
        return 19.55 * numstdres

    def _get_small_chain_correction(self, numstdres):
        """Add a constant term to the mean for the fit"""
        return 322.112
