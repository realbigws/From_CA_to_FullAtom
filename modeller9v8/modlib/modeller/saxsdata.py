"""Classes to handle SAXS (Small Angle X-ray Scattering) data"""

import _modeller
from modeller.util.modobject import modobject
from modeller.util import modlist, array

__docformat__ = "epytext en"

class SAXSList(modlist.LinkList):
    """A list of L{saxsdata} objects"""
    def __init__(self, edat):
        self.__edat = edat
        self.__list = []
        modlist.LinkList.__init__(self)

    def _insfunc(self, indx, obj):
        _modeller.mod_saxsdata_pt_new(self.__edat, indx, obj.modpt)
        self.__list.insert(indx, obj)
    def __len__(self):
        return len(self.__list)
    def _getfunc(self, indx):
        return self.__list[indx]
    def _delfunc(self, indx):
        del self.__list[indx]
        _modeller.mod_saxsdata_pt_del(self.__edat, indx)


class saxsdata(modobject):
    """Holds all SAXS (Small Angle X-ray Scattering) data"""
    __modpt = None
    env = None

    def __init__(self, env, **vars):
        self.__modpt = _modeller.mod_saxsdata_new()
        self.env = env.copy()

    def __setstate__(self, d):
        self.__dict__.update(d)
        self.__modpt = _modeller.mod_saxsdata_new()

    def __del__(self):
        if self.__modpt:
            _modeller.mod_saxsdata_free(self.__modpt)

    def __get_modpt(self):
        return self.__modpt

    def write(self, file):
        """Write SAXS data, which is currently in memory """
        fh=open(file,'w')
        for ii in range(0,self.ns):
            fh.write('%10.7f  ' % self.s[ii] + '%15.5f  ' % self.intensity[ii] +'%15.5f\n' % self.int_exp[ii]  )
        fh.close()

    def pr_score(self, mdl, maxr, filename=None, use_err=False, rfact=False):
        """Calculates P(r) score of a model by comparing model P(r) to expt
           data saxs.pr_exp.
           @param mdl: model
           @param maxr: maximum radius to score P(r) in A
           @param filename: filename of P(r) model to write (scaled to expt data)
           @param use_err: use experimental error?
           @param rfact: use rfactor as score?
           @return: (pr score,
                     scaling factor of model P(r) to match expt. P(r))"""
        from math import sqrt
        sum_e  = 0.0
        sum_mm = 0.0
        sum_em = 0.0
        sum_ee = 0.0
        mdl.saxs_pr(self, filename='None')
        imaxr = min( int(maxr/self.dr_exp)+1, len(self.p_r_exp))
        if (rfact):
            for ii in range(0, imaxr):
                sum_mm = sum_mm + self.p_r_resamp[ii]
                sum_ee = sum_ee + self.p_r_exp[ii]
            scf = sum_ee/sum_mm
            sum_ee = 0.0
            for ii in range(0, imaxr):
                sum_em = sum_em + abs(self.p_r_exp[ii]-scf*self.p_r_resamp[ii])
                sum_ee = sum_ee + abs(self.p_r_exp[ii])
            psc = sum_em/sum_ee
        else:
            if (use_err):
                for ii in range(0, imaxr):
                    sum_mm = sum_mm + self.p_r_resamp[ii]*self.p_r_resamp[ii]/self.p_r_sig[ii]
                    sum_em = sum_em + self.p_r_exp[ii]*self.p_r_resamp[ii]/self.p_r_sig[ii]
                    sum_ee = sum_ee + self.p_r_exp[ii]*self.p_r_exp[ii]/self.p_r_sig[ii]
            else:
                for ii in range(0, imaxr):
                    sum_mm = sum_mm + self.p_r_resamp[ii]*self.p_r_resamp[ii]
                    sum_em = sum_em + self.p_r_exp[ii]*self.p_r_resamp[ii]
                    sum_ee = sum_ee + self.p_r_exp[ii]*self.p_r_exp[ii]
            norm_e = sqrt(sum_ee)
            scf = sum_em / sum_mm
            scf_norm = scf / norm_e
            #psc = sum_mm*scf*scf + sum_ee - 2.*scf * sum_em
            #psc = psc / (sqrt(sum_ee))
            psc = sum_mm*scf_norm*scf_norm + 1 - 2.*scf_norm * sum_em / norm_e
        if (filename):
            fhandle=open(filename, 'w')
            for ii in range(0, len(self.p_r_exp)):
                tmp = scf*self.p_r_resamp[ii]
                fhandle.write('%10.5f ' % self.r_exp[ii] +
                  ' %15.6f\n' % tmp)
            fhandle.close()
        return (psc, scf)

    def ini_saxs(self, atmsel,
                 filename='$(LIB)/formfactors-int_tab_solvation.lib',
                 s_min=0.0, s_max=2.0, maxs=100, nmesh=100, natomtyp=15,
                 represtyp='heav', wswitch='uniform', s_hybrid=0.0,
                 s_low=0.0, s_hi=2.0, spaceflag='real', rho_solv=0.334,
                 use_lookup=True, nr=5000, dr=0.1, nr_exp=300, dr_exp=1.0,
                 use_offset=False, use_rolloff=False, use_conv=False,
                 mixflag=False, pr_smooth=False):
        """Initialize saxsdata
           @param atmsel: selection of atoms
           @param s_min: minimum frequency in reciprocal space in A^-1
           @param s_max: maximum frequency in reciprocal space in A^-1
           @param maxs: maximum number of frequencies
           @param nmesh: actual number of frequencies (<= maxs)
           @param natomtyp: number of 'atoms', i.e. scattering centers
           @param represtyp: representation: 'heav', 'allh', or 'CA'
           @param filename: filename of the library for formfactors
           @param wswitch: character for filter of scoring function options:
                  'unity', 'sq', or 'hybrid'
           @param s_hybrid: frequency above which $ s^2$ weighting is applied
                  if wswitch='hybrid'
           @param s_low: bandpass filter in A^-1 - lower cutoff
           @param s_hi: bandpass filter in A^-1 - higher cutoff
           @param spaceflag: how should I(s) be computed? 'real' space via P(r)
                  or 'reciprocal'? 'real' is more than a magnitude
                  faster but less accurate for high resolution
           @param rho_solv: electron density of solvent;
                  default=0.334 e/A^-3 (H_2O)
           @param use_lookup: use lookup tables for SINC and COS function -
                  significant increase in speed for 'reciprocal' mode
           @param nr: number of points for P(r) sampling
           @param dr: spacing (sampling) of P(r) in A
           @param nr_exp: number of points for P_exp(r) sampling
           @param dr_exp: spacing (sampling) of P(r) in A
           @param use_offset: allow for additive constant in expt. spectrum
           @param use_rolloff: allow for Gaussian rolloff of model spectrum
           @param use_conv: convolute with nitrogen formfactor to mimic hydr
                  layer
           @param mixflag: different conformations present? implemented for
                  HtpG project
           @param pr_smooth: smoothing of P(r)"""
        (inds, mdl) = atmsel.get_atom_indices()
        return _modeller.mod_saxs_ini(self.modpt, mdl.modpt, inds, s_min,
                                      s_max, maxs, nmesh, natomtyp, represtyp,
                                      filename, wswitch, s_hybrid, s_low, s_hi,
                                      spaceflag, rho_solv, use_lookup, nr, dr,
                                      nr_exp, dr_exp, use_offset, use_rolloff,
                                      use_conv, mixflag, pr_smooth)

    def saxs_read(self, filename):
        """Read in experimental SAXS data"""
        return _modeller.mod_saxs_read(self.modpt, filename)

    def read(self,
             saxsfilename,
             atmsel,
             formfacfilename='$(LIB)/formfactors-int_tab_solvation.lib',
             natomtyp=15,
             represtyp='heav', wswitch='uniform', s_hybrid=0.0,
             s_low=None, s_hi=None,
             spaceflag='real', rho_solv=0.334,
             use_lookup=True, nr=5000, dr=0.1, nr_exp=300, dr_exp=1.0,
             use_offset=False, use_rolloff=False, use_conv=False,
             mixflag=False, pr_smooth=False):
        """Read in experimental SAXS data and initialize saxsdata

           @param saxsfilename: Name of file containing SAXS spectrum
           @param atmsel: selection of atoms
           @param s_min: minimum frequency in reciprocal space in A^-1
           @param s_max: maximum frequency in reciprocal space in A^-1
           @param natomtyp: number of 'atoms', i.e. scattering centers
           @param represtyp: representation: 'heav', 'allh', or 'CA'
           @param formfacfilename: filename of the library for formfactors
           @param wswitch: character for filter of scoring function options:
                  'unity', 'sq', or 'hybrid'
           @param s_hybrid: frequency above which $ s^2$ weighting is applied
                  if wswitch='hybrid'
           @param s_low: bandpass filter in A^-1 - lower cutoff
           @param s_hi: bandpass filter in A^-1 - higher cutoff
           @param spaceflag: how should I(s) be computed? 'real' space via P(r)
                  or 'reciprocal'? 'real' is more than a magnitude
                  faster but less accurate for high resolution
           @param rho_solv: electron density of solvent;
                  default=0.334 e/A^-3 (H_2O)
           @param use_lookup: use lookup tables for SINC and COS function -
                  significant increase in speed for 'reciprocal' mode
           @param nr: number of points for P(r) sampling
           @param dr: spacing (sampling) of P(r) in A
           @param nr_exp: number of points for P_exp(r) sampling
           @param dr_exp: spacing (sampling) of P(r) in A
           @param use_offset: allow for additive constant in expt. spectrum
           @param use_rolloff: allow for Gaussian rolloff of model spectrum
           @param use_conv: convolute with nitrogen formfactor to mimic hydr
                  layer
           @param mixflag: different conformations present? implemented for
                  HtpG project
           @param pr_smooth: smoothing of P(r)"""
        try:
            fh = open(saxsfilename,'r')
        except:
            print "file "+saxsfilename+" not found :("
            return
        fh.close()
        ns = 0
        s_min = 10.
        s_max = 0.
        for line in open(saxsfilename,'r'):
            s = line.split()
            #  '#' is comment
            if (not s[0][0] == '#'):
                ns = ns +1
                if ( float(s[0]) > s_max):
                    s_max = float(s[0])
                if ( float(s[0]) < s_min):
                    s_min = float(s[0])
        if (not s_low):
            s_low = s_min - .001
        if (not s_hi):
            s_hi = s_max + .001
        print "s_min="+str(s_min)+", s_max="+str(s_max)
        print "s_low="+str(s_low)+", s_hi="+str(s_hi)
        self.ini_saxs(atmsel,
                 filename=formfacfilename,
                 s_min=s_min, s_max=s_max, maxs=ns, nmesh=ns, natomtyp=natomtyp,
                 represtyp=represtyp, wswitch=wswitch, s_hybrid=s_hybrid,
                 s_low=s_low, s_hi=s_hi, spaceflag=spaceflag, rho_solv=rho_solv,
                 use_lookup=use_lookup, nr=nr, dr=dr, nr_exp=nr_exp, dr_exp=dr_exp,
                 use_offset=use_offset, use_rolloff=use_rolloff, use_conv=use_conv,
                 mixflag=mixflag, pr_smooth=pr_smooth)
        self.saxs_read(saxsfilename)


    def saxs_pr_read(self, filename):
        """Read in experimental P(r)"""
        return _modeller.mod_saxs_pr_read(self.modpt, filename)

    def __get_s_hybrid(self):
        return _modeller.mod_saxsdata_s_hybrid_get(self.modpt)
    def __set_s_hybrid(self, val):
        return _modeller.mod_saxsdata_s_hybrid_set(self.modpt, val)

    def __get_s_max(self):
        return _modeller.mod_saxsdata_s_max_get(self.modpt)
    def __set_s_max(self, val):
        return _modeller.mod_saxsdata_s_max_set(self.modpt, val)

    def __get_s_min(self):
        return _modeller.mod_saxsdata_s_min_get(self.modpt)
    def __set_s_min(self, val):
        return _modeller.mod_saxsdata_s_min_set(self.modpt, val)

    def __get_s_low(self):
        return _modeller.mod_saxsdata_s_low_get(self.modpt)
    def __set_s_low(self, val):
        return _modeller.mod_saxsdata_s_low_set(self.modpt, val)

    def __get_s_hi(self):
        return _modeller.mod_saxsdata_s_hi_get(self.modpt)
    def __set_s_hi(self, val):
        return _modeller.mod_saxsdata_s_hi_set(self.modpt, val)

    def __get_normsq_exp(self):
        return _modeller.mod_saxsdata_normsq_exp_get(self.modpt)
    def __set_normsq_exp(self, val):
        return _modeller.mod_saxsdata_normsq_exp_set(self.modpt, val)

    def __get_ns(self):
        return _modeller.mod_saxsdata_ns_get(self.modpt)
    def __get_nr(self):
        return _modeller.mod_saxsdata_nr_get(self.modpt)
    def __get_nr_exp(self):
        return _modeller.mod_saxsdata_nr_exp_get(self.modpt)
    def __set_nr_exp(self):
        return _modeller.mod_saxsdata_nr_exp_set(self.modpt, val)
    def __get_dr(self):
        return _modeller.mod_saxsdata_dr_get(self.modpt)
    def __get_dr_exp(self):
        return _modeller.mod_saxsdata_dr_exp_get(self.modpt)
    def __set_dr_exp(self):
        return _modeller.mod_saxsdata_dr_exp_set(self.modpt, val)
    def __get_c(self):
        return _modeller.mod_saxsdata_c_get(self.modpt)
    def __set_c(self, val):
        return _modeller.mod_saxsdata_c_set(self.modpt, val)
    def __get_rolloff(self):
        return _modeller.mod_saxsdata_rolloff_get(self.modpt)
    def __set_rolloff(self):
        return _modeller.mod_saxsdata_rolloff(self.modpt, val)
    def __get_bfac(self):
        return _modeller.mod_saxsdata_bfac_get(self.modpt)
    def __set_bfac(self):
        return _modeller.mod_saxsdata_bfac(self.modpt, val)

    def __get_chi_sq(self):
        return _modeller.mod_saxsdata_chi_sq_get(self.modpt)
    def __set_chi_sq(self, val):
        return _modeller.mod_saxsdata_chi_sq_set(self.modpt, val)

    def __get_rho_solv(self):
        return _modeller.mod_saxsdata_rho_solv_get(self.modpt)
    def __set_rho_solv(self, val):
        return _modeller.mod_saxsdata_rho_solv_set(self.modpt, val)

    def __get_offset(self):
        return _modeller.mod_saxsdata_offset_get(self.modpt)
    def __set_offset(self):
        return _modeller.mod_saxsdata_offset(self.modpt, val)
    def __get_intensity(self):
        ptarr = _modeller.mod_saxsdata_intensity_get(self.modpt)
        return array.Double1DArray(ptarr, self.__get_ns)
    def __set_intensity(self, val):
        modlist.set_fixlist(self.intensity, val)
    def __get_int_exp(self):
        ptarr = _modeller.mod_saxsdata_int_exp_get(self.modpt)
        return array.Double1DArray(ptarr, self.__get_ns)
    def __set_int_exp(self, val):
        modlist.set_fixlist(self.int_exp, val)
    def __get_s(self):
        ptarr = _modeller.mod_saxsdata_s_get(self.modpt)
        return array.Double1DArray(ptarr, self.__get_ns)
    def __set_s(self, val):
        modlist.set_fixlist(self.s, val)
    def __get_sigma_exp(self):
        ptarr = _modeller.mod_saxsdata_sigma_exp_get(self.modpt)
        return array.Double1DArray(ptarr, self.__get_ns)
    def __set_sigma_exp(self, val):
        modlist.set_fixlist(self.sigma_exp, val)
    def __get_p_r(self):
        ptarr = _modeller.mod_saxsdata_p_r_get(self.modpt)
        return array.Double1DArray(ptarr, self.__get_nr)
    def __set_p_r(self):
        modlist.set_fixlist(self.p_r, val)
    def __get_p_r_exp(self):
        ptarr = _modeller.mod_saxsdata_p_r_exp_get(self.modpt)
        return array.Double1DArray(ptarr, self.__get_nr_exp)
    def __set_p_r_exp(self):
        modlist.set_fixlist(self.p_r_exp, val)
    def __get_r_exp(self):
        ptarr = _modeller.mod_saxsdata_r_exp_get(self.modpt)
        return array.Double1DArray(ptarr, self.__get_nr_exp)
    def __set_r_exp(self):
        modlist.set_fixlist(self.r_exp, val)
    def __get_p_r_resamp(self):
        ptarr = _modeller.mod_saxsdata_p_r_resamp_get(self.modpt)
        return array.Double1DArray(ptarr, self.__get_nr_exp)
    def __set_p_r_resamp(self):
        modlist.set_fixlist(self.p_r_resamp, val)
    def __get_p_r_sig(self):
        ptarr = _modeller.mod_saxsdata_p_r_sig_get(self.modpt)
        return array.Double1DArray(ptarr, self.__get_nr_exp)
    def __set_p_r_sig(self):
        modlist.set_fixlist(self.p_r_sig, val)

    modpt = property(__get_modpt)
    s_hybrid = property(__get_s_hybrid, __set_s_hybrid)
    s_max = property(__get_s_max, __set_s_max)
    s_min = property(__get_s_min, __set_s_min)
    s_low = property(__get_s_low, __set_s_low)
    s_hi  = property(__get_s_hi,  __set_s_hi )
    normsq_exp = property(__get_normsq_exp, __set_normsq_exp)
    c     = property(__get_c, __set_c)
    ns = property(__get_ns)
    nr = property(__get_nr)
    nr_exp = property(__get_nr_exp, __set_nr_exp)
    dr = property(__get_dr)
    dr_exp = property(__get_dr_exp, __set_dr_exp)
    intensity = property(__get_intensity, __set_intensity)
    s = property(__get_s, __set_s)
    int_exp = property(__get_int_exp, __set_int_exp)
    sigma_exp = property(__get_sigma_exp, __set_sigma_exp)
    p_r = property(__get_p_r, __set_p_r)
    p_r_exp = property(__get_p_r_exp, __set_p_r_exp)
    r_exp = property(__get_r_exp, __set_r_exp)
    p_r_resamp = property(__get_p_r_resamp, __set_p_r_resamp)
    p_r_sig = property(__get_p_r_sig, __set_p_r_sig)
    p_r_sig = property(__get_p_r_sig, __set_p_r_sig)
    chi_sq = property(__get_chi_sq, __set_chi_sq)
    rho_solv = property(__get_rho_solv, __set_rho_solv)
    rolloff= property(__get_rolloff, __set_rolloff)
    bfac   = property(__get_bfac, __set_bfac)
    offset = property(__get_offset, __set_offset)
