import modeller

class SalignData(object):
    """Data returned from the 'alignment.salign' method"""

    def __init__(self, aln_score, qscorepct):
        self.aln_score = aln_score
        self.qscorepct = qscorepct

def _salign_fw_local_gaps1(aln, feature_weights, ogp, egp, matrix_offset):
    """Local alignment with given parameters"""
    return aln.salign(rms_cutoff=3.5, normalize_pp_scores=False,
                      rr_file='$(LIB)/as1.sim.mat', overhang=0,
                      gap_penalties_1d=(ogp, egp),
                      local_alignment=True, matrix_offset=matrix_offset,
                      matrix_offset_3d=-0.5, gap_penalties_3d=(0, 3),
                      gap_gap_score=0, gap_residue_score=0,
                      alignment_type='tree', nsegm=2,
                      feature_weights=feature_weights,
                      improve_alignment=True, fit=True, write_fit=False,
                      output='ALIGNMENT QUALITY')

def _salign_fw_gaps3(aln, feature_weights, ogp3d, egp3d):
    """Global alignment with given parameters"""
    ogp = ogp3d
    egp = egp3d
    def _salign_wrap(aln, **keys):
        try:
            return aln.salign(auto_overhang=True, overhang_auto_limit=5,
                              overhang_factor=1, **keys)
        except modeller.ModellerError, detail:
            print "SALIGN with auto_overhang failed: %s" % str(detail)
            print "Retrying without auto_overhang"
            return aln.salign(**keys)

    return _salign_wrap(aln, rms_cutoff=3.5, normalize_pp_scores=False,
                        rr_file='$(LIB)/as1.sim.mat', overhang=0,
                        gap_penalties_1d=(ogp, egp),
                        local_alignment=False, matrix_offset=-0.2,
                        gap_penalties_3d=(ogp3d, egp3d), gap_gap_score=0,
                        gap_residue_score=0, alignment_type='tree',
                        nsegm=2, feature_weights=feature_weights,
                        improve_alignment=True, fit=True, write_fit=False,
                        write_whole_pdb=False, output='ALIGNMENT QUALITY')

def _frange(start, end=None, inc=1.0):
    """A range function that accepts floating point increments"""
    if end is None:
        end = float(start)
        start = 0.0
    else:
        start = float(start)
    count = int((end - start)/inc)
    if start + (count*inc) != end:
        count += 1
    for i in range(count):
        yield start + i*inc

class _TemporaryDirectory(object):
    """Create a temporary directory, and delete it when this object
       goes out of scope."""
    def __init__(self):
        import tempfile
        import shutil
        self.shutil = shutil
        self.tmpdir = tempfile.mkdtemp()
    def __del__(self):
        if hasattr(self, 'shutil'):
            self.shutil.rmtree(self.tmpdir)
    def get_path(self, path):
        """Return the name of a file in the temporary directory"""
        import os
        return os.path.join(self.tmpdir, path)

def iterative_structural_align(aln):
    """Given an alignment of structures, iterate over parameter values
       to obtain the best structural alignment."""
    for seq in aln:
        if not hasattr(seq, 'atoms'):
            raise modeller.ModellerError("This method only works for an " + \
                                         "alignment of structures.")
    tmpdir = _TemporaryDirectory()
    fil = tmpdir.get_path("inp.pir")
    aln.write(file=fil)

    opfile = tmpdir.get_path("salign_local_mid.ali")
    opfile2 = tmpdir.get_path("salign_local.ali")

    # -- Iterating over values of gap penalties and matrix offset
    qmax = 0.0
    win_ogp3d = None
    fw1=(1., 0., 0., 0., 1., 0.)
    fw2=(0., 1., 0., 0., 0., 0.)
    fw3=(0., 0., 0., 0., 1., 0.)

    # -- Iterating over gap penalties 1D to get initial alignments
    print "Iterate over 1D penalties to get initial alignments"
    for ogp in _frange(-150, 1, 50):
        for egp in _frange(-50, 1, 50):
            for mo in _frange(-3.0, -0.05, 0.3):
                aln.clear()
                aln.append(file=fil)
                try:
                    qwlty1 = _salign_fw_local_gaps1(aln,fw1,ogp,egp,mo)
                    if qwlty1.qscorepct >= qmax:
                        qmax = qwlty1.qscorepct
                        aln.write(file=opfile, alignment_format='PIR')
                        win_ogp = ogp
                        win_egp = egp
                        win_mo = mo
                    print "Qlty scrs", ogp,"\t",egp,"\t",qwlty1.qscorepct
                except modeller.ModellerError, detail:
                    print "Set of parameters",fw1,ogp,egp,"resulted in the following error\t"+str(detail)


    # -- Iterating over gap penalties 3D to get final alignments
    print "Iterate over 3D penalties to get final alignments"
    qmax3d = 0.
    for ogp3d in _frange(0, 3, 1):
        for egp3d in range (2, 5, 1):
            aln.clear()
            aln.append(file=opfile)
            try:
                qwlty2 = _salign_fw_gaps3(aln,fw2,ogp3d,egp3d)
                if qwlty2.qscorepct >= qmax3d:
                    qmax3d = qwlty2.qscorepct
                    aln.write(file=opfile2, alignment_format='PIR')
                    win_ogp3d = ogp3d
                    win_egp3d = egp3d
                print "Qlty scrs", ogp3d,"\t",egp3d,"\t",qwlty2.qscorepct
            except modeller.ModellerError,detail:
                print "Set of parameters",fw2,ogp3d,egp3d,"resulted in the following error\t"+str(detail)
    qmax = max(qmax, qmax3d)

    # try alternate initial alignments only if the qmax score is less than 70%
    qmax_old = qmax
    if qmax_old <= 70:
        print "Trying alternate initial alignments"
        for ogp in _frange(0.0, 2.2, 0.3):
            for egp in _frange(0.1, 2.3, 0.3):
                for mo in _frange(-3.0, -0.05, 0.3):
                    aln.clear()
                    aln.append(file=fil)
                    try:
                        qwlty1 = _salign_fw_local_gaps1(aln,fw3,ogp,egp,mo)
                        if qwlty1.qscorepct >= qmax:
                            qmax = qwlty1.qscorepct
                            aln.write(file=opfile, alignment_format='PIR')
                            win_ogp = ogp
                            win_egp = egp
                            win_mo = mo
                        print "Qlty scrs", ogp,"\t",egp,"\t",qwlty1.qscorepct
                    except modeller.ModellerError, detail:
                        print "Set of parameters",fw3,ogp,egp,"resulted in the following error\t"+str(detail)

        # -- Iterating over gap penalties 3D to get final alignments
        print "Trying alternate final alignments"
        qmax3d = 0.
        for ogp3d in _frange(0, 3, 1):
            for egp3d in range(2, 5, 1):

                aln.clear()
                aln.append(file=opfile)
                try:
                    qwlty2 = _salign_fw_gaps3(aln,fw2,ogp3d,egp3d)
                    if qwlty2.qscorepct >= qmax3d:
                        qmax3d = qwlty2.qscorepct
                        aln.write(file=opfile2, alignment_format='PIR')
                        win_ogp3d = ogp3d
                        win_egp3d = egp3d
                    print "Qlty scrs", ogp3d,"\t",egp3d,"\t",qwlty2.qscorepct
                except modeller.ModellerError,detail:
                    print "Set of parameters",fw2,ogp3d,egp3d,"resulted in the following error\t"+str(detail)
        qmax = max(qmax, qmax3d)

    print "final max quality = ", qmax

    if win_ogp3d is None:
        raise modeller.ModellerError("Structure alignment failed")
    else:
        aln.clear()
        aln.append(file=opfile2)
