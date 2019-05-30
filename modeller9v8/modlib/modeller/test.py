import unittest
import modeller
import math
from cStringIO import StringIO
import sys

__all__ = ['ModellerTest']

class ModellerTest(unittest.TestCase):
    _env = None

    def get_environ(self):
        """Set up the Modeller environment, and keep a reference to it so that
           we only need to do it once for all tests"""
        if not self._env:
            env = modeller.environ()
            env.io.atom_files_directory = ['../data', 'test/data']
            env.edat.dynamic_sphere = False
            env.libs.topology.read(file='${LIB}/top_heav.lib')
            env.libs.parameters.read(file='${LIB}/par.lib')
            ModellerTest._env = env
        return self._env

    def assertModelsEqual(self, mdl1, mdl2):
        """Check to make sure that models |mdl1| and |mdl2| are the same"""
        self.assertEqual(len(mdl1.atoms), len(mdl2.atoms))
        for (at1, at2) in zip(mdl1.atoms, mdl2.atoms):
            self.assertAtomsEqual(at1, at2)

    def assertAtomsEqual(self, at1, at2):
        """Check to make sure atoms |at1| and |at2| are the same"""
        self.assertEqual(at1.name, at2.name)
        self.assertAlmostEqual(at1.x, at2.x, places=3)

    def assertDistance(self, o1, o2, mean, tolerance):
        """Check to make sure the distance between o1 and o2 is as expected"""
        distsq = (o1.x-o2.x)**2 + (o1.y-o2.y)**2 + (o1.z-o2.z)**2
        dist = distsq ** 0.5
        msg = "Distance between %s and %s is %f - expected %f" % (o1, o2,
                                                                  dist, mean)
        self.assertInTolerance(dist, mean, tolerance, msg)

    def assertAlignmentsEqual(self, refaln, aln, check_meta=True):
        """Check to make sure that alignments |refaln| and |aln| are the same"""
        self.assertEqual(len(refaln), len(aln),
                         "Inconsistent number of sequences (%d vs. %d)" \
                         % (len(refaln), len(aln)))
        for (seq1, seq2) in zip(refaln, aln):
            self.assertAlnSequencesEqual(seq1, seq2, check_meta)
        self.assertEqual(len(refaln.comments), len(aln.comments))
        for (c1, c2) in zip(refaln.comments, aln.comments):
            self.assertEqual(c1, c2)

    def assertAlnSequencesEqual(self, seq1, seq2, check_meta=True):
        """Check to make sure that sequences |seq1| and |seq2| are the same"""
        self.assertEqual(len(seq1), len(seq2),
                         "Inconsistent number of residues (%d vs. %d)" \
                         % (len(seq1), len(seq2)))
        if check_meta:
            self.assertEqual(seq1.range[0], seq2.range[0])
            self.assertEqual(seq1.range[1], seq2.range[1])
            self.assertEqual(seq1.resolution, seq2.resolution)
            self.assertEqual(seq1.rfactor, seq2.rfactor)
            self.assertEqual(seq1.code, seq2.code)
            self.assertEqual(seq1.atom_file, seq2.atom_file)
            self.assertEqual(seq1.name, seq2.name)
            self.assertEqual(seq1.prottyp, seq2.prottyp)
            self.assertEqual(seq1.source, seq2.source)
        for (res1, res2) in zip(seq1.residues, seq2.residues):
            self.assertResiduesEqual(res1, res2)
            self.assertEqual(res1.get_position().num, res2.get_position().num)
            self.assertEqual(res1.get_leading_gaps(), res2.get_leading_gaps())
            self.assertEqual(res1.get_trailing_gaps(), res2.get_trailing_gaps())

    def assertResiduesEqual(self, res1, res2):
        """Check to make sure that residues |res1| and |res2| are the same"""
        self.assertEqual(res1.name, res2.name,
                         "Inconsistent residue names (%s vs. %s)" \
                         % (res1.name, res2.name))

    def assertProfilesEqual(self, prf1, prf2):
        """Check to make sure that two profile objects are the same"""
        self.assertEqual(len(prf1), len(prf2))
        for (seq1, seq2) in zip(prf1, prf2):
            self.assertProfileSeqsEqual(seq1, seq2)

    def assertProfileSeqsEqual(self, seq1, seq2):
        """Check to make sure that two profile sequences are the same"""
        self.assertEqual(len(seq1), len(seq2))
        self.assertEqual(seq1.code, seq2.code)
        self.assertEqual(seq1.prottyp, seq2.prottyp)
        self.assertEqual(seq1.iter, seq2.iter)
        self.assertEqual(seq1.neqv, seq2.neqv)
        self.assertEqual(seq1.fid, seq2.fid)
        self.assertEqual(seq1.evalue, seq2.evalue)

    def assertSeqDBsEqual(self, sdb1, sdb2):
        """Check to make sure that two sequence_db objects are the same"""
        self.assertEqual(len(sdb1), len(sdb2))
        for (seq1, seq2) in zip(sdb1, sdb2):
            self.assertSeqDBseqsEqual(seq1, seq2)

    def assertSeqDBseqsEqual(self, seq1, seq2):
        """Check to make sure that two sequence_db sequences are the same"""
        self.assertEqual(len(seq1), len(seq2))
        self.assertEqual(seq1.code, seq2.code)
        self.assertEqual(seq1.prottyp, seq2.prottyp)
        self.assertEqual(seq1.resol, seq2.resol)
        self.assertEqual(len(seq1.residues), len(seq2.residues))
        for (r1, r2) in zip(seq1.residues, seq2.residues):
            self.assertSeqDBresEqual(r1, r2)

    def assertSeqDBresEqual(self, res1, res2):
        """Check to make sure that two sequence_db residues are the same"""
        self.assertEqual(res1.type, res2.type)

    def assertInTolerance(self, num1, num2, tolerance, msg=None):
        """Assert that the difference between num1 and num2 is less than
           tolerance"""
        diff = abs(num1 - num2)
        if msg is None:
            msg = "%f != %f within %g" % (num1, num2, tolerance)
        self.assert_(diff < tolerance, msg)

    def run_capture_stdout(self, method, *args, **keys):
        """Run a method and capture its standard output. Returns both the
           method's own return values and the standard output."""
        saved_stdout = sys.stdout
        sio = StringIO()
        try:
            sys.stdout = sio
            ret = method(*args, **keys)
        finally:
            sys.stdout = saved_stdout
        val = sio.getvalue()
        sys.stdout.write(val)
        return ret, val

    def create_angle(self, angle, a1, a2, a3):
        """Set the coordinates of the given atoms to form (by construction)
           the given angle between them."""
        a1.x, a1.y, a1.z = -10.0, 0.0, 0.0
        a2.x, a2.y, a2.z = 0.0, 0.0, 0.0
        a3.x, a3.y, a3.z = -math.cos(angle), math.sin(angle), 0.0

    def create_dihedral(self, angle, a1, a2, a3, a4):
        """Set the coordinates of the given atoms to form (by construction)
           the given dihedral angle between them."""
        a1.x, a1.y, a1.z = 10.0, 0.0, -10.0
        a2.x, a2.y, a2.z = 0.0, 0.0, -10.0
        a3.x, a3.y, a3.z = 0.0, 0.0, 0.0
        a4.x, a4.y, a4.z = math.cos(angle), math.sin(angle), 0.0
