import modeller
import _modeller
import model
import profile
import util.modutil as modutil
import util.modlist as modlist
from modeller.util.modobject import modobject
from modeller.util.logger import log
from modeller import sequence, modfile, alnsequence, alnstructure

__docformat__ = "epytext en"

class alignment(modobject):
    """Holds an alignment of protein sequences"""
    __modpt = None
    __read_one_pair = None
    env = None

    def __init__(self, env, **vars):
        self.__modpt = _modeller.mod_alignment_new(self)
        self.env = env.copy()
        if len(vars) > 0:
            self.append(**vars)

    def __repr__(self):
        if len(self) == 0:
            return "Empty alignment"
        else:
            return "Alignment of " + ", ".join([repr(s) for s in self])

    def __str__(self):
        if len(self) == 0:
            return "<Empty alignment>"
        else:
            return "<Alignment of " + ", ".join([str(s) for s in self]) + ">"

    def __setstate__(self, d):
        self.__dict__.update(d)
        self.__modpt = _modeller.mod_alignment_new(self)

    def __del__(self):
        if self.__modpt:
            _modeller.mod_alignment_free(self.__modpt)

    def __get_modpt(self):
        return self.__modpt

    def __len__(self):
        return _modeller.mod_alignment_nseq_get(self.modpt)

    def clear(self):
        """Remove all sequences from the alignment"""
        return _modeller.mod_alignment_clear(self.modpt)

    def append(self, file, align_codes='all', atom_files=None, remove_gaps=True,
               alignment_format='PIR', io=None, allow_alternates=False):
        """Add sequence(s) from an alignment file"""
        if io is None:
            io = self.env.io
        new_align_codes = align_codes
        new_atom_files = atom_files
        align_codes = [seq.code for seq in self]
        atom_files = [seq.atom_file for seq in self]
        if new_align_codes is not None:
            if hasattr(new_align_codes, '__iter__'):
                align_codes.extend(new_align_codes)
            else:
                align_codes.append(new_align_codes)
        if new_atom_files is not None:
            fmt = alignment_format.upper()
            if fmt != 'PAP' and fmt != 'INSIGHT':
                raise TypeError("atom_files is only used for PAP and " + \
                                "INSIGHT alignment files")
            if hasattr(new_atom_files, '__iter__'):
                atom_files.extend(new_atom_files)
            else:
                atom_files.append(new_atom_files)
        func = _modeller.mod_alignment_read2
        return func(self.modpt, io.modpt, self.env.libs.modpt, align_codes,
                    atom_files, file, remove_gaps, alignment_format,
                    allow_alternates)

    def read(self, **vars):
        """Read sequence(s) from an alignment file"""
        self.clear()
        return self.append(**vars)

    def read_one(self, file, remove_gaps=False, alignment_format='PIR',
                 io=None, allow_alternates=False):
        """Read sequences from the file, one by one. Return False when no more
           can be read."""
        if io is None:
            io = self.env.io
        self.clear()
        func = _modeller.mod_alignment_read_one2
        if isinstance(file, str):
            # Cache one file handle per alignment:
            if self.__read_one_pair is None or self.__read_one_pair[0] != file:
                log.warning("read_one",
                            "Filename arguments to read_one() are " + \
                            "deprecated.\nPlease use a handle " + \
                            "(modfile.File object) instead.")
                self.__read_one_pair = (file, modfile.File(file, 'r'))
            file = self.__read_one_pair[1]
        return func(self.modpt, file.file_pointer, io.modpt,
                    self.env.libs.modpt, remove_gaps, alignment_format,
                    allow_alternates) == 0

    def append_model(self, mdl, align_codes, atom_files=''):
        """Add a sequence from a model"""
        return _modeller.mod_alignment_append_model(self.modpt, mdl.modpt,
                                                    self.env.libs.modpt,
                                                    atom_files, align_codes)

    def append_sequence(self, sequence):
        """Add a new sequence, as a string of one-letter codes"""
        return _modeller.mod_alignment_append_sequence(self.modpt, sequence,
                                                       self.env.libs.modpt)

    def compare_with(self, aln):
        """Compare with another alignment"""
        return _modeller.mod_alignment_compare_with(aln=self.modpt,
                                                    aln2=aln.modpt)

    def compare_structures(self, compare_mode=3, fit=True, fit_atoms='CA',
                           matrix_file='family.mat', output='LONG',
                           asgl_output=False, refine_local=True,
                           rms_cutoffs=(3.5, 3.5, 60., 60., 15., 60., 60.,
                                        60., 60., 60., 60.),
                           varatom='CA', edat=None, io=None):
        """Compare 3D structures"""
        if io is None:
            io = self.env.io
        if edat is None:
            edat = self.env.edat
        func = _modeller.mod_compare_structures
        return func(self.modpt, edat.modpt, io.modpt, self.env.libs.modpt,
                    compare_mode, fit, fit_atoms, matrix_file, output,
                    asgl_output, refine_local, rms_cutoffs, varatom)

    def compare_sequences(self, mdl, matrix_file, variability_file,
                          max_gaps_match, rr_file='${LIB}/as1.sim.mat'):
        """Compare sequences in alignment"""
        return _modeller.mod_compare_sequences(self.modpt, mdl.modpt,
                                               self.env.libs.modpt, rr_file,
                                               matrix_file, variability_file,
                                               max_gaps_match)

    def segment_matching(self, file, root_name, file_ext, file_id, align_block,
                         segment_report, segment_cutoff, segment_shifts,
                         segment_growth_n, segment_growth_c, min_loop_length,
                         rr_file='${LIB}/as1.sim.mat'):
        """Enumerates alignments between two blocks of sequences.
           More precisely, it enumerates the alignments between the segments in
           the first block, as defined by C{align_block}, and the sequences in
           the second block. The segments can be moved to the left and right as
           well as lengthened and shortened, relative to the initial alignment.
        """
        return _modeller.mod_segment_matching(self.modpt, self.env.libs.modpt,
                                              rr_file, file, align_block,
                                              min_loop_length, segment_shifts,
                                              segment_growth_n,
                                              segment_growth_c,
                                              segment_cutoff, segment_report,
                                              root_name, file_id, file_ext)

    def edit(self, overhang, edit_align_codes, base_align_codes,
             min_base_entries, io=None):
        """Edit overhangs in alignment.
           @param edit_align_codes: specifies the alignment codes for the
                  alignment entries whose overhangs are to be cut; in
                  addition, C{all} or C{last} can be used.
           @param base_align_codes: specifies the alignment codes for the
                  alignment entries that are used to determine the extent of
                  the overhangs to be cut from the edited entries; in
                  addition, C{all} or C{rest} (relative to C{edit_align_codes})
                  can be used.
        """
        if io is None:
            io = self.env.io
        return _modeller.mod_alignment_edit(self.modpt, io.modpt,
                                            self.env.libs.modpt, overhang,
                                            edit_align_codes, base_align_codes,
                                            min_base_entries)

    def append_profile(self, prf):
        """Add sequences from a profile"""
        return _modeller.mod_profile_to_aln(prf=prf.modpt, aln=self.modpt,
                                            libs=self.env.libs.modpt)

    def to_profile(self):
        """Converts the alignment to profile format"""
        prf = profile.profile(self.env)
        _modeller.mod_profile_from_aln(aln=self.modpt, prf=prf.modpt,
                                       libs=self.env.libs.modpt)
        return prf

    def check(self, io=None):
        """Check alignment for modeling"""
        self.check_structure_structure(io=io)
        self.check_sequence_structure(io=io)

    def check_structure_structure(self, eqvdst=6.0, io=None):
        """Check pairwise structural superpositions of all template
           structures."""
        if io is None:
            io = self.env.io
        f = _modeller.mod_alignment_check_structures
        return f(self.modpt, io.modpt, self.env.libs.modpt, eqvdst)

    def check_sequence_structure(self, gapdist=8.0, io=None):
        """Check the current sequence/structure alignment for sanity."""
        if io is None:
            io = self.env.io
        f = _modeller.mod_alignment_check_seqstruc
        return f(self.modpt, io.modpt, self.env.libs.modpt, gapdist)

    def describe(self, io=None):
        """Describe proteins"""
        if io is None:
            io = self.env.io
        return _modeller.mod_alignment_describe(aln=self.modpt, io=io.modpt,
                                                libs=self.env.libs.modpt)

    def consensus(self, align_block=0, gap_penalties_1d=(-900., -50.),
                  weigh_sequences=False, input_weights_file=None,
                  output_weights_file=None, weights_type='SIMILAR',
                  smooth_prof_weight=10):
        """Produce a consensus alignment"""
        (input_weights_file, read_weights) = self.__opt_file(input_weights_file)
        (output_weights_file,
         write_weights) = self.__opt_file(output_weights_file)
        func = _modeller.mod_alignment_consensus
        return func(self.modpt, self.env.libs.modpt, align_block,
                    gap_penalties_1d, read_weights, write_weights,
                    weigh_sequences, input_weights_file, output_weights_file,
                    weights_type, smooth_prof_weight)

    def __opt_file(self, fname):
        """Processing for optional files for some routines"""
        if fname is None or fname == "":
            return ('', False)
        else:
            return (fname, True)

    def align(self, off_diagonal=100, local_alignment=False, matrix_offset=0.,
              gap_penalties_1d=(-900., -50.), n_subopt=1, subopt_offset=0.,
              weigh_sequences=False, smooth_prof_weight=10, align_what='BLOCK',
              weights_type='SIMILAR', input_weights_file=None,
              output_weights_file=None, rr_file='$(LIB)/as1.sim.mat',
              overhang=0, align_block=0):
        """Align two (blocks of) sequences"""
        (input_weights_file, read_weights) = self.__opt_file(input_weights_file)
        (output_weights_file,
         write_weights) = self.__opt_file(output_weights_file)
        func = _modeller.mod_align
        return func(self.modpt, self.env.libs.modpt, off_diagonal,
                    local_alignment, matrix_offset, gap_penalties_1d,
                    read_weights, write_weights, n_subopt, subopt_offset,
                    weigh_sequences, smooth_prof_weight, align_what,
                    weights_type, input_weights_file, output_weights_file,
                    rr_file, overhang, align_block)

    def malign(self, rr_file='$(LIB)/as1.sim.mat', off_diagonal=100,
               local_alignment=False, matrix_offset=0., overhang=0,
               align_block=0, gap_penalties_1d=(-900., -50.)):
        """Align two or more sequences"""
        return _modeller.mod_malign(self.modpt, self.env.libs.modpt, rr_file,
                                    off_diagonal, local_alignment,
                                    matrix_offset, overhang, align_block,
                                    gap_penalties_1d)

    def align2d(self, overhang=0, align_block=0, rr_file='$(LIB)/as1.sim.mat',
                align_what='BLOCK', off_diagonal=100, max_gap_length=999999,
                local_alignment=False, matrix_offset=0.,
                gap_penalties_1d=(-100., 0.),
                gap_penalties_2d=(3.5, 3.5, 3.5, 0.2, 4.0, 6.5, 2.0, 0., 0.),
                surftyp=1, fit=True, fix_offsets=(0., -1., -2., -3., -4.),
                input_weights_file=None, output_weights_file=None, n_subopt=1,
                subopt_offset=0., input_profile_file=None,
                output_profile_file=None, weigh_sequences=False,
                smooth_prof_weight=10, weights_type='SIMILAR', io=None):
        """Align sequences with structures"""
        if io is None:
            io = self.env.io
        (input_weights_file, read_weights) = self.__opt_file(input_weights_file)
        (output_weights_file,
         write_weights) = self.__opt_file(output_weights_file)
        (input_profile_file, read_profile) = self.__opt_file(input_profile_file)
        (output_profile_file,
         write_profile) = self.__opt_file(output_profile_file)

        func = _modeller.mod_align2d
        return func(self.modpt, io.modpt, self.env.libs.modpt, overhang,
                    align_block, rr_file, align_what, off_diagonal,
                    max_gap_length, local_alignment, matrix_offset,
                    gap_penalties_1d, gap_penalties_2d, surftyp, fit,
                    fix_offsets, read_weights, write_weights,
                    input_weights_file, output_weights_file, n_subopt,
                    subopt_offset, read_profile, input_profile_file,
                    write_profile, output_profile_file, weigh_sequences,
                    smooth_prof_weight, weights_type)

    def align3d(self, off_diagonal=100, overhang=0, local_alignment=False,
                matrix_offset=0., gap_penalties_3d=(0.0, 1.75), fit=True,
                fit_atoms='CA', align3d_trf=False, output='LONG',
                align3d_repeat=False, io=None):
        """Align two structures"""
        if io is None:
            io = self.env.io
        func = _modeller.mod_align3d
        return func(self.modpt, io.modpt, self.env.libs.modpt, off_diagonal,
                    overhang, local_alignment, matrix_offset, gap_penalties_3d,
                    fit, fit_atoms, align3d_trf, output, align3d_repeat)

    def malign3d(self, off_diagonal=100, overhang=0, local_alignment=False,
                 matrix_offset=0., gap_penalties_3d=(0.0, 1.75), fit=True,
                 fit_atoms='CA', output='LONG', write_whole_pdb=True,
                 current_directory=True, write_fit=False,
                 edit_file_ext=('.pdb', '_fit.pdb'), io=None):
        """Align two or more structures"""
        if io is None:
            io = self.env.io
        return _modeller.mod_malign3d(self.modpt, io.modpt, self.env.libs.modpt,
                                      off_diagonal, overhang, local_alignment,
                                      matrix_offset, gap_penalties_3d, fit,
                                      fit_atoms, output, write_whole_pdb,
                                      current_directory, write_fit,
                                      edit_file_ext)

    def salign(self, residue_type2='REGULAR', no_ter=False, overhang=0,
               off_diagonal=100, matrix_offset=0.,
               gap_penalties_1d=(-900., -50.),
               gap_penalties_2d=(3.5, 3.5, 3.5, 0.2, 4.0, 6.5, 2.0, 0., 0.),
               gap_penalties_3d=(0.0, 1.75),
               feature_weights=(1., 0., 0., 0., 0., 0.), rms_cutoff=3.5,
               fit=True, surftyp=1, fit_on_first=False, gap_function=False,
               align_block=0, max_gap_length=999999, align_what='BLOCK',
               input_weights_file=None, output_weights_file=None,
               weigh_sequences=False, smooth_prof_weight=10,
               fix_offsets=(0., -1., -2., -3., -4.), substitution=False,
               comparison_type='MAT', matrix_comparison='CC',
               alignment_type='PROGRESSIVE', edit_file_ext=('.pdb', '_fit.pdb'),
               weights_type='SIMILAR', similarity_flag=False,
               bkgrnd_prblty_file='$(LIB)/blosum62_bkgrnd.prob',
               ext_tree_file=None, dendrogram_file='',
               matrix_scaling_factor=0.0069, auto_overhang=False,
               overhang_factor=0.4, overhang_auto_limit=60,
               local_alignment=False, improve_alignment=True, fit_atoms='CA',
               output='', write_whole_pdb=True, current_directory=True,
               write_fit=False, fit_pdbnam=True, rr_file='$(LIB)/as1.sim.mat',
               n_subopt=1, subopt_offset=0., align3d_trf=False,
               normalize_pp_scores=False, gap_gap_score=0.,
               gap_residue_score=0., nsegm=2, matrix_offset_3d=-0.1, io=None):
        """Align two or more proteins"""
        import salign
        if io is None:
            io = self.env.io
        (input_weights_file, read_weights) = self.__opt_file(input_weights_file)
        (output_weights_file,
         write_weights) = self.__opt_file(output_weights_file)
        if ext_tree_file is None:
            ext_tree_file = ''
        func = _modeller.mod_salign
        retval = func(self.modpt, io.modpt, self.env.libs.modpt, residue_type2,
                      no_ter, overhang, off_diagonal, matrix_offset,
                      gap_penalties_1d, gap_penalties_2d, gap_penalties_3d,
                      feature_weights, rms_cutoff, fit, surftyp,
                      fit_on_first, gap_function, align_block,
                      max_gap_length, align_what, read_weights,
                      write_weights, input_weights_file,
                      output_weights_file, weigh_sequences,
                      smooth_prof_weight, fix_offsets, substitution,
                      comparison_type, matrix_comparison,
                      alignment_type, edit_file_ext, weights_type,
                      similarity_flag, bkgrnd_prblty_file,
                      ext_tree_file, dendrogram_file, matrix_scaling_factor,
                      auto_overhang, overhang_factor,
                      overhang_auto_limit, local_alignment,
                      improve_alignment, fit_atoms, output,
                      write_whole_pdb, current_directory,
                      write_fit, fit_pdbnam, rr_file, n_subopt,
                      subopt_offset, align3d_trf, normalize_pp_scores,
                      gap_gap_score, gap_residue_score,
                      nsegm, matrix_offset_3d)
        return salign.SalignData(*retval)

    def write(self, file, alignment_format='PIR',
              alignment_features='INDICES CONSERVATION', align_block=0,
              align_alignment=False):
        """Write the alignment to a file"""
        if isinstance(file, str):
            file = modfile.File(file, 'w')
        return _modeller.mod_alignment_write(self.modpt, self.env.libs.modpt,
                                             file.file_pointer,
                                             alignment_format,
                                             alignment_features, align_block,
                                             align_alignment)

    def id_table(self, matrix_file):
        """Calculate percentage sequence identities"""
        from modeller.id_table import id_table
        return id_table(self, matrix_file)

    def get_suboptimals(self, f):
        """Read the suboptimal alignments from the Python file object f
           into the current alignment, one by one. Returns the alignment object
           itself each time."""
        if len(self) != 2:
            raise TypeError("alignment should contain exactly 2 sequences")
        while True:
            a = f.readline()
            if not a: break
            if not a.startswith('ALIGNMENT:'):
                raise modeller.FileFormatError("expecting ALIGNMENT: line " + \
                                         "in suboptimal alignment output file")
            seq1 = f.readline().split()
            pref = seq1.pop(0)
            if pref != 'SEQ1:':
                raise modeller.FileFormatError("expecting SEQ1: line " + \
                                         "in suboptimal alignment output file")
            seq1 = [int(x) for x in seq1]
            seq2 = f.readline().split()
            pref = seq2.pop(0)
            if pref != 'SEQ2:':
                raise modeller.FileFormatError("expecting SEQ2: line " + \
                                         "in suboptimal alignment output file")
            seq2 = [int(x) for x in seq2]

            _set_seq_gaps_subopt(self[0], seq1)
            _set_seq_gaps_subopt(self[1], seq2)
            if seq1[0] > 1 and seq2[0] > 1:
                raise NotImplementedError("Neither sequence starts at 1")
            elif seq1[0] > 1:
                self[1].residues[0].add_leading_gaps(seq1[0] - seq2[0])
            elif seq2[0] > 1:
                self[0].residues[0].add_leading_gaps(seq2[0] - seq1[0])
            yield self

    def __contains__(self, code):
        return _modeller.mod_alignment_find_code(self.modpt, code) >= 0

    def keys(self):
        return [seq.code for seq in self]

    def __delitem__(self, indx):
        ret = modutil.handle_seq_indx(self, indx,
                                      _modeller.mod_alignment_find_code,
                                      (self.modpt,))
        _modeller.mod_alnsequence_del(self.modpt, ret)

    def __getitem__(self, indx):
        ret = modutil.handle_seq_indx(self, indx,
                                      _modeller.mod_alignment_find_code,
                                      (self.modpt,))
        if isinstance(ret, int):
            if _modeller.mod_alnsequence_has_structure(self.modpt, ret):
                return alnstructure.Structure(self, ret)
            else:
                return alnsequence.Sequence(self, ret)
        else:
            return [self[ind] for ind in ret]

    def __get_comments(self):
        return modlist.SimpleVarList(self.modpt,
                                     _modeller.mod_alignment_ncomment_get,
                                     _modeller.mod_alignment_ncomment_set,
                                     _modeller.mod_alignment_comment_get,
                                     _modeller.mod_alignment_comment_set)
    def __set_comments(self, obj):
        modlist.set_varlist(self.comments, obj)
    def __del_comments(self):
        modlist.del_varlist(self.comments)
    def __get_positions(self):
        return PositionList(self)

    modpt = property(__get_modpt)
    comments = property(__get_comments, __set_comments, __del_comments,
                        doc="Alignment file comments")
    positions = property(__get_positions, doc="Alignment positions")


class Position(object):
    """An alignment position"""

    def __init__(self, aln, indx):
        self.__aln = aln
        self.__indx = indx

    def get_residue(self, seq):
        """Get the residue in C{seq} that is at this alignment position, or None
           if a gap is present."""
        aln = self.__aln
        if not isinstance(seq, alnsequence.Sequence):
            raise TypeError("Expected an alignment 'Sequence' object for seq")
        if seq.aln != aln:
            raise ValueError("seq must be a sequence in the same alignment")
        ialn = _modeller.mod_alignment_ialn_get(aln.modpt)
        ires = _modeller.mod_int2_get(ialn, self.__indx, seq._num)
        if ires == 0:
            return None
        else:
            return seq.residues[ires-1]

    def __get_num(self):
        return self.__indx
    def __get_prof(self, typ):
        prof = _modeller.mod_alignment_prof_get(self.__aln.modpt)
        return _modeller.mod_float2_get(prof, self.__indx, typ)
    def __get_helix(self):
        return self.__get_prof(0)
    def __get_strand(self):
        return self.__get_prof(1)
    def __get_buried(self):
        return self.__get_prof(2)
    def __get_straight(self):
        return self.__get_prof(3)
    num = property(__get_num)
    helix = property(__get_helix, doc="Helix secondary structure")
    strand = property(__get_strand, doc="Strand secondary structure")
    buried = property(__get_buried, doc="Buriedness")
    straight = property(__get_straight, doc="Straightness")


class PositionList(modlist.FixList):
    """A list of L{Position} objects."""

    def __init__(self, aln):
        self.__aln = aln
        modlist.FixList.__init__(self)

    def __len__(self):
        return _modeller.mod_alignment_naln_get(self.__aln.modpt)

    def _getfunc(self, indx):
        return Position(self.__aln, indx)


def _set_seq_gaps_subopt(seq, subopt_pos):
    "Set gaps in a sequence to match those implied by a suboptimal alignment"
    pos_n = 0
    for (n, res) in enumerate(seq.residues):
        gaps = res.get_leading_gaps()
        if pos_n >= len(subopt_pos):
            res.remove_leading_gaps(gaps)
        else:
            seq1_pos = subopt_pos[pos_n] - 1
            if n < seq1_pos:
                res.remove_leading_gaps(gaps)
            else:
                gaps_to_add = 0
                while pos_n < len(subopt_pos) and subopt_pos[pos_n] == 0:
                    pos_n += 1
                    gaps_to_add += 1
                if gaps_to_add > gaps:
                    res.add_leading_gaps(gaps_to_add - gaps)
                elif gaps_to_add < gaps:
                    res.remove_leading_gaps(gaps - gaps_to_add)
                if subopt_pos[pos_n] != n + 1:
                    raise Modeller.FileFormatError("suboptimal alignment " + \
                              "file SEQ1 and SEQ2 indices should be sequential")
                pos_n += 1
