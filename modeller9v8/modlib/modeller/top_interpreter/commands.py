"""Classes for running legacy TOP commands"""

import sys
from modeller import modfile
from modeller.util.logger import log
from modeller.top_interpreter.topcmds import AutoCommands
import _modeller

__docformat__ = "epytext en"

class Commands(AutoCommands):
    """Class to run legacy TOP commands; augments the auto-generated code
       in the L{AutoCommands} class in many cases"""

    def __init__(self):
        if hasattr(sys, 'dllhandle'):
            # Only true on Windows systems
            self.__dirsep = '\\'
        else:
            self.__dirsep = '/'

    def __default_file(self, varname, deffile_id, exts=None):
        """Generate a 'default' filename using TOP rules"""
        file = self.vars[varname]
        if file.upper().find('DEFAULT') >= 0:
            root_name = self.vars['root_name']
            file_id = self.vars['file_id']
            id1 = self.vars['id1']
            id2 = self.vars['id2']
            file_ext = self.vars['file_ext']
            if file_id.upper().find('DEFAULT') >= 0:
                file_id = deffile_id
            file = modfile.default(root_name, file_id, id1, id2, file_ext)
            self.vars[varname] = file
        if exts is not None:
            dirs = self.vars['directory'].split(':')
            file = _modeller.mod_fullfn(file, dirs, exts)
            self.vars[varname] = file

    def __add_outdir(self, file):
        """Prepend C{output_directory} to C{file}."""
        file = self.vars[file]
        outdir = self.vars['output_directory']
        if len(outdir) == 0 or outdir.endswith(self.__dirsep) \
           or file.startswith(self.__dirsep):
            return outdir + file
        else:
            return outdir + self.__dirsep + file

    def __pickall(self):
        """Create the implicit 'all atom' selections for MODEL"""
        mdl = self.get_mdl(1)
        self.vars['sel1'] = _modeller.mod_selection_all(mdl)
        self.vars['sel2'] = _modeller.mod_selection_all(mdl)
        self.vars['sel3'] = _modeller.mod_selection_all(mdl)

    def __get_sched(self):
        """Get MODEL's schedule"""
        mdl = self.get_mdl(1)
        return _modeller.mod_model_sched_get(mdl)

    def __check_schedule(self):
        """Check the schedule step"""
        _modeller.mod_schedule_check(self.__get_sched())

    def __get_optimization_method(self):
        """Get the optimization method, from the user or schedule"""
        ometh = self.vars['optimization_method']
        if ometh == -999:
            ometh = _modeller.mod_schedule_optimizer_get(self.__get_sched())
        return ometh

    def __update_schedule(self):
        """Update the schedule with user-chosen residue span range"""
        return _modeller.mod_schedule_update(self.__get_sched(),
                                             self.vars['residue_span_range'])

    def __process_seg(self, res):
        """Look for a given code in the alignment"""
        if ":" not in res and len(res) > 0 and res[0] != '!':
            aln = self.get_aln()
            iseq = _modeller.mod_alignment_find_code(aln, res)
            if iseq < 0:
                log.error('getchnrng',
                          "There is no such CODE in the alignment: %s" % res)
            return (aln, iseq)
        return (None, None)

    def __seg_from_aln(self, model_segment):
        """Possibly update file and model_segment with alignment information"""

        file = self.vars['file']
        model_segment = self.vars[model_segment]
        # Get PDB name from alignment?
        (aln, iseq) = self.__process_seg(model_segment[1])
        if aln:
            alnseq = _modeller.mod_alignment_alnsequence_get(aln, iseq)
            file = _modeller.mod_alnsequence_atom_files_get(alnseq)
        # Get segment range from alignment?
        (aln, iseq) = self.__process_seg(model_segment[0])
        if aln:
            seq = _modeller.mod_alignment_sequence_get(aln, iseq)
            model_segment = [_modeller.mod_sequence_rng_get(seq, i) \
                             for i in (0,1)]

        return file, model_segment

    def __parse_residue_ids(self, mdl, resids):
        """Given string residue IDs, return integers"""
        cd = _modeller.mod_model_cd_get(mdl)
        seq = _modeller.mod_model_seq_get(mdl)
        num_resids = [_modeller.mod_coordinates_find_residue(cd, seq, r) \
                      for r in resids]
        for (n, id) in enumerate(num_resids):
            if id <= 0:
                log.error("iresind",
                          "Residue identifier not found: %s" % resids[n])
        return num_resids

    def __make_sstruc_restraints(self, modfunc):
        """Make alpha or strand restraints"""
        mdl = self.get_mdl()
        libs = self.get_libs()
        resids = self.__parse_residue_ids(mdl, self.vars['residue_ids'])
        modfunc(mdl, resids, libs)

    def __make_sheet_restraints(self):
        """Make sheet restraints"""
        mdl = self.get_mdl()
        cd = _modeller.mod_model_cd_get(mdl)
        seq = _modeller.mod_model_seq_get(mdl)
        atids = self.vars['atom_ids']
        func = _modeller.mod_coordinates_find_atom
        num_atids = [func(cd, seq, a) for a in atids]
        for (n, id) in enumerate(num_atids):
            if id <= 0:
                log.error("indxatm2",
                          "Atom identifier not found: %s" % atids[n])
        _modeller.mod_restraints_make_sheet(mdl, num_atids,
                                            self.vars['sheet_h_bonds'])

    def write_model(self):
        self.__default_file('file', '.B')
        fh = modfile.File(self.__add_outdir('file'), 'w')
        AutoCommands.write_model(self, fh=fh.file_pointer)
        fh.close()
    def write_model2(self):
        self.__default_file('file', '.B')
        AutoCommands.write_model2(self, file=self.__add_outdir('file'))
    def build_model(self):
        AutoCommands.build_model(self)
        self.__pickall()
    def read_restraints(self):
        if not self.vars['add_restraints']:
            _modeller.mod_restraints_clear(self.get_mdl())
        self.__default_file('file', '.C', ('', '.rsr'))
        AutoCommands.read_restraints(self)
    def write_restraints(self):
        self.__default_file('file', '.C')
        AutoCommands.write_restraints(self, file=self.__add_outdir('file'))
    def read_model(self):
        self.__default_file('file', 'B')
        (file, model_segment) = self.__seg_from_aln('model_segment')
        AutoCommands.read_model(self, file=file, model_segment=model_segment)
        self.__pickall()
    def read_model2(self):
        self.__default_file('file', '.B')
        (file, model2_segment) = self.__seg_from_aln('model2_segment')
        AutoCommands.read_model2(self, file=file, model2_segment=model2_segment)
    def segment_matching(self):
        AutoCommands.segment_matching(self,
                                      root_name=self.__add_outdir('root_name'))
    def switch_trace(self):
        self.__default_file('file', '.D')
        AutoCommands.switch_trace(self)
    def read_schedule(self):
        self.__default_file('file', 'S', ('', '.sch'))
        AutoCommands.read_schedule(self)
    def write_schedule(self):
        self.__default_file('file', '.S')
        AutoCommands.write_schedule(self)
    def id_table(self):
        self.__default_file('matrix_file', '.M')
        AutoCommands.id_table(self,
                              matrix_file=self.__add_outdir('matrix_file'))
    def write_topology_model(self):
        self.__default_file('file', '.toplib')
        AutoCommands.write_topology_model(self, file=self.__add_outdir('file'))
    def write_alignment(self):
        self.__default_file('file', '.A')
        fh = modfile.File(self.__add_outdir('file'), 'w')
        AutoCommands.write_alignment(self, fh=fh.file_pointer)
        fh.close()
    def sequence_comparison(self):
        self.__default_file('rr_file', '.A', ('', '.mat'))
        self.__default_file('matrix_file', '.M')
        self.__default_file('variability_file', '.V')
        AutoCommands.sequence_comparison(self,
            matrix_file=self.__add_outdir('matrix_file'),
            variability_file=self.__add_outdir('variability_file'))
    def sequence_to_ali(self):
        aln = self.get_aln()
        if not self.vars['add_sequence']:
            _modeller.mod_alignment_clear(aln)
        nseq = _modeller.mod_alignment_nseq_get(aln)
        try:
            align_codes = self.vars['align_codes'][nseq]
        except IndexError:
            align_codes = ''
        try:
            atom_files = self.vars['atom_files'][nseq]
        except IndexError:
            atom_files = ''
        AutoCommands.sequence_to_ali(self, align_codes=align_codes,
                                     atom_files=atom_files)
    def read_topology(self):
        if not self.vars['add_topology']:
            _modeller.mod_topology_clear(self.get_tpl(), self.get_libs())
        self.__default_file('file', '.T', ('', '.lib'))
        AutoCommands.read_topology(self)
    def generate_topology(self):
        seq = self.vars['sequence']
        iseq = _modeller.mod_alignment_find_code(self.get_aln(), seq)
        if iseq < 0:
            log.error('generate_topology',
                      "Sequence '%s' not found in alignment" % seq)
        if not self.vars['add_segment']:
            _modeller.mod_model_topology_clear(self.get_mdl())
        AutoCommands.generate_topology(self, iseq=iseq)
    def energy(self):
        self.__default_file('file', 'P')
        span = self.__update_schedule()
        AutoCommands.energy(self, residue_span_range=span,
                            file=self.__add_outdir('file'))
    def optimize(self):
        span = self.__update_schedule()
        self.__check_schedule()
        ometh = self.__get_optimization_method()
        AutoCommands.optimize(self, residue_span_range=span,
                              optimization_method=ometh)
    def debug_function(self):
        span = self.__update_schedule()
        AutoCommands.debug_function(self, residue_span_range=span)
    def spline_restraints(self):
        span = self.__update_schedule()
        AutoCommands.spline_restraints(self, residue_span_range=span)
    def pick_hot_atoms(self):
        span = self.__update_schedule()
        AutoCommands.pick_hot_atoms(self, residue_span_range=span)
    def make_schedule(self):
        span = self.__update_schedule()
        AutoCommands.make_schedule(self, residue_span_range=span)
    def make_restraints(self):
        if not self.vars['add_restraints']:
            _modeller.mod_restraints_clear(self.get_mdl())
        span = self.__update_schedule()
        rsrtype = self.vars['restraint_type'].lower()
        if rsrtype == 'alpha':
            self.__make_sstruc_restraints(_modeller.mod_restraints_make_alpha)
        elif rsrtype == 'strand':
            self.__make_sstruc_restraints(_modeller.mod_restraints_make_strand)
        elif rsrtype == 'sheet':
            self.__make_sheet_restraints()
        else:
            AutoCommands.make_restraints(self, residue_span_range=span,
                                         exclude_distance=0.0)
    def pick_restraints(self):
        if not self.vars['add_restraints']:
            mdl = self.get_mdl()
            rsr = _modeller.mod_model_rsr_get(mdl)
            _modeller.mod_restraints_unpick_all(rsr)
        span = self.__update_schedule()
        AutoCommands.pick_restraints(self, residue_span_range=span)
    def align(self):
        self.__default_file('rr_file', '.A', ('', '.mat'))
        AutoCommands.align(self)
    def sequence_search(self):
        self.__default_file('rr_file', '.mat', ('', '.mat'))
        AutoCommands.sequence_search(self)
    def malign(self):
        self.__default_file('rr_file', '.A', ('', '.mat'))
        AutoCommands.malign(self)
    def principal_components(self):
        self.__default_file('matrix_file', '.G')
        self.__default_file('file', '.dat')
        AutoCommands.principal_components(self)
    def align2d(self):
        self.__default_file('rr_file', '.A', ('', '.mat'))
        AutoCommands.align2d(self)
    def write_pdb_xref(self):
        self.__default_file('file', '.X')
        (file, model_segment) = self.__seg_from_aln('model_segment')
        AutoCommands.write_pdb_xref(self, model_segment=model_segment,
                                    file=self.__add_outdir('file'))
    def build_profile(self):
        self.__default_file('rr_file', '.mat', ('', '.mat'))
        AutoCommands.build_profile(self)
    def read_sequence_db(self):
        chains = self.vars['chains_list']
        if chains.upper() != 'ALL':
            self.__default_file('chains_list', 'H', ('', '.list'))
        self.__default_file('seq_database_file', 'H', ('', '.seq'))
        AutoCommands.read_sequence_db(self)
    def write_sequence_db(self):
        self.__default_file('chains_list', 'H')
        self.__default_file('seq_database_file', 'H')
        AutoCommands.write_sequence_db(self)
    def seqfilter(self):
        self.__default_file('rr_file', '.mat', ('', '.mat'))
        AutoCommands.seqfilter(self)
    def open(self):
        self.__default_file('objects_file', 'TOP')
        AutoCommands.open(self)
    def inquire(self):
        self.__default_file('file', '.B')
        AutoCommands.inquire(self)
    def delete_alignment(self):
        self.aln[0] = None
    def read_alignment(self):
        if not self.vars['add_sequence']:
            _modeller.mod_alignment_clear(self.get_aln())
        self.__default_file('file', '.A', ('', '.ali'))
        AutoCommands.read_alignment(self, allow_alternates=False)
    def read_alignment2(self):
        if not self.vars['add_sequence']:
            _modeller.mod_alignment_clear(self.get_aln(2))
        self.__default_file('file', '.A', ('', '.ali'))
        AutoCommands.read_alignment2(self)
    def prof_to_aln(self):
        if not self.vars['append_aln']:
            _modeller.mod_alignment_clear(self.get_aln())
        AutoCommands.prof_to_aln(self)
    def read_parameters(self):
        if not self.vars['add_parameters']:
            # Clear both CHARMM parameters and group restraint parameters
            _modeller.mod_parameters_clear(self.get_prm())
            _modeller.mod_group_restraints_clear(self.get_gprsr())
        self.__default_file('file', '.toplib', ('', '.lib'))
        AutoCommands.read_parameters(self)
    def delete_restraint(self):
        atom_names = self.vars['atom_ids']
        mdl = self.get_mdl(1)
        cd = _modeller.mod_model_cd_get(mdl)
        seq = _modeller.mod_model_seq_get(mdl)
        atom_ids = [_modeller.mod_coordinates_find_atom(cd, seq, a) \
                    for a in atom_names]
        if 0 in atom_ids:
            log.warning("delete_restraint",
                        "One or more atoms absent from MODEL:  " \
                        + " ".join(atom_names))
        else:
            AutoCommands.unpick_restraints(self, atom_ids=atom_ids)
    def patch(self):
        mdl = self.get_mdl()
        resids = self.__parse_residue_ids(mdl, self.vars['residue_ids'])
        AutoCommands.patch(self, residue_ids=resids)
    def superpose(self):
        AutoCommands.superpose(self, rms_cutoff=self.vars['rms_cutoffs'][0])
    def salign(self):
        AutoCommands.salign(self, rms_cutoff=self.vars['rms_cutoffs'][0])
    def compare(self):
        AutoCommands.compare(self, varatom=self.vars['distance_atoms'][0])
    def make_chains(self):
        AutoCommands.make_chains(self, chop_nonstd_termini= \
                                       self.vars['chop_nonstd_terminii'])
    def condense_restraints(self):
        rsr = self.get_rsr()
        _modeller.mod_restraints_unpick_redundant(rsr)
        _modeller.mod_restraints_remove_unpicked(rsr)
    def check_alignment(self):
        aln = self.get_aln()
        io = self.get_io().modpt
        libs = self.get_libs()
        _modeller.mod_alignment_check_structures(aln, io, libs, 6.0)
        _modeller.mod_alignment_check_seqstruc(aln, io, libs, 8.0)
    def time_mark(self):
        _modeller.mod_time_mark()
    def system(self):
        _modeller.mod_system(self.vars['command'], '')
    def rotate_model(self):
        mdl = self.get_mdl(1)
        inds = _modeller.mod_selection_all(mdl)
        translation = self.vars['translation']
        _modeller.mod_selection_translate(mdl, inds, translation)
        matrix = self.vars['rotation_matrix']
        _modeller.mod_selection_transform(mdl, inds, matrix)
        axis = self.vars['rotation_axis']
        angle = self.vars['rotation_angle']
        _modeller.mod_selection_rotate_axis(mdl, inds, axis, angle)
