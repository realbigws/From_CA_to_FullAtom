from modeller.energy_data import energy_data
from modeller.io_data import io_data
import _modeller

class AutoCommands:
    """Auto-generated methods for many TOP commands"""
    aln = [None, None]
    mdl = [None, None]
    sdb = [None]
    prf = [None]
    den = [None]
    libs = None
    edat = None
    io = None
    gprsr = None

    def get_topvars(self, vars, **keys):
        args = []
        if not isinstance(vars, (list, tuple)):
            vars = (vars,)
        for var in vars:
            if keys.has_key(var):
                args.append(keys[var])
            else:
                args.append(self.vars[var])
        return args

    def set_topvars(self, vars, keys):
        for i in range(len(keys)):
            self.vars[keys[i]] = vars[i]

    def get_edat(self):
        if self.edat == None:
            self.edat = energy_data()
        return self.edat

    def get_io(self):
        if self.io == None:
            self.io = io_data()
        return self.io

    def get_gprsr(self):
        if self.gprsr == None:
            self.gprsr = _modeller.mod_group_restraints_new(None)
        return self.gprsr

    def get_aln(self, num=1):
        aln = self.aln[num - 1]
        if aln == None:
            aln = _modeller.mod_alignment_new(None)
            self.aln[num - 1] = aln
        return aln

    def get_mdl(self, num=1):
        mdl = self.mdl[num - 1]
        if mdl == None:
            mdl = _modeller.mod_model_new(None)
            _modeller.mod_model_gprsr_set(mdl, self.get_libs(),
                                          self.get_gprsr())
            self.mdl[num - 1] = mdl
        return mdl

    def get_rsr(self, num=1):
        return _modeller.mod_model_rsr_get(self.get_mdl(num))

    def get_sdb(self, num=1):
        sdb = self.sdb[num - 1]
        if sdb == None:
            sdb = _modeller.mod_sequence_db_new(None)
            self.sdb[num - 1] = sdb
        return sdb

    def get_prf(self, num=1):
        prf = self.prf[num - 1]
        if prf == None:
            prf = _modeller.mod_profile_new(None)
            self.prf[num - 1] = prf
        return prf

    def get_libs(self):
        if self.libs == None:
            self.libs = _modeller.mod_libraries_new(None)
        return self.libs

    def get_tpl(self):
        return _modeller.mod_libraries_tpl_get(self.get_libs())

    def get_prm(self):
        return _modeller.mod_libraries_prm_get(self.get_libs())

    def get_den(self, num=1):
        den = self.den[num - 1]
        if den == None:
            den = _modeller.mod_density_new(None)
            self.den[num - 1] = den
        return den

    def read_alignment2(self, **keys):
        aln2 = self.get_aln(2)
        io = self.get_io().modpt
        libs = self.get_libs()
        args = self.get_topvars(('align_codes2', 'atom_files2', 'file', 'remove_gaps', 'alignment_format'), **keys)
        outs = _modeller.mod_top_read_alignment2(aln2, io, libs, *args)
        self.set_topvars((outs,), ('end_of_file',))

    def read_model2(self, **keys):
        mdl2 = self.get_mdl(2)
        io = self.get_io().modpt
        libs = self.get_libs()
        args = self.get_topvars(('file', 'model_format', 'model2_segment'), **keys)
        _modeller.mod_top_read_model2(mdl2, io, libs, *args)

    def write_model2(self, **keys):
        mdl2 = self.get_mdl(2)
        libs = self.get_libs()
        args = self.get_topvars(('file', 'model_format', 'no_ter'), **keys)
        _modeller.mod_top_write_model2(mdl2, libs, *args)

    def read_alignment(self, **keys):
        aln = self.get_aln()
        io = self.get_io().modpt
        libs = self.get_libs()
        args = self.get_topvars(('align_codes', 'atom_files', 'file', 'remove_gaps', 'alignment_format', 'allow_alternates'), **keys)
        outs = _modeller.mod_alignment_read2(aln, io, libs, *args)
        self.set_topvars((outs,), ('end_of_file',))

    def compare_alignments(self, **keys):
        aln = self.get_aln()
        aln2 = self.get_aln(2)
        outs = _modeller.mod_alignment_compare_with(aln, aln2)
        self.set_topvars((outs,), ('pct_equiv',))

    def sequence_to_ali(self, **keys):
        aln = self.get_aln()
        mdl = self.get_mdl()
        libs = self.get_libs()
        args = self.get_topvars(('atom_files', 'align_codes'), **keys)
        _modeller.mod_alignment_append_model(aln, mdl, libs, *args)

    def write_alignment(self, **keys):
        aln = self.get_aln()
        libs = self.get_libs()
        args = self.get_topvars(('fh', 'alignment_format', 'alignment_features', 'align_block', 'align_alignment'), **keys)
        _modeller.mod_alignment_write(aln, libs, *args)

    def edit_alignment(self, **keys):
        aln = self.get_aln()
        io = self.get_io().modpt
        libs = self.get_libs()
        args = self.get_topvars(('overhang', 'edit_align_codes', 'base_align_codes', 'min_base_entries'), **keys)
        _modeller.mod_alignment_edit(aln, io, libs, *args)

    def describe(self, **keys):
        aln = self.get_aln()
        io = self.get_io().modpt
        libs = self.get_libs()
        _modeller.mod_alignment_describe(aln, io, libs)

    def id_table(self, **keys):
        aln = self.get_aln()
        args = self.get_topvars(('matrix_file'), **keys)
        _modeller.mod_id_table(aln, *args)

    def sequence_comparison(self, **keys):
        aln = self.get_aln()
        mdl = self.get_mdl()
        libs = self.get_libs()
        args = self.get_topvars(('rr_file', 'matrix_file', 'variability_file', 'max_gaps_match'), **keys)
        _modeller.mod_compare_sequences(aln, mdl, libs, *args)

    def align(self, **keys):
        aln = self.get_aln()
        libs = self.get_libs()
        args = self.get_topvars(('off_diagonal', 'local_alignment', 'matrix_offset', 'gap_penalties_1d', 'read_weights', 'write_weights', 'n_subopt', 'subopt_offset', 'weigh_sequences', 'smooth_prof_weight', 'align_what', 'weights_type', 'input_weights_file', 'output_weights_file', 'rr_file', 'overhang', 'align_block'), **keys)
        _modeller.mod_align(aln, libs, *args)

    def align2d(self, **keys):
        aln = self.get_aln()
        io = self.get_io().modpt
        libs = self.get_libs()
        args = self.get_topvars(('overhang', 'align_block', 'rr_file', 'align_what', 'off_diagonal', 'max_gap_length', 'local_alignment', 'matrix_offset', 'gap_penalties_1d', 'gap_penalties_2d', 'surftyp', 'fit', 'fix_offsets', 'read_weights', 'write_weights', 'input_weights_file', 'output_weights_file', 'n_subopt', 'subopt_offset', 'read_profile', 'input_profile_file', 'write_profile', 'output_profile_file', 'weigh_sequences', 'smooth_prof_weight', 'weights_type'), **keys)
        _modeller.mod_align2d(aln, io, libs, *args)

    def malign(self, **keys):
        aln = self.get_aln()
        libs = self.get_libs()
        args = self.get_topvars(('rr_file', 'off_diagonal', 'local_alignment', 'matrix_offset', 'overhang', 'align_block', 'gap_penalties_1d'), **keys)
        _modeller.mod_malign(aln, libs, *args)

    def align_consensus(self, **keys):
        aln = self.get_aln()
        libs = self.get_libs()
        args = self.get_topvars(('align_block', 'gap_penalties_1d', 'read_weights', 'write_weights', 'weigh_sequences', 'input_weights_file', 'output_weights_file', 'weights_type', 'smooth_prof_weight'), **keys)
        _modeller.mod_alignment_consensus(aln, libs, *args)

    def compare(self, **keys):
        aln = self.get_aln()
        edat = self.get_edat().modpt
        io = self.get_io().modpt
        libs = self.get_libs()
        args = self.get_topvars(('compare_mode', 'fit', 'fit_atoms', 'matrix_file', 'output', 'asgl_output', 'refine_local', 'rms_cutoffs', 'varatom'), **keys)
        _modeller.mod_compare_structures(aln, edat, io, libs, *args)

    def align3d(self, **keys):
        aln = self.get_aln()
        io = self.get_io().modpt
        libs = self.get_libs()
        args = self.get_topvars(('off_diagonal', 'overhang', 'local_alignment', 'matrix_offset', 'gap_penalties_3d', 'fit', 'fit_atoms', 'align3d_trf', 'output', 'align3d_repeat'), **keys)
        _modeller.mod_align3d(aln, io, libs, *args)

    def malign3d(self, **keys):
        aln = self.get_aln()
        io = self.get_io().modpt
        libs = self.get_libs()
        args = self.get_topvars(('off_diagonal', 'overhang', 'local_alignment', 'matrix_offset', 'gap_penalties_3d', 'fit', 'fit_atoms', 'output', 'write_whole_pdb', 'current_directory', 'write_fit', 'edit_file_ext'), **keys)
        _modeller.mod_malign3d(aln, io, libs, *args)

    def salign(self, **keys):
        aln = self.get_aln()
        io = self.get_io().modpt
        libs = self.get_libs()
        args = self.get_topvars(('residue_type2', 'no_ter', 'overhang', 'off_diagonal', 'matrix_offset', 'gap_penalties_1d', 'gap_penalties_2d', 'gap_penalties_3d', 'feature_weights', 'rms_cutoff', 'fit', 'surftyp', 'fit_on_first', 'gap_function', 'align_block', 'max_gap_length', 'align_what', 'read_weights', 'write_weights', 'input_weights_file', 'output_weights_file', 'weigh_sequences', 'smooth_prof_weight', 'fix_offsets', 'substitution', 'comparison_type', 'matrix_comparison', 'alignment_type', 'edit_file_ext', 'weights_type', 'similarity_flag', 'bkgrnd_prblty_file', 'ext_tree_file', 'dendrogram_file', 'matrix_scaling_factor', 'auto_overhang', 'overhang_factor', 'overhang_auto_limit', 'local_alignment', 'improve_alignment', 'fit_atoms', 'output', 'write_whole_pdb', 'current_directory', 'write_fit', 'fit_pdbnam', 'rr_file', 'n_subopt', 'subopt_offset', 'align3d_trf', 'normalize_pp_scores', 'gap_gap_score', 'gap_residue_score', 'nsegm', 'matrix_offset_3d'), **keys)
        outs = _modeller.mod_salign(aln, io, libs, *args)
        self.set_topvars(outs, ('aln_score', 'qscorepct'))

    def aln_to_prof(self, **keys):
        aln = self.get_aln()
        prf = self.get_prf()
        libs = self.get_libs()
        _modeller.mod_profile_from_aln(aln, prf, libs)

    def segment_matching(self, **keys):
        aln = self.get_aln()
        libs = self.get_libs()
        args = self.get_topvars(('rr_file', 'file', 'align_block', 'min_loop_length', 'segment_shifts', 'segment_growth_n', 'segment_growth_c', 'segment_cutoff', 'segment_report', 'root_name', 'file_id', 'file_ext'), **keys)
        _modeller.mod_segment_matching(aln, libs, *args)

    def read_density(self, **keys):
        den = self.get_den()
        args = self.get_topvars(('file', 'em_density_format', 'em_map_size', 'filter_type', 'voxel_size', 'resolution', 'filter_values', 'density_type', 'px', 'py', 'pz', 'cc_func_type'), **keys)
        _modeller.mod_density_read(den, *args)

    def em_grid_search(self, **keys):
        den = self.get_den()
        args = self.get_topvars(('em_density_format', 'num_structures', 'dock_order', 'start_type', 'translate_type', 'number_of_steps', 'angular_step_size', 'temperature', 'best_docked_models', 'em_fit_output_file', 'em_pdb_name', 'chains_num'), **keys)
        _modeller.mod_density_grid_search(den, *args)

    def switch_trace(self, **keys):
        mdl = self.get_mdl()
        args = self.get_topvars(('file', 'output_directory'), **keys)
        _modeller.mod_top_switch_trace(mdl, *args)

    def assess_model(self, **keys):
        mdl = self.get_mdl()
        libs = self.get_libs()
        args = self.get_topvars(('assess_method'), **keys)
        outs = _modeller.mod_assess_model(mdl, libs, *args)
        self.set_topvars((outs,), ('molpdf',))

    def dendrogram(self, **keys):
        args = self.get_topvars(('matrix_file', 'cluster_cut'), **keys)
        _modeller.mod_dendrogram(*args)

    def principal_components(self, **keys):
        args = self.get_topvars(('matrix_file', 'file'), **keys)
        _modeller.mod_principal_components(*args)

    def delete_file(self, **keys):
        args = self.get_topvars(('file'), **keys)
        _modeller.mod_file_delete(*args)

    def inquire(self, **keys):
        args = self.get_topvars(('file'), **keys)
        outs = _modeller.mod_inquire(*args)
        self.set_topvars((outs,), ('file_exists',))

    def read_atom_classes(self, **keys):
        gprsr = self.get_gprsr()
        args = self.get_topvars(('atom_classes_file'), **keys)
        _modeller.mod_atom_classes_read(gprsr, *args)

    def time_mark(self, **keys):
        _modeller.mod_time_mark()

    def read_model(self, **keys):
        mdl = self.get_mdl()
        io = self.get_io().modpt
        libs = self.get_libs()
        args = self.get_topvars(('file', 'model_format', 'model_segment'), **keys)
        _modeller.mod_model_read(mdl, io, libs, *args)

    def write_model(self, **keys):
        mdl = self.get_mdl()
        libs = self.get_libs()
        args = self.get_topvars(('sel1', 'fh', 'model_format', 'no_ter', 'write_all_atoms'), **keys)
        _modeller.mod_model_write(mdl, libs, *args)

    def generate_topology(self, **keys):
        mdl = self.get_mdl()
        aln = self.get_aln()
        io = self.get_io().modpt
        libs = self.get_libs()
        args = self.get_topvars(('iseq', 'patch_default'), **keys)
        _modeller.mod_model_topology_generate(mdl, aln, io, libs, *args)

    def patch(self, **keys):
        mdl = self.get_mdl()
        libs = self.get_libs()
        args = self.get_topvars(('residue_type', 'residue_ids'), **keys)
        _modeller.mod_model_patch(mdl, libs, *args)

    def patch_ss_templates(self, **keys):
        mdl = self.get_mdl()
        aln = self.get_aln()
        io = self.get_io().modpt
        libs = self.get_libs()
        _modeller.mod_model_patch_ss_templates(mdl, aln, io, libs)

    def patch_ss_model(self, **keys):
        mdl = self.get_mdl()
        libs = self.get_libs()
        _modeller.mod_model_patch_ss(mdl, libs)

    def build_model(self, **keys):
        mdl = self.get_mdl()
        libs = self.get_libs()
        args = self.get_topvars(('build_method', 'initialize_xyz'), **keys)
        _modeller.mod_model_build(mdl, libs, *args)

    def transfer_xyz(self, **keys):
        mdl = self.get_mdl()
        aln = self.get_aln()
        io = self.get_io().modpt
        libs = self.get_libs()
        args = self.get_topvars(('cluster_cut', 'cluster_method'), **keys)
        _modeller.mod_transfer_xyz(mdl, aln, io, libs, *args)

    def transfer_res_numb(self, **keys):
        mdl = self.get_mdl()
        mdl2 = self.get_mdl(2)
        aln = self.get_aln()
        libs = self.get_libs()
        _modeller.mod_model_res_num_from(mdl, mdl2, aln, libs)

    def rename_segments(self, **keys):
        mdl = self.get_mdl()
        args = self.get_topvars(('segment_ids', 'renumber_residues'), **keys)
        _modeller.mod_model_rename_segments(mdl, *args)

    def pick_atoms(self, **keys):
        mdl = self.get_mdl()
        aln = self.get_aln()
        libs = self.get_libs()
        args = self.get_topvars(('sel', 'pick_atoms_set', 'sphere_radius', 'selection_slab', 'gap_extension', 'minmax_loop_length', 'res_types', 'atom_types', 'selection_mode', 'selection_search', 'selection_status', 'selection_from', 'sphere_center', 'selection_segment'), **keys)
        outs = _modeller.mod_top_pick_atoms(mdl, aln, libs, *args)
        self.set_topvars((outs,), ('selout',))

    def iupac_model(self, **keys):
        mdl = self.get_mdl()
        libs = self.get_libs()
        _modeller.mod_model_to_iupac(mdl, libs)

    def reorder_atoms(self, **keys):
        mdl = self.get_mdl()
        libs = self.get_libs()
        _modeller.mod_model_reorder_atoms(mdl, libs)

    def orient_model(self, **keys):
        mdl = self.get_mdl()
        outs = _modeller.mod_model_orient(mdl)
        self.set_topvars(outs, ('rotation', 'translation'))

    def write_data(self, **keys):
        mdl = self.get_mdl()
        edat = self.get_edat().modpt
        libs = self.get_libs()
        args = self.get_topvars(('surftyp', 'neighbor_cutoff', 'accessibility_type', 'output', 'file', 'probe_radius', 'psa_integration_step', 'dnr_accpt_lib'), **keys)
        _modeller.mod_model_write_data(mdl, edat, libs, *args)

    def write_pdb_xref(self, **keys):
        mdl = self.get_mdl()
        libs = self.get_libs()
        args = self.get_topvars(('file', 'model_segment'), **keys)
        _modeller.mod_top_write_pdb_xref(mdl, libs, *args)

    def make_region(self, **keys):
        mdl = self.get_mdl()
        libs = self.get_libs()
        args = self.get_topvars(('atom_accessibility', 'region_size'), **keys)
        _modeller.mod_model_make_region(mdl, libs, *args)

    def color_aln_model(self, **keys):
        mdl = self.get_mdl()
        aln = self.get_aln()
        _modeller.mod_model_color(mdl, aln)

    def make_chains(self, **keys):
        mdl = self.get_mdl()
        libs = self.get_libs()
        args = self.get_topvars(('file', 'structure_types', 'minimal_resolution', 'minimal_chain_length', 'max_nonstdres', 'chop_nonstd_termini', 'minimal_stdres', 'alignment_format'), **keys)
        _modeller.mod_make_chains(mdl, libs, *args)

    def optimize(self, **keys):
        mdl = self.get_mdl()
        edat = self.get_edat().modpt
        libs = self.get_libs()
        args = self.get_topvars(('sel1', 'optimization_method', 'init_velocities', 'temperature', 'md_time_step', 'min_atom_shift', 'cap_atom_shift', 'trace_output', 'equilibrate', 'max_iterations', 'residue_span_range', 'md_return', 'output'), **keys)
        outs = _modeller.mod_top_optimize(mdl, edat, libs, *args)
        self.set_topvars((outs,), ('molpdf',))

    def read_profile(self, **keys):
        prf = self.get_prf()
        libs = self.get_libs()
        args = self.get_topvars(('file', 'profile_format'), **keys)
        _modeller.mod_profile_read(prf, libs, *args)

    def write_profile(self, **keys):
        prf = self.get_prf()
        libs = self.get_libs()
        args = self.get_topvars(('file', 'profile_format'), **keys)
        _modeller.mod_profile_write(prf, libs, *args)

    def prof_to_aln(self, **keys):
        prf = self.get_prf()
        aln = self.get_aln()
        libs = self.get_libs()
        _modeller.mod_profile_to_aln(prf, aln, libs)

    def profile_profile_scan(self, **keys):
        prf = self.get_prf()
        libs = self.get_libs()
        args = self.get_topvars(('profile_list_file', 'matrix_offset', 'profile_format', 'rr_file', 'gap_penalties_1d', 'matrix_scaling_factor', 'max_aln_evalue', 'aln_base_filename', 'score_statistics', 'output_alignments', 'output_scores', 'output_score_file', 'pssm_weights_type', 'write_summary', 'summary_file', 'ccmatrix_offset', 'score_type'), **keys)
        _modeller.mod_profile_scan(prf, pssmdb, libs, *args)

    def build_profile(self, **keys):
        prf = self.get_prf()
        sdb = self.get_sdb()
        libs = self.get_libs()
        args = self.get_topvars(('gap_penalties_1d', 'matrix_offset', 'rr_file', 'n_prof_iterations', 'max_aln_evalue', 'matrix_scaling_factor', 'check_profile', 'output_scores', 'output_score_file', 'gaps_in_target', 'score_statistics', 'pssm_weights_type', 'write_pssm', 'pssm_file'), **keys)
        _modeller.mod_profile_build(prf, sdb, libs, *args)

    def make_restraints(self, **keys):
        mdl = self.get_mdl()
        edat = self.get_edat().modpt
        aln = self.get_aln()
        io = self.get_io().modpt
        libs = self.get_libs()
        args = self.get_topvars(('sel1', 'sel2', 'sel3', 'residue_span_range', 'restraint_type', 'restraint_sel_atoms', 'restraint_group', 'basis_pdf_weight', 'distance_rsr_model', 'maximal_distance', 'basis_relative_weight', 'spline_on_site', 'residue_span_sign', 'intersegment', 'dih_lib_only', 'spline_dx', 'spline_min_points', 'spline_range', 'mnch_lib', 'accessibility_type', 'restraint_stdev', 'restraint_stdev2', 'surftyp', 'distngh', 'exclude_distance'), **keys)
        _modeller.mod_restraints_make(mdl, edat, aln, io, libs, *args)

    def pick_restraints(self, **keys):
        mdl = self.get_mdl()
        args = self.get_topvars(('sel1', 'residue_span_range', 'restraint_sel_atoms', 'restraints_filter'), **keys)
        _modeller.mod_restraints_pick(mdl, *args)

    def add_restraint(self, **keys):
        mdl = self.get_mdl()
        args = self.get_topvars(('restraint_parameters', 'atom_ids'), **keys)
        _modeller.mod_restraints_top_add(mdl, *args)

    def unpick_restraints(self, **keys):
        mdl = self.get_mdl()
        args = self.get_topvars(('atom_ids'), **keys)
        _modeller.mod_restraints_unpick(mdl, *args)

    def reindex_restraints(self, **keys):
        mdl = self.get_mdl()
        mdl2 = self.get_mdl(2)
        _modeller.mod_restraints_reindex(mdl, mdl2)

    def spline_restraints(self, **keys):
        mdl = self.get_mdl()
        edat = self.get_edat().modpt
        libs = self.get_libs()
        args = self.get_topvars(('spline_dx', 'spline_range', 'spline_select', 'spline_min_points', 'output'), **keys)
        _modeller.mod_restraints_spline(mdl, edat, libs, *args)

    def read_restraints(self, **keys):
        mdl = self.get_mdl()
        args = self.get_topvars(('file'), **keys)
        _modeller.mod_restraints_read(mdl, *args)

    def write_restraints(self, **keys):
        mdl = self.get_mdl()
        args = self.get_topvars(('file'), **keys)
        _modeller.mod_restraints_write(mdl, *args)

    def define_symmetry(self, **keys):
        mdl = self.get_mdl()
        args = self.get_topvars(('sel2', 'sel3', 'add_symmetry', 'symmetry_weight'), **keys)
        _modeller.mod_symmetry_define(mdl, *args)

    def make_schedule(self, **keys):
        mdl = self.get_mdl()
        args = self.get_topvars(('library_schedule', 'residue_span_range', 'schedule_scale'), **keys)
        _modeller.mod_top_make_schedule(mdl, *args)

    def read_schedule(self, **keys):
        mdl = self.get_mdl()
        args = self.get_topvars(('file', 'schedule_scale'), **keys)
        _modeller.mod_top_read_schedule(mdl, *args)

    def write_schedule(self, **keys):
        mdl = self.get_mdl()
        args = self.get_topvars(('file', 'output_directory'), **keys)
        _modeller.mod_top_write_schedule(mdl, *args)

    def mutate_model(self, **keys):
        mdl = self.get_mdl()
        libs = self.get_libs()
        args = self.get_topvars(('sel1', 'residue_type'), **keys)
        _modeller.mod_selection_mutate(mdl, libs, *args)

    def randomize_xyz(self, **keys):
        mdl = self.get_mdl()
        libs = self.get_libs()
        args = self.get_topvars(('sel1', 'deviation'), **keys)
        _modeller.mod_randomize_xyz(mdl, libs, *args)

    def superpose(self, **keys):
        mdl = self.get_mdl()
        mdl2 = self.get_mdl(2)
        aln = self.get_aln()
        libs = self.get_libs()
        args = self.get_topvars(('sel1', 'swap_atoms_in_res', 'reference_atom', 'reference_distance', 'superpose_refine', 'fit', 'refine_local', 'rms_cutoff'), **keys)
        outs = _modeller.mod_superpose(mdl, mdl2, aln, libs, *args)
        self.set_topvars(outs, ('rms_initial', 'rms_final', 'drms_final', 'rotation', 'translation', 'num_equiv_pos', 'num_equiv_dist', 'num_equiv_cutoff_pos', 'num_equiv_cutoff_dist', 'rms_cutoff_final', 'drms_cutoff_final'))

    def rotate_dihedrals(self, **keys):
        mdl = self.get_mdl()
        libs = self.get_libs()
        args = self.get_topvars(('sel1', 'deviation', 'change', 'dihedrals'), **keys)
        _modeller.mod_rotate_dihedrals(mdl, libs, *args)

    def unbuild_model(self, **keys):
        mdl = self.get_mdl()
        args = self.get_topvars(('sel1'), **keys)
        _modeller.mod_selection_unbuild(mdl, *args)

    def pick_hot_atoms(self, **keys):
        mdl = self.get_mdl()
        edat = self.get_edat().modpt
        libs = self.get_libs()
        args = self.get_topvars(('sel1', 'residue_span_range', 'pick_hot_cutoff', 'viol_report_cut', 'schedule_scale'), **keys)
        outs = _modeller.mod_selection_hot_atoms(mdl, edat, libs, *args)
        self.set_topvars((outs,), ('sel1out',))

    def energy(self, **keys):
        mdl = self.get_mdl()
        edat = self.get_edat().modpt
        libs = self.get_libs()
        args = self.get_topvars(('sel1', 'asgl_output', 'normalize_profile', 'residue_span_range', 'output', 'file', 'smoothing_window', 'viol_report_cut', 'viol_report_cut2', 'schedule_scale'), **keys)
        outs = _modeller.mod_energy(mdl, edat, libs, *args)
        self.set_topvars(outs, ('molpdf', 'violtyp'))

    def debug_function(self, **keys):
        mdl = self.get_mdl()
        edat = self.get_edat().modpt
        libs = self.get_libs()
        args = self.get_topvars(('sel1', 'residue_span_range', 'debug_function_cutoff', 'detailed_debugging', 'schedule_scale'), **keys)
        _modeller.mod_debug_function(mdl, edat, libs, *args)

    def read_sequence_db(self, **keys):
        sdb = self.get_sdb()
        libs = self.get_libs()
        args = self.get_topvars(('chains_list', 'seq_database_file', 'seq_database_format', 'clean_sequences', 'minmax_db_seq_len'), **keys)
        _modeller.mod_sequence_db_read(sdb, libs, *args)

    def write_sequence_db(self, **keys):
        sdb = self.get_sdb()
        libs = self.get_libs()
        args = self.get_topvars(('chains_list', 'seq_database_file', 'seq_database_format'), **keys)
        _modeller.mod_sequence_db_write(sdb, libs, *args)

    def sequence_search(self, **keys):
        sdb = self.get_sdb()
        aln = self.get_aln()
        io = self.get_io().modpt
        libs = self.get_libs()
        args = self.get_topvars(('search_randomizations', 'search_top_list', 'off_diagonal', 'overhang', 'gap_penalties_1d', 'signif_cutoff', 'rr_file', 'matrix_offset', 'fast_search_cutoff', 'data_file', 'search_group_list', 'search_sort', 'output', 'alignment_features', 'seq_database_file', 'local_alignment', 'fast_search'), **keys)
        _modeller.mod_sequence_db_search(sdb, aln, io, libs, *args)

    def seqfilter(self, **keys):
        sdb = self.get_sdb()
        libs = self.get_libs()
        args = self.get_topvars(('gap_penalties_1d', 'matrix_offset', 'rr_file', 'seqid_cut', 'max_diff_res', 'output_grp_file', 'output_cod_file'), **keys)
        _modeller.mod_sequence_db_filter(sdb, libs, *args)

    def write(self, **keys):
        args = self.get_topvars(('io_unit', 'objects'), **keys)
        _modeller.mod_top_write(*args)

    def read(self, **keys):
        args = self.get_topvars(('io_unit'), **keys)
        outs = _modeller.mod_top_read(*args)
        self.set_topvars((outs,), ('record',))

    def open(self, **keys):
        args = self.get_topvars(('io_unit', 'objects_file', 'output_directory', 'file_status', 'file_access'), **keys)
        outs = _modeller.mod_top_open(*args)
        self.set_topvars((outs,), ('number_lines',))

    def close(self, **keys):
        args = self.get_topvars(('io_unit'), **keys)
        _modeller.mod_top_close(*args)

    def stop(self, **keys):
        _modeller.mod_top_stop()

    def expand_alignment(self, **keys):
        aln = self.get_aln()
        args = self.get_topvars(('expand_control', 'root_name', 'file_id', 'file_ext'), **keys)
        _modeller.mod_top_expand_alignment(aln, *args)

    def read_parameters(self, **keys):
        prm = self.get_prm()
        gprsr = self.get_gprsr()
        libs = self.get_libs()
        args = self.get_topvars(('file'), **keys)
        _modeller.mod_top_read_parameters(prm, gprsr, libs, *args)

    def read_restyp_lib(self, **keys):
        libs = self.get_libs()
        args = self.get_topvars(('restyp_lib_file'), **keys)
        _modeller.mod_read_restyp_lib(libs, *args)

    def read_topology(self, **keys):
        tpl = self.get_tpl()
        libs = self.get_libs()
        args = self.get_topvars(('file'), **keys)
        _modeller.mod_topology_read(tpl, libs, *args)

    def make_topology_model(self, **keys):
        tpl = self.get_tpl()
        _modeller.mod_topology_model_make(tpl)

    def write_topology_model(self, **keys):
        tpl = self.get_tpl()
        libs = self.get_libs()
        args = self.get_topvars(('file'), **keys)
        _modeller.mod_topology_model_write(tpl, libs, *args)
