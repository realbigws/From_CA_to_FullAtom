/** \file mod_core.h       Core Modeller functions.
 *
 *             Part of MODELLER, Copyright(c) 1989-2008 Andrej Sali
 */

/* This file is auto-generated: any changes will be lost */

#ifndef MOD_CORE_H
#define MOD_CORE_H

#include "mod_types.h"

#ifdef __cplusplus
extern "C" {
#endif

/**  Initialize the Modeller system. */
void mod_start(int *ierr);

/**  Initialize Modeller, check the number of command line arguments
 *   (print a usage message if incorrect) and return the directories
 *   containing Modeller library and dynamic library files. */
void mod_init_cmdline(int numarg, char **libdir2, char **dynlibdir2, int *ierr);

/**  Return the directory containing Modeller .h files. */
char * mod_incdir_get(void);

/**  Return the directory containing Modeller dynamic library files. */
char * mod_dynlibdir_get(void);

/**  Shut down Modeller, at the end of a job. */
void mod_end(void);

/**  Create and return a new libraries object, tied to the given
 *   scripting language object. */
struct mod_libraries * mod_libraries_new(const void *scriptobj);

/**  Release the memory used by a libraries object. */
void mod_libraries_free(struct mod_libraries *libs);

/**  Get the initial random seed. */
int mod_libraries_rand_seed_get(const struct mod_libraries *libs);

/**  Set the initial random seed. */
void mod_libraries_rand_seed_set(struct mod_libraries *libs, int rand_seed);

/**  Read in most libraries (e.g. residue types, classes) */
void mod_libraries_read_libs(struct mod_libraries *libs, const char *restyp_lib_file, int *ierr);

/**  Get the next random number between 0 and 1. */
float mod_random_number(struct mod_libraries *libs);

/**  Perturb the random number generator by some amount, n. */
void mod_random_perturb(struct mod_libraries *libs, int n);

/**  Clustering. */
void mod_principal_components(const char *matrix_file, const char *file, int *ierr);

/**  Clustering. */
void mod_dendrogram(const char *matrix_file, float cluster_cut, int *ierr);

/**  Read residue type file; deprecated, for TOP interpreter. */
void mod_read_restyp_lib(struct mod_libraries *libs, const char *restyp_lib_file, int *ierr);

/** Make a database of PSSMs */
void mod_pssmdb_make(const struct mod_libraries *libs, const char *profile_list_file, const char *profile_format, const char *rr_file, float matrix_offset, float matrix_scaling_factor, const char *pssmdb_name, const char *pssm_weights_type, int *ierr);

/**  Create and return a new group_restraints object, tied to the given
 *   scripting language object. */
struct mod_group_restraints * mod_group_restraints_new(const void *scriptobj);

/**  Release the memory used by a group_restraints object. */
void mod_group_restraints_free(struct mod_group_restraints *gprsr);

/**  Read group restraint parameters. */
void mod_group_restraints_read(struct mod_group_restraints *gprsr, const struct mod_libraries *libs, const char *file, int *ierr);

/**  Clear group restraint parameters. */
void mod_group_restraints_clear(struct mod_group_restraints *gprsr);

/**  Read atom classes for group restraints. */
void mod_atom_classes_read(struct mod_group_restraints *gprsr, const char *atom_classes_file, int *ierr);

/**  Create and return a new alignment object, tied to the given
 *   scripting language object. */
struct mod_alignment * mod_alignment_new(const void *scriptobj);

/**  Release the memory used by an alignment object. */
void mod_alignment_free(struct mod_alignment *aln);

/**  Get the alignment code of the indx'th sequence. */
char * mod_alignment_codes_get(const struct mod_alignment *aln, int indx);

/**  Set the alignment code of the indx'th sequence. */
void mod_alignment_codes_set(struct mod_alignment *aln, int indx, const char *code);

/**  Get the indx'th comment. */
char * mod_alignment_comment_get(const struct mod_alignment *aln, int indx);

/**  Set the indx'th comment. */
void mod_alignment_comment_set(struct mod_alignment *aln, int indx, const char *val);

/**  Set the number of comments. */
void mod_alignment_ncomment_set(struct mod_alignment *aln, int val);

/**  Get the index of the named sequence. */
int mod_alignment_find_code(const struct mod_alignment *aln, const char *code);

/**  Move one structure to another. */
/**  Move one sequence to another. */
/**  Move one alnsequence to another. */
/**  Move one alignment sequence to another. */
/**  Delete the sequence(s) from the alignment. */
void mod_alnsequence_del(struct mod_alignment *aln, const int indices[], int n_indices);

/**  Set the type of of the ires'th residue in the iseq'th sequence. */
void mod_alnresidue_type_set(struct mod_alignment *aln, int ires, int iseq, int irestyp, const struct mod_libraries *libs);

/**  Add gaps after the ires'th residue in the iseq'th sequence. */
void mod_alnresidue_add_gaps(struct mod_alignment *aln, int ires, int iseq, int ngap);

/**  Remove gaps after the ires'th residue in the iseq'th sequence. */
void mod_alnresidue_remove_gaps(struct mod_alignment *aln, int ires, int iseq, int ngap);

/**  Check to make sure the seq'th sequence in the alignment matches the
 *   given model. */
void mod_alnsequence_check_model(const struct mod_alignment *aln, int seq, const struct mod_model *mdl, const struct mod_libraries *libs, int *ierr);

/**  Should structural information be read in for this sequence? */
gboolean mod_alnsequence_has_structure(struct mod_alignment *aln, int iseq);

/**  Read a single template structure from file. */
void mod_alnstructure_read(struct mod_alignment *aln, int iseq, const struct mod_io_data *io, struct mod_libraries *libs, int *ierr);

/**  Write a template structure's current coordinates to a PDB file. */
void mod_alnstructure_write(const struct mod_coordinates *cd, const struct mod_sequence *seq, struct mod_file *fh, const struct mod_libraries *libs, int *ierr);

/**  Invalidate all structural information for a template. The next time
 *   this information is requested, it will be re-read from the atom file. */
void mod_alnstructure_invalidate(struct mod_alignment *aln, int iseq);

/**  Get the residue range of the next insertion in the last sequence in
 *   the alignment. */
void mod_alignment_next_insert(const struct mod_alignment *aln, int pos, int *pos_out, int *ins_st, int *ins_end, int minlength, int maxlength, int extension);

/**  Get the residue range of the next deletion in the last sequence in
 *   the alignment. */
void mod_alignment_next_delete(const struct mod_alignment *aln, int pos, int *pos_out, int *ins_st, int *ins_end, int extension);

/**  Edit overhangs in alignment. */
void mod_alignment_edit(struct mod_alignment *aln, const struct mod_io_data *io, const struct mod_libraries *libs, int overhang, const char edit_align_codes[], int edit_align_codes_len, int n_edit_align_codes, const char base_align_codes[], int base_align_codes_len, int n_base_align_codes, int min_base_entries, int *ierr);

/**  Align two or more sequences/structures of proteins. */
void mod_salign(struct mod_alignment *aln, const struct mod_io_data *io, struct mod_libraries *libs, const char *residue_type2, gboolean no_ter, int overhang, int off_diagonal, float matrix_offset, const float gap_penalties_1d[2], const float gap_penalties_2d[9], const float gap_penalties_3d[2], const float feature_weights[6], float rms_cutoff, gboolean fit, int surftyp, gboolean fit_on_first, gboolean gap_function, int align_block, int max_gap_length, const char *align_what, gboolean read_weights, gboolean write_weights, const char *input_weights_file, const char *output_weights_file, gboolean weigh_sequences, float smooth_prof_weight, const float fix_offsets[5], gboolean substitution, const char *comparison_type, const char *matrix_comparison, const char *alignment_type, const char edit_file_ext[2], int edit_file_ext_len, const char *weights_type, float *aln_score, gboolean similarity_flag, const char *bkgrnd_prblty_file, const char *ext_tree_file, const char *dendrogram_file, float matrix_scaling_factor, gboolean auto_overhang, float overhang_factor, int overhang_auto_limit, gboolean local_alignment, gboolean improve_alignment, const char *fit_atoms, const char *output, gboolean write_whole_pdb, gboolean current_directory, gboolean write_fit, gboolean fit_pdbnam, const char *rr_file, int n_subopt, float subopt_offset, gboolean align3d_trf, gboolean normalize_pp_scores, float gap_gap_score, float gap_residue_score, float *qscorepct, int nsegm, float matrix_offset_3d, int *ierr);

/**  Align two or more sequences. */
void mod_malign(struct mod_alignment *aln, const struct mod_libraries *libs, const char *rr_file, int off_diagonal, gboolean local_alignment, float matrix_offset, int overhang, int align_block, const float gap_penalties_1d[2], int *ierr);

/**  Compare structures. */
void mod_compare_structures(struct mod_alignment *aln, const struct mod_energy_data *edat, const struct mod_io_data *io, struct mod_libraries *libs, int compare_mode, gboolean fit, const char *fit_atoms, const char *matrix_file, const char *output, gboolean asgl_output, gboolean refine_local, const float rms_cutoffs[11], const char *varatom, int *ierr);

/**  Read alignment. Compatibility function to avoid breaking the API. */
void mod_alignment_read(struct mod_alignment *aln, const struct mod_io_data *io, const struct mod_libraries *libs, const char align_codes[], int align_codes_len, int n_align_codes, const char atom_files[], int atom_files_len, int n_atom_files, const char *file, gboolean remove_gaps, const char *alignment_format, int *end_of_file, int *ierr);

/**  Read alignment. To be renamed mod_alignment_read at next API version
 *   bump. */
void mod_alignment_read2(struct mod_alignment *aln, const struct mod_io_data *io, const struct mod_libraries *libs, const char align_codes[], int align_codes_len, int n_align_codes, const char atom_files[], int atom_files_len, int n_atom_files, const char *file, gboolean remove_gaps, const char *alignment_format, int *end_of_file, gboolean allow_alternates, int *ierr);

/**  Read a single sequence from an alignment file.
 *   Compatibility function to avoid breaking the API. */
void mod_alignment_read_one(struct mod_alignment *aln, struct mod_file *fh, const struct mod_io_data *io, const struct mod_libraries *libs, gboolean remove_gaps, const char *alignment_format, int *end_of_file, int *ierr);

/**  Read a single sequence from an alignment file (PIR or FASTA).
 *   To be renamed mod_alignment_read_one at next API version bump. */
void mod_alignment_read_one2(struct mod_alignment *aln, struct mod_file *fh, const struct mod_io_data *io, const struct mod_libraries *libs, gboolean remove_gaps, const char *alignment_format, gboolean allow_alternates, int *end_of_file, int *ierr);

/**  Clear alignment. */
void mod_alignment_clear(struct mod_alignment *aln);

/**  Write alignment to file. */
void mod_alignment_write(struct mod_alignment *aln, const struct mod_libraries *libs, struct mod_file *fh, const char *alignment_format, const char *alignment_features, int align_block, gboolean align_alignment, int *ierr);

/**  Align two (blocks of) sequences. */
void mod_align(struct mod_alignment *aln, const struct mod_libraries *libs, int off_diagonal, gboolean local_alignment, float matrix_offset, const float gap_penalties_1d[2], gboolean read_weights, gboolean write_weights, int n_subopt, float subopt_offset, gboolean weigh_sequences, float smooth_prof_weight, const char *align_what, const char *weights_type, const char *input_weights_file, const char *output_weights_file, const char *rr_file, int overhang, int align_block, int *ierr);

/**  Align two structures. */
void mod_align3d(struct mod_alignment *aln, const struct mod_io_data *io, struct mod_libraries *libs, int off_diagonal, int overhang, gboolean local_alignment, float matrix_offset, const float gap_penalties_3d[2], gboolean fit, const char *fit_atoms, gboolean align3d_trf, const char *output, gboolean align3d_repeat, int *ierr);

/**  Describe proteins. */
void mod_alignment_describe(struct mod_alignment *aln, const struct mod_io_data *io, struct mod_libraries *libs, int *ierr);

/**  Compare sequences in alignment. */
void mod_compare_sequences(const struct mod_alignment *aln, struct mod_model *mdl, const struct mod_libraries *libs, const char *rr_file, const char *matrix_file, const char *variability_file, int max_gaps_match, int *ierr);

/**  Align two or more structures. */
void mod_malign3d(struct mod_alignment *aln, const struct mod_io_data *io, struct mod_libraries *libs, int off_diagonal, int overhang, gboolean local_alignment, float matrix_offset, const float gap_penalties_3d[2], gboolean fit, const char *fit_atoms, const char *output, gboolean write_whole_pdb, gboolean current_directory, gboolean write_fit, const char edit_file_ext[2], int edit_file_ext_len, int *ierr);

/**  Copy model sequence and coordinates to alignment. */
void mod_alignment_append_model(struct mod_alignment *aln, struct mod_model *mdl, const struct mod_libraries *libs, const char *atom_files, const char *align_codes);

/**  Compare two alignments */
void mod_alignment_compare_with(const struct mod_alignment *aln, const struct mod_alignment *aln2, float *pct_equiv, int *ierr);

/**  Consensus sequence alignment. */
void mod_alignment_consensus(struct mod_alignment *aln, const struct mod_libraries *libs, int align_block, const float gap_penalties_1d[2], gboolean read_weights, gboolean write_weights, gboolean weigh_sequences, const char *input_weights_file, const char *output_weights_file, const char *weights_type, float smooth_prof_weight, int *ierr);

/**  Check pairwise structural superpositions of all template structures. */
void mod_alignment_check_structures(struct mod_alignment *aln, const struct mod_io_data *io, struct mod_libraries *libs, float eqvdst, int *n_exceed, int *ierr);

/**  Check the current sequence/structure alignment for sanity. */
void mod_alignment_check_seqstruc(struct mod_alignment *aln, const struct mod_io_data *io, struct mod_libraries *libs, float gapdist, int *n_exceed, int *ierr);

/**  Align sequences with structures. */
void mod_align2d(struct mod_alignment *aln, const struct mod_io_data *io, struct mod_libraries *libs, int overhang, int align_block, const char *rr_file, const char *align_what, int off_diagonal, int max_gap_length, gboolean local_alignment, float matrix_offset, const float gap_penalties_1d[2], const float gap_penalties_2d[9], int surftyp, gboolean fit, const float fix_offsets[5], gboolean read_weights, gboolean write_weights, const char *input_weights_file, const char *output_weights_file, int n_subopt, float subopt_offset, gboolean read_profile, const char *input_profile_file, gboolean write_profile, const char *output_profile_file, gboolean weigh_sequences, float smooth_prof_weight, const char *weights_type, int *ierr);

/**  Align segments. */
void mod_segment_matching(struct mod_alignment *aln, const struct mod_libraries *libs, const char *rr_file, const char *file, int align_block, const int min_loop_length[], int n_min_loop_length, const int segment_shifts[], int n_segment_shifts, const int segment_growth_n[], int n_segment_growth_n, const int segment_growth_c[], int n_segment_growth_c, float segment_cutoff, int segment_report, const char *root_name, const char *file_id, const char *file_ext, int *ierr);

/**  Calculate percentage sequence identities. */
void mod_id_table(const struct mod_alignment *aln, const char *matrix_file, int *ierr);

/**  Transfer residue properties of predicted secondary structure
 *   from one sequence to all other sequences. */
void mod_transfer_res_prop(struct mod_alignment *aln, int iseq);

/**  Add a sequence (given by one-letter codes) to an alignment. */
void mod_alignment_append_sequence(struct mod_alignment *aln, const char *newseq, const struct mod_libraries *libs, int *ierr);

/**  Create and return a new model object, tied to the given
 *   scripting language object. */
struct mod_model * mod_model_new(const void *scriptobj);

/**  Release the memory used by a model object. */
void mod_model_free(struct mod_model *mdl);

/**  Set the group_restraints object used by the model. */
void mod_model_gprsr_set(struct mod_model *mdl, const struct mod_libraries *libs, const struct mod_group_restraints *val);

/**  Remove any group_restraints object from the model. */
void mod_model_gprsr_unset(struct mod_model *mdl, const struct mod_libraries *libs);

/**  Return the next set of atom indices which matches the given residue
 *   type, atom names, and residue offsets. */
void mod_model_find_atoms(const struct mod_model *mdl, const char *resnam, const char atms[], int atms_len, int n_atms, const int resoffs[], int n_resoffs, int start, int *ires, int **indatm, int *n_indatm, const struct mod_libraries *libs);

/**  Return the next set of atom indices which defines the given dihedral
 *   type (phi, psi, omega, chi12345). */
void mod_model_find_dihedrals(const struct mod_model *mdl, const char *resnam, int idiht, int start, int *ires, int **indatm, int *n_indatm, const struct mod_libraries *libs);

/** Calculate the mainchain RMSDs */
void mod_model_fast_rmsd(const struct mod_model *mdl, const struct mod_model *mdl2, float *rmsd, int *ierr);

/**  Define a random surface patch of atoms. */
void mod_model_make_region(struct mod_model *mdl, struct mod_libraries *libs, float atom_accessibility, int region_size, int *ierr);

/**  Guess model disulfides from model structure. */
void mod_model_patch_ss(struct mod_model *mdl, const struct mod_libraries *libs, int *ierr);

/**  Write derivative model data. */
void mod_model_write_data(struct mod_model *mdl, const struct mod_energy_data *edat, const struct mod_libraries *libs, int surftyp, float neighbor_cutoff, int accessibility_type, const char *output, const char *file, float probe_radius, float psa_integration_step, const char *dnr_accpt_lib, int *ierr);

/**  Guess model disulfides from templates. */
void mod_model_patch_ss_templates(struct mod_model *mdl, struct mod_alignment *aln, const struct mod_io_data *io, struct mod_libraries *libs, int *ierr);

/**  Read the atom model file; READ_MODEL */
void mod_model_read(struct mod_model *mdl, const struct mod_io_data *io, struct mod_libraries *libs, const char *file, const char *model_format, const char model_segment[2], int model_segment_len, int *ierr);

/**  Write the model */
void mod_model_write(const struct mod_model *mdl, const struct mod_libraries *libs, const int sel1[], int n_sel1, struct mod_file *fh, const char *model_format, gboolean no_ter, gboolean write_all_atoms, int *ierr);

/**  Patch the MODEL */
void mod_model_patch(struct mod_model *mdl, const struct mod_libraries *libs, const char *residue_type, const int residue_ids[], int n_residue_ids, int *ierr);

/**  Transfer residue numbers */
void mod_model_res_num_from(struct mod_model *mdl, struct mod_model *mdl2, const struct mod_alignment *aln, const struct mod_libraries *libs, int *ierr);

/**  Build model coordinates from topology. */
void mod_model_build(struct mod_model *mdl, struct mod_libraries *libs, const char *build_method, gboolean initialize_xyz, int *ierr);

/**  Build model coordinates from internal coordinates. */
void mod_model_build_ic(struct mod_model *mdl, struct mod_libraries *libs, gboolean initialize_xyz, int *ierr);

/**  Clear model topology. */
void mod_model_topology_clear(struct mod_model *mdl);

/**  Generate model topology. */
void mod_model_topology_generate(struct mod_model *mdl, struct mod_alignment *aln, const struct mod_io_data *io, struct mod_libraries *libs, int iseq, gboolean patch_default, int *ierr);

/**  Standardize order of model atoms. */
void mod_model_reorder_atoms(struct mod_model *mdl, struct mod_libraries *libs, int *ierr);

/**  Copy templates' coordinates to the model. */
void mod_transfer_xyz(struct mod_model *mdl, struct mod_alignment *aln, const struct mod_io_data *io, struct mod_libraries *libs, float cluster_cut, const char *cluster_method, int *ierr);

/**  Center and orient the model. */
void mod_model_orient(struct mod_model *mdl, float rotation[9], float translation[3], int *ierr);

/**  Rename model segments and residues. */
void mod_model_rename_segments(struct mod_model *mdl, const char segment_ids[], int segment_ids_len, int n_segment_ids, const int renumber_residues[], int n_renumber_residues, int *ierr);

/**  Color model according to alignment. */
void mod_model_color(struct mod_model *mdl, const struct mod_alignment *aln, int *ierr);

/**  Standardize certain dihedral angles. */
void mod_model_to_iupac(struct mod_model *mdl, const struct mod_libraries *libs, int *ierr);

/**  Define similar segments. */
void mod_symmetry_define(struct mod_model *mdl, const int sel2[], int n_sel2, const int sel3[], int n_sel3, const gboolean add_symmetry[2], float symmetry_weight, int *ierr);

/**  Calculate protein volume. */
void mod_volume(struct mod_model *mdl, const struct mod_energy_data *edat, const struct mod_libraries *libs, const char *file, const char *output, float grid_unit, int number_of_steps, gboolean orient, int *ierr);

/**  Calculate protein cavity volume. */
void mod_volume_cavity(struct mod_model *mdl, const struct mod_energy_data *edat, const struct mod_libraries *libs, const char *file, const char *output, float grid_unit, int number_of_steps, float probe_radius, float rcutp, float rcutl, float vmin, float rlink, gboolean orient, gboolean accuracy_border, int *ierr);

/**  Create and return a new sequence_db object, tied to the given
 *   scripting language object. */
struct mod_sequence_db * mod_sequence_db_new(const void *scriptobj);

/**  Release the memory used by a sequence_db object. */
void mod_sequence_db_free(struct mod_sequence_db *sdb);

/**  Get the alignment code of the indx'th sequence. */
char * mod_sequence_db_code_get(const struct mod_sequence_db *sdb, int indx);

/**  Set the alignment code of the indx'th sequence. */
void mod_sequence_db_code_set(struct mod_sequence_db *sdb, int indx, const char *code);

/**  Get the protein type of the indx'th sequence. */
char * mod_sequence_db_prottyp_get(const struct mod_sequence_db *sdb, int indx);

/**  Set the protein type of the indx'th sequence. */
void mod_sequence_db_prottyp_set(struct mod_sequence_db *sdb, int indx, const char *prottyp);

/**  Get the residue integer type. */
int mod_sequence_db_restype_get(const struct mod_sequence_db *sdb, int seq, int indx);

/**  Set the residue integer type. */
void mod_sequence_db_restype_set(struct mod_sequence_db *sdb, int seq, int indx, int val);

/**  Cluster sequences by sequence-identity. */
void mod_sequence_db_filter(const struct mod_sequence_db *sdb, const struct mod_libraries *libs, const float gap_penalties_1d[2], float matrix_offset, const char *rr_file, int seqid_cut, int max_diff_res, const char *output_grp_file, const char *output_cod_file, int *ierr);

/**  Read a database of sequences. */
void mod_sequence_db_read(struct mod_sequence_db *sdb, const struct mod_libraries *libs, const char *chains_list, const char *seq_database_file, const char *seq_database_format, gboolean clean_sequences, const int minmax_db_seq_len[2], int *ierr);

/**  Write a database of sequences. */
void mod_sequence_db_write(struct mod_sequence_db *sdb, const struct mod_libraries *libs, const char *chains_list, const char *seq_database_file, const char *seq_database_format, int *ierr);

/**  Convert a database of sequences to binary format. */
void mod_sequence_db_convert(struct mod_sequence_db *sdb, const struct mod_libraries *libs, const char *chains_list, const char *seq_database_file, const char *seq_database_format, const char *outfile, gboolean clean_sequences, const int minmax_db_seq_len[2], int *ierr);

/**  Search for similar sequences. */
void mod_sequence_db_search(struct mod_sequence_db *sdb, struct mod_alignment *aln, const struct mod_io_data *io, struct mod_libraries *libs, int search_randomizations, int search_top_list, int off_diagonal, int overhang, const float gap_penalties_1d[2], const float signif_cutoff[2], const char *rr_file, float matrix_offset, float fast_search_cutoff, gboolean data_file, const char *search_group_list, const char *search_sort, const char *output, const char *alignment_features, const char *seq_database_file, gboolean local_alignment, gboolean fast_search, int *ierr);

/**  Create and return a new profile object, tied to the given
 *   scripting language object. */
struct mod_profile * mod_profile_new(const void *scriptobj);

/**  Release the memory used by a profile object. */
void mod_profile_free(struct mod_profile *prf);

/**  Get the profile sequence alignment code. */
char * mod_profile_code_get(const struct mod_profile *prf, int indx);

/**  Set the profile sequence alignment code. */
void mod_profile_code_set(struct mod_profile *prf, int indx, const char *code);

/**  Get the profile sequence protein type. */
char * mod_profile_prottyp_get(const struct mod_profile *prf, int indx);

/**  Set the profile sequence protein type. */
void mod_profile_prottyp_set(struct mod_profile *prf, int indx, const char *prottyp);

/**  Build a profile for a given sequence or alignment. */
void mod_profile_build(struct mod_profile *prf, const struct mod_sequence_db *sdb, const struct mod_libraries *libs, const float gap_penalties_1d[2], float matrix_offset, const char *rr_file, int n_prof_iterations, float max_aln_evalue, float matrix_scaling_factor, gboolean check_profile, gboolean output_scores, const char *output_score_file, gboolean gaps_in_target, gboolean score_statistics, const char *pssm_weights_type, gboolean write_pssm, const char *pssm_file, int *ierr);

/**  Compare a target profile against a database of profiles. */
void mod_profile_scan(struct mod_profile *prf, const struct mod_pssmdbobj *pssmdb, const struct mod_libraries *libs, const char *profile_list_file, float matrix_offset, const char *profile_format, const char *rr_file, const float gap_penalties_1d[2], float matrix_scaling_factor, float max_aln_evalue, const char *aln_base_filename, gboolean score_statistics, gboolean output_alignments, gboolean output_scores, const char *output_score_file, const char *pssm_weights_type, gboolean write_summary, const char *summary_file, float ccmatrix_offset, const char *score_type, int *ierr);

/**  Read a profile of a sequence. */
void mod_profile_read(struct mod_profile *prf, const struct mod_libraries *libs, const char *file, const char *profile_format, int *ierr);

/**  Write a profile. */
void mod_profile_write(const struct mod_profile *prf, const struct mod_libraries *libs, const char *file, const char *profile_format, int *ierr);

/**  Convert alignment to profile format. */
void mod_profile_from_aln(const struct mod_alignment *aln, struct mod_profile *prf, const struct mod_libraries *libs, int *ierr);

/**  Convert profile to alignment format. */
void mod_profile_to_aln(const struct mod_profile *prf, struct mod_alignment *aln, const struct mod_libraries *libs);

/**  Delete all restraints. */
void mod_restraints_clear(struct mod_model *mdl);

/**  Unselect all restraints. */
void mod_restraints_unpick_all(struct mod_restraints *rsr);

/**  Add a single restraint to the list. */
void mod_restraints_add(struct mod_restraints *rsr, int iform, int group, const int iftyp[], int n_iftyp, const int modal[], int n_modal, const int natm[], int n_natm, const int indatm[], int n_indatm, const float pcsr[], int n_pcsr, gboolean use_array, int *array_index, int *ierr);

/**  Get one restraint. */
void mod_restraints_get(const struct mod_restraints *rsr, int indx, int *group, int **modal, int *n_modal, int **indatm, int *n_indatm, int **feat, int *n_feat, float **param, int *n_param, int *ierr);

/**  Read spatial restraints. */
void mod_restraints_read(struct mod_model *mdl, const char *file, int *ierr);

/**  Write current restraints */
void mod_restraints_write(struct mod_model *mdl, const char *file, int *ierr);

/**  Make restraints. */
void mod_restraints_make(struct mod_model *mdl, struct mod_energy_data *edat, struct mod_alignment *aln, const struct mod_io_data *io, struct mod_libraries *libs, const int sel1[], int n_sel1, const int sel2[], int n_sel2, const int sel3[], int n_sel3, const int residue_span_range[2], const char *restraint_type, int restraint_sel_atoms, int restraint_group, const char *basis_pdf_weight, int distance_rsr_model, float maximal_distance, float basis_relative_weight, gboolean spline_on_site, gboolean residue_span_sign, gboolean intersegment, gboolean dih_lib_only, float spline_dx, int spline_min_points, float spline_range, int mnch_lib, int accessibility_type, const float restraint_stdev[2], const float restraint_stdev2[3], int surftyp, float distngh, float exclude_distance, int *ierr);

/**  Make alpha (helix) secondary structure restraints. */
void mod_restraints_make_alpha(struct mod_model *mdl, const int residue_ids[2], const struct mod_libraries *libs, int *ierr);

/**  Make strand secondary structure restraints. */
void mod_restraints_make_strand(struct mod_model *mdl, const int residue_ids[2], const struct mod_libraries *libs, int *ierr);

/**  Make sheet secondary structure restraints. */
void mod_restraints_make_sheet(struct mod_model *mdl, const int atom_ids[2], int sheet_h_bonds, int *ierr);

/**  Pick restraints for selected atoms. */
void mod_restraints_pick(struct mod_model *mdl, const int sel1[], int n_sel1, const int residue_span_range[2], int restraint_sel_atoms, const float restraints_filter[], int n_restraints_filter, int *ierr);

/**  Remove unselected restraints. */
void mod_restraints_remove_unpicked(struct mod_restraints *rsr);

/**  Unpick redundant cosine dihedral restraints. */
void mod_restraints_unpick_redundant(struct mod_restraints *rsr, int *nunpick);

/**  Unselect restraints. */
void mod_restraints_unpick(struct mod_model *mdl, const int atom_ids[], int n_atom_ids, int *ierr);

/**  Add a restraint; deprecated interface for TOP. */
void mod_restraints_top_add(struct mod_model *mdl, const float restraint_parameters[], int n_restraint_parameters, const char atom_ids[], int atom_ids_len, int n_atom_ids, int *ierr);

/**  Renumber mdl2 restraints for mdl. */
void mod_restraints_reindex(struct mod_model *mdl, const struct mod_model *mdl2, int *ierr);

/**  Approximate restraints by splines. */
void mod_restraints_spline(struct mod_model *mdl, struct mod_energy_data *edat, const struct mod_libraries *libs, float spline_dx, float spline_range, const int spline_select[3], int spline_min_points, const char *output, int *ierr);

/**  Create and return a new density object, tied to the given
 *   scripting language object. */
struct mod_density * mod_density_new(const void *scriptobj);

/**  Release the memory used by a density object. */
void mod_density_free(struct mod_density *den);

/**  Calculate the forces on the selected atoms in the given model,
 *   resulting purely from the fit of the model to the density. */
void mod_density_forces(struct mod_density *den, struct mod_model *mdl, float *ccf);

/**  Read an EM density file */
void mod_density_read(struct mod_density *den, const char *file, const char *em_density_format, int em_map_size, const char *filter_type, float voxel_size, float resolution, const float filter_values[2], const char *density_type, float px, float py, float pz, const char *cc_func_type, int *ierr);

/**  EM fitting by grid search */
void mod_density_grid_search(const struct mod_density *den, const char *em_density_format, int num_structures, const char *dock_order, const char *start_type, const char *translate_type, int number_of_steps, float angular_step_size, float temperature, int best_docked_models, const char *em_fit_output_file, const char em_pdb_name[], int em_pdb_name_len, int n_em_pdb_name, const int chains_num[], int n_chains_num, int *ierr);

/**  Set residue span range from user-specified value or schedule definition. */
void mod_schedule_update(struct mod_schedule *sched, const int residue_span_range[2], int *nrang1, int *nrang2, int *ierr);

/**  Check whether there is a schedule in memory. */
void mod_schedule_check(struct mod_schedule *sched, int *ierr);

/**  Get the schedule's current optimizer. */
int mod_schedule_optimizer_get(const struct mod_schedule *sch);

/**  Create and return a new energy_data object. */
struct mod_energy_data * mod_energy_data_new(void);

/**  Release the memory used by an energy_data object. */
void mod_energy_data_free(struct mod_energy_data *ene);

/**  Set an energy_data object's associated EM density. */
void mod_energy_data_density_set(struct mod_energy_data *ene, const struct mod_density *den);

/**  Unset an energy_data object's associated EM density. */
void mod_energy_data_density_unset(struct mod_energy_data *ene);

/**  Create and return a new io_data object. */
struct mod_io_data * mod_io_data_new(void);

/**  Release the memory used by an io_data object. */
void mod_io_data_free(struct mod_io_data *io);

/**  Clear all parameter information. */
void mod_parameters_clear(struct mod_parameters *prm);

/**  Return True if parameters have been read in already. */
gboolean mod_parameters_in_memory(const struct mod_parameters *prm);

/**  Read parameters from file. */
void mod_parameters_read(struct mod_parameters *prm, const struct mod_libraries *libs, const char *file, int *ierr);

/**  Clear all topology information. */
void mod_topology_clear(struct mod_topology *tpl, const struct mod_libraries *libs);

/**  Return True if the topology library has already been read. */
gboolean mod_topology_in_memory(const struct mod_topology *tpl);

/**  Read residue topology library from file. */
void mod_topology_read(struct mod_topology *tpl, const struct mod_libraries *libs, const char *file, int *ierr);

/**  Write residue topology library to file. */
void mod_topology_model_write(const struct mod_topology *tpl, const struct mod_libraries *libs, const char *file, int *ierr);

/**  Make a subset topology library. */
void mod_topology_model_make(struct mod_topology *tpl, int *ierr);

/**  Get the element symbol. */
char * mod_topology_element_get(const struct mod_topology *tpl, int indx);

/**  Set the element symbol. */
void mod_topology_element_set(struct mod_topology *tpl, int indx, const char *val);

/**  Get the atom type. */
char * mod_coordinates_atmnam_get(const struct mod_coordinates *cd, int indx);

/**  Get the element symbol. */
char * mod_coordinates_element_get(const struct mod_coordinates *cd, int indx);

/**  Set the element symbol. */
void mod_coordinates_element_set(struct mod_coordinates *cd, int indx, const char *val);

/**  Get the PDB residue number of the given residue. */
char * mod_coordinates_resnum_get(const struct mod_coordinates *cd, int indx);

/**  Set the PDB residue number of the given residue. */
void mod_coordinates_resnum_set(struct mod_coordinates *cd, int indx, const char *val);

/**  Get the current value of a dihedral (in degrees) of the given type.
 *   \note Does not currently account for disulfide bridges. */
float mod_coordinates_dihedral_get(const struct mod_coordinates *cd, const struct mod_sequence *seq, const struct mod_libraries *libs, int idihtyp, int indx, int *ierr);

/**  Given a dihedral value in degrees, return the dihedral class. */
int mod_sequence_dihclass_get(const struct mod_sequence *seq, const struct mod_libraries *libs, int idihtyp, float dih, int indx);

/**  Return True iff the given residue has the given dihedral type. */
gboolean mod_coordinates_has_dihedral(const struct mod_coordinates *cd, const struct mod_sequence *seq, const struct mod_libraries *libs, int idihtyp, int indx, int *ierr);

/**  Get the atom indices of the given dihedral type. */
void mod_coordinates_dihatoms_get(const struct mod_coordinates *cd, const struct mod_sequence *seq, const struct mod_libraries *libs, int idihtyp, int indx, int **indatm, int *n_indatm, int *ierr);

/**  Get the integer type of a given atom. */
int mod_coordinates_atom_type_get(const struct mod_coordinates *cd, const struct mod_sequence *seq, int indx, const struct mod_libraries *libs, int *ierr);

/**  Get the index of a named atom. */
int mod_coordinates_find_atom(const struct mod_coordinates *cd, const struct mod_sequence *seq, const char *atm, int *ierr);

/**  Get the index of a named residue. */
int mod_coordinates_find_residue(const struct mod_coordinates *cd, const struct mod_sequence *seq, const char *res);

/**  Find the given atom name in the specified residue, honoring +/-
 *   specifiers if given. */
int mod_residue_find_atom(const struct mod_coordinates *cd, const struct mod_sequence *seq, int ir, const char *atmnam);

/**  Calculate the average sidechain Biso for the given residue.
 *   The number of sidechain atoms is also returned. */
void mod_residue_sidechain_biso(const struct mod_coordinates *cd, const struct mod_sequence *seq, int ir, float *biso, int *n_sdch_atoms);

/**  Create and return a new saxsdata object. */
struct mod_saxsdata * mod_saxsdata_new(void);

/**  Release the memory used by a saxsdata object. */
void mod_saxsdata_free(struct mod_saxsdata *saxsd);

/**  Calculate SAXS intensity from model.
 *   last change 10/39/06 FF - now works again on model */
void mod_saxs_intens(const struct mod_model *mdl, struct mod_saxsdata *saxsd, const char *filename, gboolean fitflag, int *ierr);

/**  Calculate SAXS score from model. */
void mod_saxs_chifun(struct mod_energy_data *enedata, const struct mod_model *mdl, gboolean transfer_is, int *ierr);

/**   Calculate and P(r) */
void mod_saxs_pr(const struct mod_model *mdl, struct mod_saxsdata *saxsd, const char *filename, int *ierr);

/**  Initialization of SAXS data.
 *   In this routine the SAXS parameters are set, the sampling in reciprocal
 *   space is determined and the resulting scattering factors are computed.
 *   last change 04/10/06
 *   10/10/06 FF - note: sel1 and n_sel1 have to be called like that to be
 *                 passed correctly ...
 *   01/28/07 FF - set nr_exp and dr_exp in ini_saxs
 *   03/14/07 FF - included rolloff */
void mod_saxs_ini(struct mod_saxsdata *saxsd, const struct mod_model *mdl, const int sel1[], int n_sel1, double s_min, double s_max, int maxs, int nmesh, int natomtyp, const char *represtyp, const char *filename, const char *wswitch, double s_hybrid, double s_low, double s_hi, const char *spaceflag, double rho_solv, gboolean use_lookup, int nr, double dr, int nr_exp, double dr_exp, gboolean use_offset, gboolean use_rolloff, gboolean use_conv, gboolean mixflag, gboolean pr_smooth, int *ierr);

/**  Read in experimental (or simulated) SAXS data */
void mod_saxs_read(struct mod_saxsdata *saxsd, const char *filename, int *ierr);

/**  Read in experimental P(r) to saxsstructure */
void mod_saxs_pr_read(struct mod_saxsdata *saxsd, const char *filename, int *ierr);

/**  Insert a new saxsdata pointer before |indx|. */
void mod_saxsdata_pt_new(struct mod_energy_data *edat, int indx, const struct mod_saxsdata *saxsd);

/**  Remove the saxsdata pointer at position |indx| in the list. */
void mod_saxsdata_pt_del(struct mod_energy_data *edat, int indx);

/**  Create and return a new PSSM DB object, tied to the given scripting
 *   language object. */
struct mod_pssmdbobj * mod_pssmdbobj_new(const void *scriptobj);

/**  Release the memory used by a PSSM DB object. */
void mod_pssmdbobj_free(struct mod_pssmdbobj *psm);

/** Read in a database of PSSMs */
void mod_pssmdbobj_read(struct mod_pssmdbobj *psm, const char *pssmdb_name, const char *pssmdb_format, int *ierr);

/**  Return True iff the given feature type is an angle. */
gboolean mod_feature_isangle(int ifeat, int *ierr);

/**  Evaluate the given feature type on the given atoms. */
float mod_feature_eval(const struct mod_model *mdl, const int indatm[], int n_indatm, int iftyp, int *ierr);

/**  Evaluate the first derivatives of the given feature. */
void mod_feature_deriv(const struct mod_model *mdl, const int indatm[], int n_indatm, int iftyp, float feat, float dervx[], float dervy[], float dervz[], int *ierr);

/**  Calculate the distance between the current feature value and the
 *   mean (dealing correctly with periodic features). */
float mod_feature_delta(float feat, float rmean, int iftyp, int *ierr);

/**  Return the center of mass (unweighted) of the given selection. */
void mod_selection_com_get(const struct mod_model *mdl, const int indatm[], int n_indatm, float *cx, float *cy, float *cz);

/**  Translate all atoms in the given selection by vector. */
void mod_selection_translate(struct mod_model *mdl, const int indatm[], int n_indatm, const float vector[3]);

/**  Rotate the given selection about an axis through the origin,
 *   by the given angle (in degrees). */
void mod_selection_rotate_axis(struct mod_model *mdl, const int indatm[], int n_indatm, const float axis[3], float angle);

/**  Transform the selection coordinates with the given matrix. */
void mod_selection_transform(struct mod_model *mdl, const int indatm[], int n_indatm, const float matrix[9]);

/**  Return a new selection containing only atoms of the given
 *   space-separated types. */
void mod_selection_atom_types(const struct mod_model *mdl, const int indatm[], int n_indatm, const char *atom_types, int **newindatm, int *n_newindatm);

/**  Return a new selection containing only atoms in residues of the
 *   given space-separated types. */
void mod_selection_residue_types(const struct mod_model *mdl, const int indatm[], int n_indatm, const char *residue_types, int **newindatm, int *n_newindatm, const struct mod_libraries *libs);

/**  Return a new selection containing only atoms in HETATM residues. */
void mod_selection_het_residues(const struct mod_model *mdl, const int indatm[], int n_indatm, int **newindatm, int *n_newindatm, const struct mod_libraries *libs);

/**  Return a new selection containing only atoms in water residues. */
void mod_selection_water_residues(const struct mod_model *mdl, const int indatm[], int n_indatm, int **newindatm, int *n_newindatm, const struct mod_libraries *libs);

/**  Return a new selection containing only atoms in residue types with
 *   no defined topology. */
void mod_selection_no_topology(const struct mod_model *mdl, const int indatm[], int n_indatm, int **newindatm, int *n_newindatm, const struct mod_libraries *libs);

/**  Return a new selection containing only atoms in standard residue
 *   types. */
void mod_selection_std_residues(const struct mod_model *mdl, const int indatm[], int n_indatm, int **newindatm, int *n_newindatm, const struct mod_libraries *libs);

/**  Return a new selection containing only atoms with defined coordinates. */
void mod_selection_defined(const struct mod_model *mdl, const int indatm[], int n_indatm, int **newindatm, int *n_newindatm);

/**  Return a new selection containing all atoms within the given radius
 *   of the given point. */
void mod_selection_sphere(const struct mod_model *mdl, float x, float y, float z, float radius, int **newindatm, int *n_newindatm);

/**  Return a selection containing all atoms from the given model. */
void mod_selection_all(const struct mod_model *mdl, int **newindatm, int *n_newindatm);

/**  Superpose MODEL2 on top of MODEL */
void mod_superpose(struct mod_model *mdl, struct mod_model *mdl2, const struct mod_alignment *aln, const struct mod_libraries *libs, const int sel1[], int n_sel1, const char *swap_atoms_in_res, const char *reference_atom, float reference_distance, gboolean superpose_refine, gboolean fit, gboolean refine_local, float rms_cutoff, float *rms_initial, float *rms_final, float *drms_final, float rotation[9], float translation[3], int *num_equiv_pos, int *num_equiv_dist, int *num_equiv_cutoff_pos, int *num_equiv_cutoff_dist, float *rms_cutoff_final, float *drms_cutoff_final, int *ierr);

/**  Report energy */
void mod_energy(struct mod_model *mdl, struct mod_energy_data *edat, const struct mod_libraries *libs, const int sel1[], int n_sel1, gboolean asgl_output, gboolean normalize_profile, const int residue_span_range[2], const char *output, const char *file, int smoothing_window, float *molpdf, const float viol_report_cut[], int n_viol_report_cut, const float viol_report_cut2[], int n_viol_report_cut2, const float schedule_scale[], int n_schedule_scale, float **violtyp, int *n_violtyp, int *ierr);

/**  Get just the energy (no derivatives or printouts). */
void mod_selection_objfunc(struct mod_model *mdl, struct mod_energy_data *edat, const struct mod_libraries *libs, const int sel1[], int n_sel1, const int residue_span_range[2], float *molpdf, const float schedule_scale[], int n_schedule_scale, int *ierr);

/**  Calculate energy RMS and profile for restraints of the given physical
 *   type. */
void mod_rms_profile(struct mod_model *mdl, struct mod_energy_data *edat, const struct mod_libraries *libs, const int sel1[], int n_sel1, const int residue_span_range[2], gboolean energy_profile, gboolean viol_profile, int iphytyp, const float schedule_scale[], int n_schedule_scale, float **profres, int *n_profres, int **nprofres, int *n_nprofres, float *min_rms, float *heavy_rms, int *ierr);

/**  Change model structure by rotation around dihedral angles */
void mod_rotate_dihedrals(struct mod_model *mdl, struct mod_libraries *libs, const int sel1[], int n_sel1, float deviation, const char *change, const char dihedrals[], int dihedrals_len, int n_dihedrals, int *ierr);

/**  Randomize selected coordinates. */
void mod_randomize_xyz(struct mod_model *mdl, struct mod_libraries *libs, const int sel1[], int n_sel1, float deviation, int *ierr);

/**  Test code self-consistency. The number of derivatives which exceed
 *   the cutoffs (n_exceed) is returned. */
void mod_debug_function(struct mod_model *mdl, struct mod_energy_data *edat, const struct mod_libraries *libs, const int sel1[], int n_sel1, const int residue_span_range[2], const float debug_function_cutoff[3], gboolean detailed_debugging, const float schedule_scale[], int n_schedule_scale, int *n_exceed, int *ierr);

/**  Select atoms violating restraints. */
void mod_selection_hot_atoms(struct mod_model *mdl, struct mod_energy_data *edat, const struct mod_libraries *libs, const int sel1[], int n_sel1, int **sel1out, int *n_sel1out, const int residue_span_range[2], float pick_hot_cutoff, const float viol_report_cut[], int n_viol_report_cut, const float schedule_scale[], int n_schedule_scale, int *ierr);

/**  Reset model information. */
void mod_selection_unbuild(struct mod_model *mdl, const int sel1[], int n_sel1);

/**  Create and return a new md_optimizer object, tied to the given
 *   scripting language object. */
struct mod_md_optimizer * mod_md_optimizer_new(const void *scriptobj);

/**  Release the memory used by an md_optimizer object. */
void mod_md_optimizer_free(struct mod_md_optimizer *opt);

/**  One step of MD optimization */
void mod_md_optimize(struct mod_md_optimizer *mdopt, struct mod_model *mdl, struct mod_energy_data *edat, struct mod_libraries *libs, const int sel1[], int n_sel1, gboolean init_velocities, float temperature, float md_time_step, float cap_atom_shift, int equilibrate, int max_iterations, const int residue_span_range[2], const char *md_return, float guide_factor, float guide_time, float friction, const char *output, float *molpdf, int *ierr);

/**  Get the guide force on the indx'th atom of the model. */
void mod_md_optimizer_guide_f_get(const struct mod_md_optimizer *mdopt, int indx, float *x, float *y, float *z);

/**  Create and return a new cg_optimizer object, tied to the given
 *   scripting language object. */
struct mod_cg_optimizer * mod_cg_optimizer_new(const void *scriptobj);

/**  Release the memory used by a cg_optimizer object. */
void mod_cg_optimizer_free(struct mod_cg_optimizer *opt);

/**  One step of CG optimization */
void mod_cg_optimize(struct mod_cg_optimizer *cgopt, struct mod_model *mdl, struct mod_energy_data *edat, struct mod_libraries *libs, const int sel1[], int n_sel1, float min_atom_shift, int max_iterations, const int residue_span_range[2], const char *output, float *molpdf, int *ierr);

/**  Create and return a new state_optimizer object, tied to the given
 *   scripting language object. */
struct mod_state_optimizer * mod_state_optimizer_new(const void *scriptobj);

/**  Release the memory used by a state_optimizer object. */
void mod_state_optimizer_free(struct mod_state_optimizer *opt);

/**  Prepare for an optimization using this state_optimizer object. */
void mod_state_optimizer_start(struct mod_state_optimizer *opt, struct mod_model *mdl, const struct mod_libraries *libs, struct mod_energy_data *enedata, const int inds[], int n_inds, int *ierr);

/**  Get the current state (e.g. x/y/z coordinates) of the optimizer. */
void mod_state_optimizer_state_get(const struct mod_state_optimizer *opt, float **state, int *n_state);

/**  Set the current state (e.g. x/y/z coordinates) of the optimizer. */
void mod_state_optimizer_state_set(struct mod_state_optimizer *opt, struct mod_energy_data *enedata, const float state[], int n_state, int *ierr);

/**  Calculate the energy and first derivatives of the given state. */
void mod_state_optimizer_energy(struct mod_state_optimizer *opt, struct mod_energy_data *enedata, const float state[], int n_state, float *ene, float **dstate, int *n_dstate, const struct mod_libraries *libs, int *ierr);

/**  Calculate just the energy (no first derivatives) of the given state. */
float mod_state_optimizer_energy_only(struct mod_state_optimizer *opt, struct mod_energy_data *enedata, const float state[], int n_state, const struct mod_libraries *libs, int *ierr);

/**  Increment the state_optimizer object's step counter. */
void mod_state_optimizer_next_step(struct mod_state_optimizer *opt, int *ierr);

/**  Do any cleanup at the end of optimization. */
void mod_state_optimizer_finish(struct mod_state_optimizer *opt, int *ierr);

/**  Create and return a new qn_optimizer object, tied to the given
 *   scripting language object. */
struct mod_qn_optimizer * mod_qn_optimizer_new(const void *scriptobj);

/**  Release the memory used by a qn_optimizer object. */
void mod_qn_optimizer_free(struct mod_qn_optimizer *opt);

/**  One step of quasi Newton optimization */
void mod_qn_optimize(struct mod_qn_optimizer *qnopt, struct mod_model *mdl, struct mod_energy_data *edat, struct mod_libraries *libs, const int sel1[], int n_sel1, float min_atom_shift, float max_atom_shift, int max_iterations, const int residue_span_range[2], const char *output, float *molpdf, int *ierr);

/**  Get the atom indices of a given rigid body. */
void mod_rigid_body_get(const struct mod_model *mdl, int indx, int **indatm, int *n_indatm, float *scale_factor);

/**  Set the atom indices of a given rigid body. */
void mod_rigid_body_set(struct mod_model *mdl, int indx, const int indatm[], int n_indatm, float scale_factor);

/**  Get the number of rigid bodies. */
int mod_model_nrigid_get(const struct mod_model *mdl);

/**  Set the number of rigid bodies. */
void mod_model_nrigid_set(struct mod_model *mdl, int num);

/**  Get an alignment sequence object. */
struct mod_sequence * mod_alignment_sequence_get(const struct mod_alignment *aln, int indx);

/**  Get the ID of the given segment. */
char * mod_sequence_segid_get(const struct mod_sequence *seq, int indx);

/**  Set the ID of the given segment. */
void mod_sequence_segid_set(struct mod_sequence *seq, int indx, const char *val);

/**  Get the chain index for the given residue index. */
int mod_sequence_chain_for_res(const struct mod_sequence *seq, int indx);

/**  Get the index of a named chain, or 0 if not found. */
int mod_sequence_find_chain(const char *segid, const struct mod_sequence *seq);

/**  Get the range of the sequence. */
char * mod_sequence_rng_get(const struct mod_sequence *seq, int indx);

/**  Set the range of the given segment. */
void mod_sequence_rng_set(struct mod_sequence *seq, int indx, const char *val);

/**  Does the given chain pass all filter criteria? */
gboolean mod_chain_filter(const struct mod_sequence *seq, int indx, const char *structure_types, float minimal_resolution, int minimal_chain_length, int max_nonstdres, gboolean chop_nonstd_termini, int minimal_stdres);

/**  Write the given chain out to an alignment file. */
void mod_chain_write(const struct mod_sequence *seq, const struct mod_coordinates *cd, int indx, struct mod_file *fh, const char *atom_file, const char *align_code, const char *comment, const char *format, gboolean chop_nonstd_termini, const struct mod_libraries *libs, int *ierr);

/**  Return suitable atom_file and align_codes for this chain, given
 *   a model filename. */
void mod_chain_atom_file_and_code(const struct mod_sequence *seq, int indx, const char *filename, char **atom_file, char **align_code);

/**  Remove all chain breaks between chains indx1 and indx2; all chains
 *   become part of the indx1 chain. */
void mod_chain_join(struct mod_sequence *seq, int indx1, int indx2);

/**  Given a residue type, get the internal (CHARMM) residue name
 *   (usually 4-letter). */
char * mod_residue_name_from_type(int restyp, const struct mod_libraries *libs);

/**  Given a residue type, get the PDB residue name (3-letter).
 *   If multiple synonyms are defined, return the first one. */
char * mod_residue_pdbnam_from_type(int restyp, const struct mod_libraries *libs);

/**  Given a residue type, get the one letter residue code. */
char * mod_residue_code_from_type(int restyp, const struct mod_libraries *libs);

/**  Given a residue type, get the residue group (as given in resgrp.lib).
 *   If the residue is in no group, return 0. */
int mod_residue_group_from_type(int restyp, int residue_grouping, const struct mod_libraries *libs);

/**  Return TRUE iff this residue is a PDB HETATM residue. */
gboolean mod_residue_is_hetatm(int restyp, const struct mod_libraries *libs);

/**  Map a residue type to a 'standard' one (20 standard amino acids,
 *   plus gap, as defined in the STD column in modlib/restyp.lib) */
int mod_residue_standard_type(int restyp, const struct mod_libraries *libs);

/**  Return a residue type from the internal CHARMM name. */
int mod_residue_type_from_name(const char *resnam, const struct mod_libraries *libs);

/**  Return a residue type from a PDB name. */
int mod_residue_type_from_pdbnam(const char *resnam, const struct mod_libraries *libs);

/**  Return a residue type from the one-letter code. */
int mod_residue_type_from_code(const char *resnam, const struct mod_libraries *libs);

/**  Write the header of a CHARMM-style trajectory file. */
void mod_traj_write_header(const struct mod_model *mdl, const char *file, int skip, const int sel[], int n_sel, int *unitnum, int *ierr);

/**  Write one set of coordinates to a CHARMM trajectory file. */
void mod_traj_write_set(const struct mod_model *mdl, const int sel[], int n_sel, int unitnum);

/**  Close a trajectory file. */
void mod_traj_close(int unitnum);

/**  Update the position of a single pseudo atom. */
void mod_pseudo_atom_update(struct mod_model *mdl, int indx, int *ierr);

/**  Get the type and real atom indices of a given pseudo atom. */
void mod_pseudo_atom_get(const struct mod_model *mdl, int indx, int *psdtype, int **indatm, int *n_indatm);

/**  Append a single pseudo atom to the end of the list. */
int mod_pseudo_atom_append(struct mod_model *mdl, int typ, const int indatm[], int n_indatm);

/**  Given an atom type, get the CHARMM name. */
char * mod_atom_name_from_type(int atmtyp, const struct mod_libraries *libs, int *ierr);

/**  Get the atom indices of a given excluded pair. */
void mod_excluded_pair_get(const struct mod_model *mdl, int indx, int **indatm, int *n_indatm);

/**  Set the atom indices of a given excluded pair. */
void mod_excluded_pair_set(struct mod_model *mdl, int indx, const int indatm[2]);

/**  Get the number of excluded pairs. */
int mod_model_nexcl_get(const struct mod_model *mdl);

/**  Set the number of excluded pairs. */
void mod_model_nexcl_set(struct mod_model *mdl, int num);

/**  Get the number of nonbonded pairs. */
int mod_model_npairs_get(const struct mod_model *mdl);

/**  Get the indx'th element of the 1D Fortran INTEGER array. */
int mod_int1_get(const struct mod_int1_array *ptarr, int indx);

/**  Get a C pointer to the 1D Fortran INTEGER array. */
int * mod_int1_pt(const struct mod_int1_array *ptarr);

/**  Get a C pointer to the 1D Fortran INTEGER(8) array. */
gint64 * mod_int641_pt(const struct mod_int641_array *ptarr);

/**  Set the indx'th element of the 1D Fortran INTEGER array. */
void mod_int1_set(struct mod_int1_array *ptarr, int indx, int val);

/**  Get a C pointer to the 1D Fortran REAL array. */
float * mod_float1_pt(const struct mod_float1_array *ptarr);

/**  Get the indx'th element of the 1D Fortran REAL array. */
float mod_float1_get(const struct mod_float1_array *ptarr, int indx);

/**  Set the indx'th element of the 1D Fortran REAL array. */
void mod_float1_set(struct mod_float1_array *ptarr, int indx, float val);

/**  Get a C pointer to the 1D Fortran REAL(8) array. */
double * mod_double1_pt(const struct mod_double1_array *ptarr);

/**  Get the indx'th element of the 1D Fortran REAL(8) array. */
double  mod_double1_get(const struct mod_double1_array *ptarr, int indx);

/**  Set the indx'th element of the 1D Fortran REAL(8) array. */
void mod_double1_set(struct mod_double1_array *ptarr, int indx, double val);

/**  Get a C pointer to the 1D Fortran INTEGER(1) array. */
unsigned char * mod_uchar1_pt(const struct mod_uchar1_array *ptarr);

/**  Get the indx1,indx2'th element of the 2D Fortran INTEGER array. */
int mod_int2_get(const struct mod_int2_array *ptarr, int indx1, int indx2);

/**  Set the indx1,indx2'th element of the 2D Fortran INTEGER array. */
void mod_int2_set(struct mod_int2_array *ptarr, int indx1, int indx2, int val);

/**  Get the indx1,indx2'th element of the 2D Fortran REAL array. */
float mod_float2_get(const struct mod_float2_array *ptarr, int indx1, int indx2);

/**  Set the indx1,indx2'th element of the 2D Fortran REAL array. */
void mod_float2_set(struct mod_float2_array *ptarr, int indx1, int indx2, float val);

/**  Get the indx1,indx2'th element of the 2D Fortran REAL(8) array. */
double  mod_double2_get(const struct mod_double2_array *ptarr, int indx1, int indx2);

/**  Set the indx1,indx2'th element of the 2D Fortran REAL(8) array. */
void mod_double2_set(struct mod_double2_array *ptarr, int indx1, int indx2, double val);

/**  Get the indx1..3'th element of the 3D Fortran INTEGER array. */
int mod_int3_get(const struct mod_int3_array *ptarr, int indx1, int indx2, int indx3);

/**  Set the indx1..3'th element of the 3D Fortran INTEGER array. */
void mod_int3_set(struct mod_int3_array *ptarr, int indx1, int indx2, int indx3, int val);

/**  Get the indx1..3'th element of the 3D Fortran REAL array. */
float mod_float3_get(const struct mod_float3_array *ptarr, int indx1, int indx2, int indx3);

/**  Set the indx1..3'th element of the 3D Fortran REAL array. */
void mod_float3_set(struct mod_float3_array *ptarr, int indx1, int indx2, int indx3, float val);

/**  Redimension the 3D Fortran REAL array, optionally retaining existing
 *   data. */
void mod_float3_dimension(struct mod_float3_array *ptarr, int dim1, int dim2, int dim3, gboolean retain);

/**  Get a C pointer to the 3D Fortran REAL array. */
float * mod_float3_pt(const struct mod_float3_array *ptarr);

/**  Get an alignment structure object. */
struct mod_structure * mod_alignment_structure_get(const struct mod_alignment *aln, int indx);

/**  Is the given atom a hydrogen? */
gboolean mod_atom_is_hydrogen(const char *atmnam);

/**  Return a list of all pairs of residues linked by disulfide bridges. */
void mod_find_ss(int **iss, int *nss, const struct mod_structure *struc, const struct mod_sequence *seq, int *ierr);

/**  Report violations of symmetry restraints. */
void mod_symmetry_report(const struct mod_model *mdl, float deviation);

/**  Get an alignment alnsequence object. */
struct mod_alnsequence * mod_alignment_alnsequence_get(const struct mod_alignment *aln, int indx);

/**  Initialize an mdt_library object. */
void mod_mdt_library_init(struct mod_mdt_library *mlib);

/**  Deallocate an mdt_library object. */
void mod_mdt_library_dealloc(struct mod_mdt_library *mlib);

/**  Initialize a mod_mdt object. */
void mod_mdt_init(struct mod_mdt *mdt);

/**  Deallocate a mod_mdt object. */
void mod_mdt_dealloc(struct mod_mdt *mdt);

/**  Setup and check a new MDT. */
void mod_mdt_setup_check(struct mod_mdt *mdt, const struct mod_mdt_library *mlib, int *ierr);

/**  Read all template data needed for this MDT. */
void mod_mdt_getdata(const struct mod_mdt *mdt, int *nseqacc, struct mod_alignment *aln, float distngh, gboolean sdchngh, int surftyp, int iacc1typ, const struct mod_io_data *io, struct mod_libraries *libs, int *ierr);

/**  Calculates the data needed to create or use the currently
 *   specified or read in MDT pdf. It makes the bin indices where possible
 *   and does some other things to speed up the tabulation or use of MDT. */
void mod_mdt_precalc(const struct mod_mdt *mdt, const struct mod_mdt_library *mlib, struct mod_alignment *aln, const struct mod_libraries *libs, int *ierr);

/**  Print out TOP errors, if any. */
void mod_top_error(int error_status, int stop_on_error);

/**  Write statement */
void mod_top_write(int io_unit, const char objects[], int objects_len, int n_objects, int *ierr);

/**  Read statement */
void mod_top_read(int io_unit, char **record, int *ierr);

/**  Open statement */
void mod_top_open(int io_unit, const char *objects_file, const char *output_directory, const char *file_status, const char *file_access, int *number_lines, int *ierr);

/**  Close a file */
void mod_top_close(int io_unit);

/**  Terminate Modeller */
void mod_top_stop(void);

/**  Read schedule file */
void mod_top_read_schedule(struct mod_model *mdl, const char *file, const float schedule_scale[], int n_schedule_scale, int *ierr);

/**  One variable target function step of optimization */
void mod_top_optimize(struct mod_model *mdl, struct mod_energy_data *edat, struct mod_libraries *libs, const int sel1[], int n_sel1, int optimization_method, gboolean init_velocities, float temperature, float md_time_step, float min_atom_shift, float cap_atom_shift, int trace_output, int equilibrate, int max_iterations, const int residue_span_range[2], const char *md_return, const char *output, float *molpdf, int *ierr);

/**  Pick atoms */
void mod_top_pick_atoms(struct mod_model *mdl, const struct mod_alignment *aln, const struct mod_libraries *libs, const int sel[], int n_sel, int **selout, int *n_selout, int pick_atoms_set, float sphere_radius, const float selection_slab[5], const int gap_extension[2], const int minmax_loop_length[2], const char *res_types, const char *atom_types, const char *selection_mode, const char *selection_search, const char *selection_status, const char *selection_from, const char sphere_center[2], int sphere_center_len, const char selection_segment[2], int selection_segment_len, int *ierr);

/**  Switch the trace file */
void mod_top_switch_trace(struct mod_model *mdl, const char *file, const char *output_directory, int *ierr);

/**  Write model 2 */
void mod_top_write_model2(struct mod_model *mdl2, const struct mod_libraries *libs, const char *file, const char *model_format, gboolean no_ter, int *ierr);

/**  Read the model 2 file */
void mod_top_read_model2(struct mod_model *mdl2, const struct mod_io_data *io, struct mod_libraries *libs, const char *file, const char *model_format, const char model2_segment[2], int model2_segment_len, int *ierr);

/**  Get the default variable target function schedule; DEFAULT_SCHEDULE */
void mod_top_make_schedule(struct mod_model *mdl, int library_schedule, const int residue_span_range[2], const float schedule_scale[], int n_schedule_scale, int *ierr);

/**  Write schedule file */
void mod_top_write_schedule(const struct mod_model *mdl, const char *file, const char *output_directory, int *ierr);

/**  Read parameters. */
void mod_top_read_parameters(struct mod_parameters *prm, struct mod_group_restraints *gprsr, const struct mod_libraries *libs, const char *file, int *ierr);

/**  Read the second alignment; no possibility to add sequences here. */
void mod_top_read_alignment2(struct mod_alignment *aln2, const struct mod_io_data *io, const struct mod_libraries *libs, const char align_codes2[], int align_codes2_len, int n_align_codes2, const char atom_files2[], int atom_files2_len, int n_atom_files2, const char *file, gboolean remove_gaps, const char *alignment_format, int *end_of_file, int *ierr);

/**  Put all models into alignment. */
void mod_top_expand_alignment(struct mod_alignment *aln, const int expand_control[5], const char *root_name, const char *file_id, const char *file_ext, int *ierr);

/**  Write residue number/index correspondence. */
void mod_top_write_pdb_xref(struct mod_model *mdl, const struct mod_libraries *libs, const char *file, const char model_segment[2], int model_segment_len, int *ierr);

/**  Open a file with Fortran routines.
 *
 *   Deprecated; use openf() instead, which uses C routines. */
/**  Fetch sequences from PDB file. */
void mod_make_chains(struct mod_model *mdl, const struct mod_libraries *libs, const char *file, const char *structure_types, float minimal_resolution, int minimal_chain_length, int max_nonstdres, gboolean chop_nonstd_termini, int minimal_stdres, const char *alignment_format, int *ierr);


#ifdef __cplusplus
}
#endif

#endif  /* MOD_CORE_H */
