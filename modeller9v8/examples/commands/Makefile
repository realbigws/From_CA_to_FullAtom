INSTALL=/usr/bin/install -c
CLEANFILES=*.log *.rsr segmatch.* *.ini* *.lib CHAINS*.bin *.mat *.sch *.ssm* \
           *.var toxin-* *.xref* *.clr *.fit* *.trc *.crd *_fit.pdb *.reg \
           *.iup *.fasta *.princ *.B *.ca *.noSS toxin.pap align2d.ali *.cif \
           align2d.pap 1fas-1.atm  color_aln_model.pap define_symmetry-* \
           2ctx.aln.pdb seqfilt.* pdb95.bin alntoprof.prf alntoprof.bin \
           readprofile.pir buildprofile.prf buildprofile.ali \
	   1fas.sol 1fas.psa 1fas.ngh 1fas.dih best* final.xplor first.xplor \
           ppscan.log T3lzt-ppscan* 1lzdA.chn overhang-1.* make_pssmdb.log \
           profiles.pssm *.MRC 1cuk_init_1.pdb 1cuk_1_1.pdb \
           *.psa *.sol tmp.pdb sequence.* c-example extended-chain.pdb \
           alpha-helix.pdb fer*.ali 2abx-*.pdb

DISTFILES=1aan.seg 1fdn.swap.atm \
	  add_restraint.py align2d.py align3d.py align_consensus.py \
	  align.py all_hydrogen.py build_model.py build_profile.py \
	  color_aln_model.py compare_alignments.py define_symmetry.py \
	  describe.py edit_alignment.py ednf2.pap \
	  expand_alignment.py id_table.py iupac_model.py make_chains.py \
	  Makefile make_region.py make_restraints.py make_schedule.py \
	  make_topology_model.py malign3d.py malign.py mutate_selection.py \
	  orient_model.py patch_disulfides.py patch.py \
          pdb95.fsa c-example.c pick_atoms.py \
	  pick_hot_atoms.py pick_restraints.py randomize_xyz.py \
	  read_alignment.py group_restraints.py read_model.py \
	  read_sequence_db.py reindex_restraints.py convert_sequence_db.py \
	  rename_segments.py reorder_atoms.py \
	  rotate_dihedrals.py rotate_model.py secondary_structure.py \
	  segment_matching.py sequence_search.py shuff_seqlib_500.bin \
	  seqfilter.py join_chains.py \
	  shuff_seqlib_500.fsa superpose.py swap_atoms_in_res.py \
	  toxin.ali transfer_res_numb.py transfer_xyz.py sequences.pir \
	  very-short-for-test.cod write_data.py write_pdb_xref.py \
	  em_grid_search.py 1cuk-a2.pdb 1cuk-a2.mrc \
          aln_to_prof.py toxin.prf read_profile.py \
          ppscan.py 1hwa-uniprot90.prf 1v7sA-uniprot90.prf \
          T3lzt-uniprot90.prf 1i56A-uniprot90.prf 3lzt-uniprot90.prf \
          profiles.list overhang.ali make_pssmdb.py saxs.py build_sequence.py \
          aln_append_model.py alignment_read_one.py
INSTDIR=${MODINSTALL}/examples/commands
TESTS="mod add_restraint.py; mod align2d.py; mod align3d.py; \
       mod aln_append_model.py; mod alignment_read_one.py; \
       mod align_consensus.py; mod align.py; mod all_hydrogen.py; \
       mod aln_to_prof.py; mod build_model.py; mod build_profile.py; \
       mod color_aln_model.py; mod compare_alignments.py; \
       mod define_symmetry.py; mod describe.py; \
       mod edit_alignment.py; mod em_grid_search.py; \
       mod expand_alignment.py; mod id_table.py; \
       mod iupac_model.py; mod make_chains.py; mod make_region.py; \
       mod make_restraints.py; mod make_schedule.py; \
       mod make_topology_model.py; mod malign3d.py; mod malign.py; \
       mod mutate_selection.py; mod orient_model.py; \
       mod patch_disulfides.py; mod patch.py; mod pick_atoms.py; \
       mod pick_hot_atoms.py; mod pick_restraints.py; mod ppscan.py; \
       mod randomize_xyz.py; mod read_profile.py; \
       mod read_alignment.py; mod group_restraints.py; mod read_model.py; \
       mod read_sequence_db.py; mod convert_sequence_db.py; \
       mod reindex_restraints.py; mod rename_segments.py; \
       mod reorder_atoms.py; mod rotate_dihedrals.py; mod rotate_model.py; \
       mod secondary_structure.py; mod segment_matching.py; mod seqfilter.py; \
       mod superpose.py; mod swap_atoms_in_res.py; \
       mod transfer_res_numb.py; mod transfer_xyz.py; \
       mod write_data.py; mod write_pdb_xref.py; mod make_pssmdb.py; \
       mod build_sequence.py; mod join_chains.py"
INCDIR=..

include ${INCDIR}/Makefile.include
