from modeller.sequence_db import sequence_db
from modeller.alignment import alignment

def sequence_srch(env, sequence, segfile, chains_list,
                  toplib='${LIB}/top_heav.lib', search_randomizations=0,
                  signif_cutoff=(4.0, 5.0)):
    """For a given target sequence, find all template chains in PDB.
       - align all of them together and write the alignment to a file.
       - calculate identity matrix and write it to a file.
       - calculate a dendrogram and write it to the .log file."""

    penalties = (-600, -400) # This should be fine for significance tests with
                             # default matrix as1.sim.mat
    overhang=20
    seqfile = '$(LIB)/pdball.pir'

    sdb = sequence_db(env, seq_database_file=seqfile, chains_list=chains_list,
                      seq_database_format='PIR')
    aln = alignment(env, file=segfile, align_codes=sequence)
    sdb.search(aln, search_top_list=30, off_diagonal=9999,
               search_group_list = '$(LIB)/pdb_40.grp',
               seq_database_file=seqfile, gap_penalties_1d=penalties,
               output='SHORT', signif_cutoff=signif_cutoff,
               search_randomizations=search_randomizations)

    if len(aln) > 0:
        aln.write(file='alignment.tmp')
        aln.malign(overhang=overhang, off_diagonal=150,
                   gap_penalties_1d=penalties)
        aln.malign3d(overhang=overhang, gap_penalties_3d=(0, 3),
                     off_diagonal=150)
        blk = len(aln)
        aln.append(file=segfile, align_codes=sequence)

        # Only needed for ALIGN2D's PSA run
        env.libs.topology.read(file=toplib)

        aln.align2d(align_block=blk, gap_penalties_1d=(-450, 0),
                    gap_penalties_2d=(0.35, 1.2, 0.9, 1.2, 0.6, 8.6, 1.2, 0, 0),
                    max_gap_length=50, off_diagonal=150, overhang=overhang)
        aln.write(file=sequence + '.ali')
        aln.write(file=sequence + '.pap', alignment_format='PAP',
                  alignment_features='HELIX BETA ACCESSIBILITY ' + \
                                     'STRAIGHTNESS CONSERVATION INDICES')
        mat = sequence + '.mat'
        aln.id_table(matrix_file=mat)
        env.dendrogram(matrix_file=mat)
