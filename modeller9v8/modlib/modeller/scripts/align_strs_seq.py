from modeller.alignment import alignment

def align_strs_seq(env, segfile, alnfile, knowns, sequence, matrix_file,
                   overhang=0, write_fit=False):
    """Align a single sequence with several structures"""

    # Read the sequences of structures from the specs in SEGFILE:
    aln = alignment(env, file=segfile, align_codes=knowns)

    # Only align structures if there's more than one:
    if len(aln) > 1:
        # do a multiple sequence alignment of structures:
        aln.malign(gap_penalties_1d=(-600, -400), overhang=overhang)

        # do a multiple structural alignment of structures:
        aln.malign3d(gap_penalties_3d=(0.0, 2.0), fit_atoms='CA',
                     overhang=overhang, write_fit=write_fit)

    # remember the number of structures
    align_block = len(aln)

    # add the sequence of the unknown to the sequence/alignment arrays
    # containing the aligned structures:
    aln.append(file=segfile, align_codes=sequence)

    # align the last sequence with the fixed alignment of structures:
    aln.align(align_block=align_block, gap_penalties_1d=(-600, -400),
              overhang=overhang)

    # write the alignment to a file ALNFILE
    aln.write(file=alnfile)

    # do some sequence comparisons:
    aln.id_table(matrix_file)

    return aln
