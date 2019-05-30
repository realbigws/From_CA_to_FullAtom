"""Functions to calculate percentage sequence identities"""

import sys

def id_table(aln, matrix_file):
    """Calculate percentage sequence identities"""
    coder, mat = get_seqid_matrix(aln)

    print """
Sequence identity comparison (ID_TABLE):

   Diagonal       ... number of residues;
   Upper triangle ... number of identical residues;
   Lower triangle ... % sequence identity, id/min(length).
"""

    write_seqid_matrix(sys.stdout, coder, mat)

    seqid_to_dissimilarity(mat)

    fh = open(matrix_file, "w")
    write_phylip_matrix(fh, coder, mat)


def get_seqid_matrix(aln):
    """Get a matrix of sequence identities and equivalent positions"""
    nseq = len(aln)
    mat = [[0] * nseq for n in range(nseq)]
    coder = []
    for (n1, seq1) in enumerate(aln):
        coder.append("%-6s@%3.1f" % (seq1.code[:6], abs(seq1.resolution)))
        mat[n1][n1] = len(seq1.residues)
        for (n2, seq2) in enumerate(aln):
            if n1 < n2:
                mat[n1][n2] = seq1.get_num_equiv(seq2)
                mat[n2][n1] = seq1.get_sequence_identity(seq2)
    return coder, mat

def write_seqid_matrix(fh, coder, mat):
    """Write a sequence identity matrix to a file"""
    nseq = len(mat)
    print " " * 8, "".join([ s[:8] for s in coder ])
    for n1 in range(nseq):
        print >> fh, coder[n1][:8], \
                     "".join(["%8d" % (mat[n1][n2]+0.5) for n2 in range(nseq)])

def seqid_to_dissimilarity(mat):
    """Convert a sequence identity matrix into one of dissimilarities"""
    nseq = len(mat)
    for n1 in range(nseq):
        mat[n1][n1] = 0
        for n2 in range(n1+1, nseq):
            mat[n2][n1] = 100.0 - mat[n2][n1]
            mat[n1][n2] = mat[n2][n1]

def write_phylip_matrix(fh, coder, mat):
    """Write a sequence dissimilarity matrix to a file in PHYLIP format"""
    nseq = len(mat)
    print >> fh, "%4d" % nseq
    for n1 in range(nseq):
        print >> fh, "%-10s" % coder[n1][:10] + \
                     "".join(["%4d" % (mat[n1][n2]+0.5) for n2 in range(nseq)])
