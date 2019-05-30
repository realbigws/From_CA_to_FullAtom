from modeller.alignment import alignment

def principal_components(env, family, cluster_cut):
    """Does principal components on family"""
    aln = alignment(env, file=family+'.ali')

    mat = family + '.mat'
    aln.id_table(matrix_file=mat)

    env.principal_components(matrix_file=mat, file=family+'.dat')

    env.dendrogram(matrix_file=mat, cluster_cut=cluster_cut)
