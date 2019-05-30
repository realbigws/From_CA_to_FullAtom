from modeller.alignment import alignment
from modeller.selection import selection
from modeller.model import model

def fit(env, model, code, model2, code2, alnfile, model2_fit):
    """Superposes model2 on model, and writes out a file with model2 superposed
       on model."""
    m1 = model(env, file=model)
    m2 = model(env, file=model2)

    aln = alignment(env, file=alnfile, align_codes=(code, code2))
    atmsel = selection(m1).only_atom_types('CA')

    atmsel.superpose(m2, aln)

    m2.write(file=model2_fit)
