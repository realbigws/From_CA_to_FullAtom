from modeller import *
log.verbose()
env = environ()

sdb = sequence_db(env)
sdb.convert(seq_database_file='pdb95.fsa', seq_database_format='FASTA',
            chains_list='ALL', minmax_db_seq_len=[1, 40000],
            clean_sequences=True, outfile='pdb95.bin')
