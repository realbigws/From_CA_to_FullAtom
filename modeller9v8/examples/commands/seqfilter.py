from modeller import *

log.verbose()
env = environ()

sdb = sequence_db(env, seq_database_file='sequences.pir',
                  seq_database_format='PIR',
                  chains_list='ALL', minmax_db_seq_len=[30, 3000],
                  clean_sequences=True)

sdb.filter(rr_file='${LIB}/id.sim.mat', gap_penalties_1d=[-3000, -1000],
           max_diff_res=30, seqid_cut=95, output_grp_file='seqfilt.grp',
           output_cod_file='seqfilt.cod')
