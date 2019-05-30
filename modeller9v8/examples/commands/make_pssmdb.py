from modeller import *

env = environ()

env.make_pssmdb(profile_list_file = 'profiles.list',
                matrix_offset     = -450,
                rr_file           = '${LIB}/blosum62.sim.mat',
                pssmdb_name       = 'profiles.pssm',
                profile_format    = 'TEXT',
                pssm_weights_type = 'HH1')
