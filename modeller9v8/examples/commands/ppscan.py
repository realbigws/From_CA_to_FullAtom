# Example for: profile.scan()

from modeller import *

env = environ()

# First create a database of PSSMs
env.make_pssmdb(profile_list_file = 'profiles.list',
                matrix_offset     = -450,
                rr_file           = '${LIB}/blosum62.sim.mat',
                pssmdb_name       = 'profiles.pssm',
                profile_format    = 'TEXT',
                pssm_weights_type = 'HH1')

# Read in the target profile
prf = profile(env, file='T3lzt-uniprot90.prf', profile_format='TEXT')

# Read the PSSM database
psm = pssmdb(env, pssmdb_name = 'profiles.pssm', pssmdb_format = 'text')

# Scan against all profiles in the 'profiles.list' file
# The score_statistics flag is set to false since there are not
# enough database profiles to calculate statistics.
prf.scan(profile_list_file = 'profiles.list',
         psm               = psm,
         matrix_offset     = -450,
         ccmatrix_offset   = -100,
         rr_file           = '${LIB}/blosum62.sim.mat',
         gap_penalties_1d  = (-700, -70),
         score_statistics  = False,
         output_alignments = True,
         output_score_file = None,
         profile_format    = 'TEXT',
         max_aln_evalue    = 1,
         aln_base_filename = 'T3lzt-ppscan',
         pssm_weights_type = 'HH1',
         summary_file      = 'T3lzt-ppscan.sum')
