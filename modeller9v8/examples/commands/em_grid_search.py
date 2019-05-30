from modeller import *

log.verbose()
env = environ()

den = density(env, file='1cuk-a2.mrc', em_density_format='MRC',
              voxel_size=1., resolution=8., em_map_size=40,
              cc_func_type='CCF', density_type='SPHERE')

den.grid_search(em_density_format='MRC', num_structures=1,
                em_pdb_name=['1cuk-a2.pdb'], chains_num=[1],
                start_type='CENTER', number_of_steps=1, angular_step_size=30.,
                temperature=0., best_docked_models=1,
                em_fit_output_file='test-cr.log')
