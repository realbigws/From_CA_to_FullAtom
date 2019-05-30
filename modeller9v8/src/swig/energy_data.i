%inline %{

/** Set an energy_data object's members */
static void mod_energy_data_set(struct mod_energy_data *ene,
                                float contact_shell, float relative_dielectric,
                                float radii_factor, float sphere_stdv,
                                float update_dynamic, int nonbonded_sel_atoms,
                                int nlogn_use, gboolean covalent_cys,
                                gboolean dynamic_pairs, gboolean dynamic_sphere,
                                gboolean dynamic_coulomb,
                                gboolean dynamic_lennard,
                                gboolean dynamic_modeller,
                                gboolean dynamic_access,
                                const float lennard_jones_switch[2],
                                const float coulomb_switch[2],
                                const gboolean excl_local[4])
{
  int i;
  ene->contact_shell = contact_shell;
  ene->relative_dielectric = relative_dielectric;
  ene->radii_factor = radii_factor;
  ene->sphere_stdev = sphere_stdv;
  ene->update_dynamic = update_dynamic;
  ene->nonbonded_sel_atoms = nonbonded_sel_atoms;
  ene->nlogn_use = nlogn_use;
  ene->covalent_cys = covalent_cys;
  ene->dynamic_pairs = dynamic_pairs;
  ene->dynamic_sphere = dynamic_sphere;
  ene->dynamic_coulomb = dynamic_coulomb;
  ene->dynamic_lennard = dynamic_lennard;
  ene->dynamic_modeller = dynamic_modeller;
  ene->dynamic_access = dynamic_access;
  for (i = 0; i < 2; ++i) {
    ene->lennard_jones_switch[i] = lennard_jones_switch[i];
    ene->coulomb_switch[i] = coulomb_switch[i];
  }
  for (i = 0; i < 4; ++i) {
    ene->excl_local[i] = excl_local[i];
  }
}

%}
