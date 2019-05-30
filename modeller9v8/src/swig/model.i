struct mod_model {
%mutable;
  float seq_id;
  float last_energy;
  GString *remark, *header;
%immutable;
  struct mod_coordinates cd;
  struct mod_sequence seq;
  struct mod_float1_array charge;
  struct mod_float1_array vx, vy, vz;
  struct mod_float1_array dvx, dvy, dvz;
  struct mod_int1_array iatta, iattyp;
  struct mod_pseudo_atoms psd;
  struct mod_schedule sched;
  struct mod_symmetry sym;
  struct mod_restraints rsr;
  struct mod_model_topology mtp;
};
