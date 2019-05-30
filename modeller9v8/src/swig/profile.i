struct mod_profile {
%immutable;
  int nseq, npos;
  struct mod_int1_array iter;
  struct mod_int1_array neqv;
  struct mod_int1_array nres;
  struct mod_float1_array fid;
  struct mod_double1_array evalue;
  struct mod_int2_array sprofile;
  char *filename;
};
