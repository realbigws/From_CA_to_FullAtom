struct mod_sequence {
%mutable;
  char *name, *source, *prottyp;
  float resol;
  float rfactr;
%immutable;
  int nres;
  int nseg;
  struct mod_int1_array irestyp;
  struct mod_int1_array iress1;
  struct mod_int1_array iress2;
};
