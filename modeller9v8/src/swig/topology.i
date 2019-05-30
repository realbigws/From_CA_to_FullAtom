struct mod_topology {
%immutable;
  int nresrt;
  int nrtl;
  int ntopmod;
  gboolean autang;
  gboolean autdih;
  struct mod_float1_array amrt;
%mutable;
  int submodel;
};
