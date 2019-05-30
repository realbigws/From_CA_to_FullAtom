struct mod_alignment {
%immutable;
  int nseq;
  int naln;
  int ncomment;
  struct mod_int2_array ialn;
  struct mod_int2_array invaln;
  struct mod_float2_array prof;
};
