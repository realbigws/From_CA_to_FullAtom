struct mod_saxsdata {
%immutable;
  int ns;
  int maxs;
  int nr;
  int natomtyp;
  int nscatts;
  double s_max;
  double s_min;
  double dr;
  struct mod_double1_array intensity;
  struct mod_double1_array int_exp;
  struct mod_double1_array sigma_exp;
  struct mod_double1_array s;
  struct mod_double1_array p_r;
  struct mod_double1_array p_r_exp;
  struct mod_double1_array r_exp;
  struct mod_double1_array p_r_resamp;
  struct mod_double1_array p_r_sig;
%mutable;
  double s_low;
  double s_hi;
  double s_hybrid;
  double c;
  double rolloff;
  double bfac;
  double offset;
  double chi_sq;
  double rho_solv;
  double dr_exp;
  int nr_exp;
};
