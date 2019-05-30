struct mod_optimizer {
%immutable;
  int step;
  int funcs;
  float init_e;
  float current_e;
  float shiftavr;
  float shiftmax;
};

struct mod_cg_optimizer {
%immutable;
  struct mod_optimizer opt;
  float gsq;
};

struct mod_md_optimizer {
%immutable;
  struct mod_optimizer opt;
  float aver_e;
  float max_e;
  float min_e;
  float stdev_e;
  float smdavr;
  float kinetic_e;
  float kinetic_temp;
  int imdmax;
  int imdmin;
%mutable;
  float temperature;
};

struct mod_state_optimizer {
%immutable;
  struct mod_optimizer opt;
};

struct mod_qn_optimizer {
%immutable;
  struct mod_optimizer opt;
};
