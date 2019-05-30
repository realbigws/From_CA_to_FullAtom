struct mod_density {
%immutable;
  float vox_size;
  int nx, ny, nz;
  float px, py, pz;
  float grid_sqr_sum;
  int filter_type;
  int function_type;
  float cc;
  double norm_factor;
  struct mod_float3_array grid;
%mutable;
  float resolution;
};
