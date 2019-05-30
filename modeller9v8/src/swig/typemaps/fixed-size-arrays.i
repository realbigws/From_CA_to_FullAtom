/* Convert fixed-length Perl/Python lists to/from C arrays */

#ifdef SWIGPERL
%typemap(in) int[ANY] {
    $1 = perl_to_int_array($input, $1_dim0, NULL, "$1_name");
}

%typemap(in) float[ANY] {
    $1 = perl_to_float_array($input, $1_dim0, NULL, "$1_name");
}
#endif

#ifdef SWIGPYTHON
%typemap(argout,fragment="t_output_helper") float[ANY] {
    PyObject *o = python_from_float_array($1, $1_dim0);
    $result = t_output_helper($result, o);
}
%typemap(in,numinputs=0) float[ANY] {
  $1 = malloc(sizeof(float) * $1_dim0);
}
%typemap(argout,fragment="t_output_helper") int[ANY] {
    PyObject *o = python_from_int_array($1, $1_dim0);
    $result = t_output_helper($result, o);
}
%typemap(in,numinputs=0) int[ANY] {
  $1 = malloc(sizeof(int) * $1_dim0);
}

%typemap(in) const int[ANY] {
    $1 = python_to_int_array($input, $1_dim0, NULL, NULL, "$1_name");
    if (!$1) SWIG_fail;
}
%typemap(argout) const int[ANY] {
}

%typemap(in) const float[ANY] {
    $1 = python_to_float_array($input, $1_dim0, NULL, NULL, "$1_name");
    if (!$1) SWIG_fail;
}
%typemap(argout) const float[ANY] {
}

%typemap(in) const gboolean[ANY] {
    $1 = python_to_gboolean_array($input, $1_dim0, NULL, NULL, "$1_name");
    if (!$1) SWIG_fail;
}
%typemap(argout) const gboolean[ANY] {
}
#endif

/* Free the C arrays after calling the Modeller function */
%typemap(freearg) int[ANY] {
    if ($1) free($1);
}
%typemap(freearg) float[ANY] = int[ANY];
%typemap(freearg) gboolean[ANY] = int[ANY];

/* Initialize the C array pointer to avoid crashes in freearg (SWIG does not
   init this pointer for us) */
%typemap(arginit) int[ANY] {
  $1 = NULL;
}
%typemap(arginit) float[ANY] = int[ANY];
%typemap(arginit) gboolean[ANY] = int[ANY];

