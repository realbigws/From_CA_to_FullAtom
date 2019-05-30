/* Generic typemaps for handling input and output of variable-sized arrays.
   The actual typemaps are %apply'd in swig-arrays.i */

#ifdef SWIGPYTHON
%typemap(in) (const char STRARRAY[], int STRLEN) {
    $1 = python_to_f_str_array($input, $1_dim0, NULL, &$2, "$1_name");
    if (!$1) SWIG_fail;
}
#endif
#ifdef SWIGPERL
%typemap(in) (const char STRARRAY[], int STRLEN) {
    $1 = perl_to_f_str_array($input, $1_dim0, NULL, &$2, "$1_name");
}
#endif
%typemap(arginit) (const char STRARRAY[], int STRLEN) {
  $1 = NULL;
}
%typemap(freearg) (const char STRARRAY[], int STRLEN) {
    if ($1) free($1);
}

#ifdef SWIGPYTHON
%typemap(in) (const float VARARRAY[], int N_VARARRAY) {
    $1 = python_to_float_array($input, 0, &$2, NULL, "$1_name");
    if (!$1) SWIG_fail;
}
#endif
#ifdef SWIGPERL
%typemap(in) (const float VARARRAY[], int N_VARARRAY) {
    $1 = perl_to_float_array($input, 0, &$2, "$1_name");
}
#endif
%typemap(freearg) (const float VARARRAY[], int N_VARARRAY) {
    if ($1) free($1);
}
%typemap(arginit) (const float VARARRAY[], int N_VARARRAY) {
  $1 = NULL;
}
#ifdef SWIGPYTHON
%typemap(in) (const int VARARRAY[], int N_VARARRAY) {
    $1 = python_to_int_array($input, 0, &$2, NULL, "$1_name");
    if (!$1) SWIG_fail;
}
#endif
#ifdef SWIGPERL
%typemap(in) (const int VARARRAY[], int N_VARARRAY) {
    $1 = perl_to_int_array($input, 0, &$2, "$1_name");
}
#endif
%typemap(freearg) (const int VARARRAY[], int N_VARARRAY) {
    if ($1) free($1);
}
%typemap(arginit) (const int VARARRAY[], int N_VARARRAY) {
  $1 = NULL;
}
#ifdef SWIGPYTHON
%typemap(argout,fragment="t_output_helper") (int **OUTVARARRAY, int *N_OUTVARARRAY) {
    PyObject *o = python_from_int_array(*$1, *$2);
    $result = t_output_helper($result, o);
}
#endif
%typemap(freearg) (int **OUTVARARRAY, int *N_OUTVARARRAY) {
    if (*$1) free(*$1);
}
#ifdef SWIGPYTHON
%typemap(in,numinputs=0) (int **OUTVARARRAY, int *N_OUTVARARRAY) (int *temp1, int temp2) {
  temp1 = NULL;
  $1 = &temp1;
  $2 = &temp2;
}
%typemap(argout,fragment="t_output_helper") (float **OUTVARARRAY, int *N_OUTVARARRAY) {
    PyObject *o = python_from_float_array(*$1, *$2);
    $result = t_output_helper($result, o);
}
#endif
%typemap(freearg) (float **OUTVARARRAY, int *N_OUTVARARRAY) {
    if (*$1) free(*$1);
}
#ifdef SWIGPYTHON
%typemap(in,numinputs=0) (float **OUTVARARRAY, int *N_OUTVARARRAY) (float *temp1, int temp2) {
  temp1 = NULL;
  $1 = &temp1;
  $2 = &temp2;
}
#endif

#ifdef SWIGPYTHON
%typemap(in) (const struct mod_model *VARARRAY[], int N_VARARRAY) {
    $1 = python_to_model_array($input, 0, &$2, NULL, "$1_name");
    if (!$1) SWIG_fail;
}
#endif
%typemap(freearg) (const struct mod_model *VARARRAY[], int N_VARARRAY) {
    if ($1) free($1);
}
%typemap(arginit) (const struct mod_model *VARARRAY[], int N_VARARRAY) {
  $1 = NULL;
}

#ifdef SWIGPYTHON
%typemap(in) (const char VARSTRARRAY[], int STRLEN, int N_VARSTRARRAY) {
    $1 = python_to_f_str_array($input, 0, &$3, &$2, "$1_name");
    if (!$1) SWIG_fail;
}
#endif
%typemap(freearg) (const char VARSTRARRAY[], int STRLEN, int N_VARSTRARRAY) {
    if ($1) free($1);
}
%typemap(arginit) (const char VARSTRARRAY[], int STRLEN, int N_VARSTRARRAY) {
  $1 = NULL;
}
