/* Return Fortran LOGICAL (C gboolean) variables as Python bools */

#ifdef SWIGPYTHON
%typemap(out) gboolean {
  $result = PyBool_FromLong($1);
}
#endif
