/* Convert Modeller errors (ierr returned != 0) to Python or Perl exceptions */

%typemap(in, numinputs=0) int *ierr (int temp) {
  $1 = &temp;
}

%typemap(argout) int *ierr {
  if (*$1 != 0 && check_for_error()) {
#ifdef SWIGPYTHON
    Py_DECREF(resultobj);
#endif
    SWIG_fail;
  }
}

/* When the initial value is used, initialize it with the current
   error status (if known). */
%apply int *ierr { int *ierr_inout };
%typemap(in, numinputs=0) int *ierr_inout (int temp) {
#ifdef SWIGPYTHON
  temp = (PyErr_Occurred() ? 1 : 0);
#else
  temp = 0;
#endif
  $1 = &temp;
}
