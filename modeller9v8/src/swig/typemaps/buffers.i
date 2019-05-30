/* Convert Python strings to/from C fixed-size buffers */


#ifdef SWIGPYTHON
%typemap(in) (const char *buffer, size_t bufsize) (Py_ssize_t temp) {
  if (PyString_AsStringAndSize($input, &$1, &temp) == -1) {
    SWIG_fail;
  } else {
    $2 = (size_t)temp;
  }
}

%typemap(in,fragment=SWIG_AsVal_frag(size_t)) (char *buffer, size_t bufsize) (int res) {
  res = SWIG_AsVal(size_t)($input, &$2);
  if (!SWIG_IsOK(res)) {
    %argument_fail(res, "(size_t bufsize)", $symname, $argnum);
  }
  $1 = g_malloc($2);
}

%typemap(freearg,match="in") (char *buffer, size_t bufsize) {
  g_free($1);
}

%typemap(argout) (char *buffer, size_t bufsize) {
  %append_output(PyString_FromStringAndSize($1, (Py_ssize_t)$2));
}
#endif
