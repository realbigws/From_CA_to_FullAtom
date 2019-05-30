/* Conversions between Python strings and C GString objects */

%typemap(in) GString * (char *temp) {
  $1 = NULL;
  if (PyString_Check($input)) {
    temp = PyString_AsString($input);
  } else {
    PyErr_Format(PyExc_TypeError, "Expected a string for argument number %d",
                 $argnum);
    SWIG_fail;
  }
}

%typemap(memberin) GString * {
  g_string_assign($1, temp2);
}

%typemap(out) GString * {
  $result = PyString_FromString($1->str);
}
