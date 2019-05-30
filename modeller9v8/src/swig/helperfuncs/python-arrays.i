/* Convert a Python list to a C array */
%define PYTOARRAY(type, name, chkfunc, convfunc, errmsg)
%{
static type *python_to_ ## name ## _array(PyObject *input, int size, int *sizevar, type *destlist, char *varname) {
  int i, seqlen;
  Py_ssize_t pyseqlen;
  type *result;
  /* Return a 1-element list if the user supplied a lone value */
  if (sizevar && chkfunc(input)) {
    result = malloc(sizeof(type));
    *sizevar = 1;
    result[0] = (type)convfunc(input);
    return result;
  }
  if (!PySequence_Check(input) || PyString_Check(input)) {
    PyErr_Format(PyExc_ValueError, "Expected a sequence for %s", varname);
    return NULL;
  }
  /* Get sequence length, and make sure it fits into an int (not necessarily
     the case with Python 2.5, which uses 64-bit lengths on 64-bit platforms) */
  pyseqlen = PySequence_Length(input);
  if (pyseqlen > INT_MAX) {
    PyErr_Format(PyExc_ValueError,
                 "Length of sequence for %s exceeds maximum size", varname);
    return NULL;
  }
  seqlen = (int)pyseqlen;
  if (sizevar) {
    *sizevar = size = seqlen;
  } else if (seqlen != size) {
    PyErr_Format(PyExc_ValueError,
                 "Expected a sequence of length %d for %s; got %d", size,
                 varname, seqlen);
    return NULL;
  }
  if (destlist) {
    result = destlist;
  } else if (size == 0) {
    result = malloc(sizeof(type));
  } else {
    result = malloc(size * sizeof(type));
  }
  for (i = 0; i < size; i++) {
    PyObject *o = PySequence_GetItem(input,i);
    if (chkfunc(o)) {
      result[i] = (type)convfunc(o);
      Py_DECREF(o);
    } else {
      Py_XDECREF(o);
      PyErr_Format(PyExc_ValueError, errmsg, varname, i);
      if (!destlist) {
        free(result);
      }
      return NULL;
    }
  }
  return result;
}
%}
%enddef

PYTOARRAY(int, int, PyInt_Check, PyInt_AsLong, "%s[%d] must be an integer")
PYTOARRAY(float, float, PyNumber_Check, PyFloat_AsDouble, "%s[%d] must be a number")
PYTOARRAY(gboolean, gboolean, PyBool_Check, PyInt_AsLong, "%s[%d] must be a bool")


/* Convert a C array into a Python list */
%define PYFROMARRAY(type, convfunc)
%{
static PyObject *python_from_ ## type ## _array(type *inarr, int len)
{
  PyObject *tuple = PyTuple_New(len);
  if (tuple) {
    int i;
    for (i = 0; i < len; i++) {
      PyObject *obj = convfunc(inarr[i]);
      if (!obj) {
        Py_DECREF(tuple);
        tuple = NULL;
        break;
      }
      PyTuple_SetItem(tuple, i, obj);
    }
  }
  return tuple;
}
%}
%enddef

PYFROMARRAY(int, PyInt_FromLong)
PYFROMARRAY(float, PyFloat_FromDouble)

%{
/* Convert a Python list into a Fortran string array */
static char *python_to_f_str_array(PyObject *input, int size, int *sizevar,
                                   int *strlenvar, char *varname) {
  int i, memsize, seqlen;
  Py_ssize_t pyseqlen;
  char *result;
  *strlenvar = 1;
  /* Return a one-element list if the user provided a lone string */
  if (sizevar && PyString_Check(input)) {
    int len;
    char *pt = PyString_AsString(input);
    len = strlen(pt);
    *strlenvar = len;
    *sizevar = 1;
    if (len == 0) {
      result = malloc(1);
    } else {
      result = malloc(len);
      memcpy(result, pt, len);
    }
    return result;
  }
  if (!PySequence_Check(input) || PyString_Check(input)) {
    PyErr_Format(PyExc_ValueError, "Expected a sequence for %s", varname);
    return NULL;
  }
  /* Get sequence length, and make sure it fits into an int (not necessarily
     the case with Python 2.5, which uses 64-bit lengths on 64-bit platforms) */  pyseqlen = PySequence_Length(input);
  if (pyseqlen > INT_MAX) {
    PyErr_Format(PyExc_ValueError,
                 "Length of sequence for %s exceeds maximum size", varname);
    return NULL;
  }
  seqlen = (int)pyseqlen;
  if (sizevar) {
    *sizevar = size = seqlen;
  } else if (seqlen != size) {
    PyErr_Format(PyExc_ValueError,
                 "Expected a sequence of length %d for %s; got %d", size,
                 varname, seqlen);
    return NULL;
  }
  for (i = 0; i < size; i++) {
    PyObject *o = PySequence_GetItem(input,i);
    if (PyString_Check(o)) {
      int len;
      char *pt = PyString_AsString(o);
      len = strlen(pt);
      if (len > *strlenvar) *strlenvar = len;
      Py_DECREF(o);
    } else {
      Py_XDECREF(o);
      PyErr_Format(PyExc_ValueError, "%s[%d] must be a string", varname, i);
      return NULL;
    }
  }
  memsize = *strlenvar;
  if (memsize == 0) {
    memsize = 1;
  }
  if (size > 0)  {
    memsize *= size;
  }
  result = malloc(memsize);
  memset(result, ' ', memsize);
  for (i = 0; i < size; i++) {
    PyObject *o = PySequence_GetItem(input,i);
    char *pt = PyString_AsString(o);
    memcpy(result + i * (*strlenvar), pt, strlen(pt));
    Py_DECREF(o);
  }
  return result;
}
%}

%{
/* Convert a Python list into a C string array (NULL-terminated char **) */
static char **python_to_c_str_array(PyObject *input, char *varname) {
  int i, seqlen;
  Py_ssize_t pyseqlen;
  char **result;
  /* Return a one-element list if the user provided a lone string */
  if (PyString_Check(input)) {
    result = g_malloc0(2 * sizeof(char *));
    result[0] = g_strdup(PyString_AsString(input));
    return result;
  } else if (!PySequence_Check(input)) {
    PyErr_Format(PyExc_ValueError, "Expected a sequence for %s", varname);
    return NULL;
  }
  /* Get sequence length, and make sure it fits into an int (not necessarily
     the case with Python 2.5, which uses 64-bit lengths on 64-bit platforms) */  pyseqlen = PySequence_Length(input);
  if (pyseqlen > INT_MAX) {
    PyErr_Format(PyExc_ValueError,
                 "Length of sequence for %s exceeds maximum size", varname);
    return NULL;
  }
  seqlen = (int)pyseqlen;
  result = g_malloc0((seqlen + 1) * sizeof(char *));
  for (i = 0; i < seqlen; i++) {
    PyObject *o = PySequence_GetItem(input,i);
    if (PyString_Check(o)) {
      result[i] = g_strdup(PyString_AsString(o));
      Py_DECREF(o);
    } else {
      Py_XDECREF(o);
      PyErr_Format(PyExc_ValueError, "%s[%d] must be a string", varname, i);
      g_strfreev(result);
      return NULL;
    }
  }
  return result;
}
%}
