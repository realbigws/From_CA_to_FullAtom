/* Special handling for exceptions raised by mod_file_open */
%exception mod_file_open {
  $action
  if (!result && check_for_error()) {
    SWIG_fail;
  }
}

/* Don't wrap the mod_file structure; should be treated as opaque */
%ignore mod_file;

/* Provide our own implementation for mod_file_open_stream */
%ignore mod_file_open_stream;

%{
/* Adapter to use cStringIO read with mod_file objects */
static size_t cstr_readfunc(void *data, char *buffer, size_t bufsize, int *ierr)
{
  char *buf;
  int len;
  PyObject *obj = (PyObject *)data;
  len = PycStringIO->cread(obj, &buf, (Py_ssize_t)bufsize);
  if (PyErr_Occurred()) {
    *ierr = 1;
    return 0;
  } else {
    *ierr = 0;
    if ((size_t)len == bufsize) {
      memcpy(buffer, buf, len);
    }
  }
  return (size_t)len;
}

/* Adapter to use cStringIO write with mod_file objects */
static void cstr_writefunc(void *data, const char *buffer, size_t bufsize,
                           int *ierr)
{
  PyObject *obj = (PyObject *)data;
  PycStringIO->cwrite(obj, (char *)buffer, (Py_ssize_t)bufsize);
  *ierr = (PyErr_Occurred() != NULL);
}

/* Release the reference to the underlying cStringIO object */
static void cstr_freefunc(void *data)
{
  Py_DECREF((PyObject *)data);
}
%}

%inline %{
/* Open a mod_file object that accesses a Python cStringIO object */
static struct mod_file *mod_file_open_cstring(PyObject *obj, int *ierr)
{
  if (!PycStringIO) {
    PyErr_SetString(PyExc_ImportError, "Could not import cStringIO module");
    *ierr = 1;
    return NULL;
  } else if (!PycStringIO_OutputCheck(obj) && !PycStringIO_InputCheck(obj)) {
    PyErr_SetString(PyExc_TypeError, "Expecting a cStringIO object");
    *ierr = 1;
    return NULL;
  } else {
    *ierr = 0;
    Py_INCREF(obj);
    return mod_file_open_stream(cstr_readfunc, cstr_writefunc, cstr_freefunc,
                                obj);
  }
}  
%} /* inline */

%include "../include/mod_file.h"
