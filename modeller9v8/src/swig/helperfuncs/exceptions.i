%{
#ifdef SWIGPYTHON
/* Generic Modeller error */
static PyObject *moderror;
/* File format error */
static PyObject *file_format_error;
/* Statistics error */
static PyObject *statistics_error;
/* Sequence mismatch error */
static PyObject *sequence_mismatch_error;
#endif

/* Raise a Python or Perl exception if Modeller returned an error */
static int check_for_error(void)
{
  GError *err;
  static const char *interr = "INTERNAL error: Error code returned, but no "
                              "error information set. This is a bug in "
                              "Modeller. Please report it to the Modeller "
                              "developers, including any input files "
                              "necessary to reproduce the problem.";
  if ((err = mod_error_get())) {
#ifdef SWIGPYTHON
    if (err->domain == MOD_ERROR) {
      switch(err->code) {
      case MOD_ERROR_ZERODIV:
        PyErr_SetString(PyExc_ZeroDivisionError, err->message); break;
      case MOD_ERROR_IO:
        PyErr_SetString(PyExc_IOError, err->message); break;
      case MOD_ERROR_MEMORY:
        PyErr_SetString(PyExc_MemoryError, err->message); break;
      case MOD_ERROR_EOF:
        PyErr_SetString(PyExc_EOFError, err->message); break;
      case MOD_ERROR_TYPE:
        PyErr_SetString(PyExc_TypeError, err->message); break;
      case MOD_ERROR_NOTIMP:
        PyErr_SetString(PyExc_NotImplementedError, err->message); break;
      case MOD_ERROR_INDEX:
        PyErr_SetString(PyExc_IndexError, err->message); break;
      case MOD_ERROR_VALUE:
        PyErr_SetString(PyExc_ValueError, err->message); break;
      case MOD_ERROR_FILE_FORMAT:
        PyErr_SetString(file_format_error, err->message); break;
      case MOD_ERROR_OVERFLOW:
        PyErr_SetString(PyExc_OverflowError, err->message); break;
      case MOD_ERROR_STATISTICS:
        PyErr_SetString(statistics_error, err->message); break;
      case MOD_ERROR_SEQUENCE_MISMATCH:
        PyErr_SetString(sequence_mismatch_error, err->message); break;
      case MOD_ERROR_FAILED:
        PyErr_SetString(moderror, err->message);
      }
    } else {
      PyErr_SetString(moderror, err->message);
    }
#endif
#ifdef SWIGPERL
    croak(err->message);
#endif
    g_error_free(err);
    mod_error_clear();
    return 1;
#ifdef SWIGPYTHON
  /* Propagate a Python exception already raised (e.g. by a callback) */
  } else if (PyErr_Occurred()) {
    return 1;
#endif
  /* An error code was returned, but no error is set - this is also an error! */
  } else {
#ifdef SWIGPYTHON
    PyErr_SetString(moderror, interr);
#endif
#ifdef SWIGPERL
    croak(interr);
#endif
    return 1;
  }
}
%}
