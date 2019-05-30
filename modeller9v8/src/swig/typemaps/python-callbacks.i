/* Convert Python callable objects into C callback routines, so that
   we can call Python routines from within Modeller C or Fortran code. */

/* Evaluate a user-defined feature. */
%typemap(in) (cb_feval evalfunc, void *evaldata) {
  if (set_python_callback((void **)&$1, &$2, $input, python_cb_feval)) {
    SWIG_fail;
  }
}

/* Evaluate the derivatives of a user-defined feature. */
%typemap(in) (cb_fderiv derivfunc, void *derivdata) {
  if (set_python_callback((void **)&$1, &$2, $input, python_cb_fderiv)) {
    SWIG_fail;
  }
}

/* Identify whether a user-defined feature is an 'angle' */
%typemap(in) (cb_fangle anglefunc, void *angledata) {
  if (set_python_callback((void **)&$1, &$2, $input, python_cb_fangle)) {
    SWIG_fail;
  }
}

/* Evaluate a user-defined energy term */
%typemap(in) (cb_teval evalfunc, cb_free freefunc, void *evaldata) {
  if (set_python_callback_free((void **)&$1, &$2, &$3, $input,
                               python_cb_teval)) {
    SWIG_fail;
  }
}

/* Return the user-defined local similarity. */
%typemap(in) (cb_simloc evalfunc, cb_free freefunc, void *evaldata) {
  if (set_python_callback_free((void **)&$1, &$2, &$3, $input,
                               python_cb_simloc)) {
    SWIG_fail;
  }
}

/* Evaluate a user-defined functional form. */
%typemap(in) (cb_reval evalfunc, void *evaldata) {
  if (set_python_callback((void **)&$1, &$2, $input, python_cb_reval)) {
    SWIG_fail;
  }
}

/* Evaluate the violation of a user-defined restraint form */
%typemap(in) (cb_vmin vminfunc, void *vmindata) {
  if (set_python_callback((void **)&$1, &$2, $input, python_cb_violation)) {
    SWIG_fail;
  }
}
%apply (cb_vmin vminfunc, void *vmindata) { (cb_vheavy vheavyfunc, void *vheavydata) };
%apply (cb_vmin vminfunc, void *vmindata) { (cb_rvmin rvminfunc, void *rvmindata) };
%apply (cb_vmin vminfunc, void *vmindata) { (cb_rvheavy rvheavyfunc, void *rvheavydata) };

/* Evaluate the mean of a user-defined restraint form */
%typemap(in) (cb_rsrmin rsrminfunc, void *rsrmindata) {
  if (set_python_callback((void **)&$1, &$2, $input, python_cb_mean)) {
    SWIG_fail;
  }
}
%apply (cb_rsrmin rsrminfunc, void *rsrmindata) { (cb_rsrheavy rsrheavyfunc, void *rsrheavydata) };

/* Return dervx/y/z arrays (e.g. from featder). Note that we hardcode arg3,
   which is n_indatm */
%typemap(in,numinputs=0) float dervx[] {
}
%typemap(check) float dervx[] {
  $1 = malloc(arg3 * sizeof(float));
}

%typemap(argout,fragment="t_output_helper") float dervx[] {
  PyObject *o = python_from_float_array($1, arg3);
  $result = t_output_helper($result, o);
}
%apply float dervx[] { float dervy[] };
%apply float dervx[] { float dervz[] };

/* Carry out an optimizer action */
%typemap(in) (cb_action actionfunc, cb_free freefunc, void *actiondata) {
  if (set_python_callback_free((void **)&$1, &$2, &$3, $input,
                               python_cb_action)) {
    SWIG_fail;
  }
}

/* Log an output message */
%typemap(in) (cb_log logfunc, cb_free freefunc, void *logdata) {
  if (set_python_callback_free((void **)&$1, &$2, &$3, $input,
                               python_cb_log)) {
    SWIG_fail;
  }
}
