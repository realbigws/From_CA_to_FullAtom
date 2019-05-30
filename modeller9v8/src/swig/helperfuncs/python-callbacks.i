/* Callback functions to allow calling Python functions from within Modeller
   C or Fortran code. The arguments to these C functions should match the
   definitions for callback functions in user_features.h, etc., and they
   convert C values and errors to/from Python. */

%{
/* Identify whether a user-defined feature is an 'angle' */
static int python_cb_fangle(void *data, gboolean *out1)
{
  PyObject *func, *arglist, *result;

  func = (PyObject *)data;
  if (!(arglist = Py_BuildValue("()"))) {
    return 1;
  }
  result = PyEval_CallObject(func, arglist);
  Py_DECREF(arglist);
  if (result && PyBool_Check(result)) {
    *out1 = PyInt_AsLong(result);
    Py_DECREF(result);
    return 0;
  } else {
    if (result) {
      Py_DECREF(result);
      PyErr_SetString(PyExc_TypeError,
                      "Callback function should return a bool");
    }
    return 1;
  }
}

/* Evaluate a user-defined feature. */
static int python_cb_feval(void *data, const struct mod_model *model, int *in1,
                           int in2, float *out1)
{
  PyObject *func, *arglist, *result, *tuple;

  func = (PyObject *)data;
  tuple = python_from_int_array(in1, in2);
  arglist = Py_BuildValue("(OO)", model->scriptobj, tuple);
  Py_XDECREF(tuple);
  if (!arglist) {
    return 1;
  }
  result = PyEval_CallObject(func, arglist);
  Py_DECREF(arglist);
  if (result && PyNumber_Check(result)) {
    *out1 = PyFloat_AsDouble(result);
    Py_DECREF(result);
    return 0;
  } else {
    if (result) {
      Py_DECREF(result);
      PyErr_SetString(PyExc_TypeError,
                      "Callback function should return a number");
    }
    return 1;
  }
}

/* Evaluate the derivatives of a user-defined feature. */
static int python_cb_fderiv(void *data, const struct mod_model *model, int *in1,
                            int in2, float in3, float *out1, float *out2,
                            float *out3)
{
  PyObject *func, *arglist, *result, *tuple;

  func = (PyObject *)data;
  tuple = python_from_int_array(in1, in2);
  arglist = Py_BuildValue("(OOf)", model->scriptobj, tuple, in3);
  Py_XDECREF(tuple);
  if (!arglist) {
    return 1;
  }
  result = PyEval_CallObject(func, arglist);
  Py_DECREF(arglist);
  if (result && PySequence_Check(result) && PySequence_Size(result) == 3) {
    PyObject *o1, *o2, *o3;
    int ret;
    o1 = PySequence_GetItem(result, 0);
    o2 = PySequence_GetItem(result, 1);
    o3 = PySequence_GetItem(result, 2);
    ret = (!python_to_float_array(o1, in2, NULL, out1, "dervx")
           || !python_to_float_array(o2, in2, NULL, out2, "dervy")
           || !python_to_float_array(o3, in2, NULL, out3, "dervz"));
    Py_DECREF(o1);
    Py_DECREF(o2);
    Py_DECREF(o3);
    Py_DECREF(result);
    return ret;
  } else {
    if (result) {
      Py_DECREF(result);
      PyErr_SetString(PyExc_TypeError,
                      "Callback function should return (dervx, dervy, dervz)");
    }
    return 1;
  }
}

/* Evaluate a user-defined energy term */
static int python_cb_teval(void *data, const struct mod_model *model,
                           gboolean deriv, int atind[], int n_atind,
                           float *e_term, float dvx[], float dvy[],
                           float dvz[], const struct mod_energy_data *enedata,
                           const struct mod_libraries *libs)
{
  PyObject *func, *arglist, *result, *tuple, *derobj;

  func = (PyObject *)data;
  tuple = python_from_int_array(atind, n_atind);
  derobj = PyBool_FromLong(deriv);
  arglist = Py_BuildValue("(OOO)", model->scriptobj, derobj, tuple);
  Py_XDECREF(derobj);
  Py_XDECREF(tuple);
  if (!arglist) {
    return 1;
  }
  result = PyEval_CallObject(func, arglist);
  Py_DECREF(arglist);
  if (result && deriv && PySequence_Check(result)
      && PySequence_Size(result) == 4) {
    PyObject *o1, *o2, *o3, *o4;
    int ret;
    o1 = PySequence_GetItem(result, 0);
    o2 = PySequence_GetItem(result, 1);
    o3 = PySequence_GetItem(result, 2);
    o4 = PySequence_GetItem(result, 3);
    if (!PyNumber_Check(o1)) {
      PyErr_SetString(PyExc_ValueError, "Expected a number");
      ret = 1;
    } else {
      *e_term = PyFloat_AsDouble(o1);
      ret = !(python_to_float_array(o2, n_atind, NULL, dvx, "dvx")
              && python_to_float_array(o3, n_atind, NULL, dvy, "dvy")
              && python_to_float_array(o4, n_atind, NULL, dvz, "dvz"));
    }
    Py_DECREF(o1);
    Py_DECREF(o2);
    Py_DECREF(o3);
    Py_DECREF(o4);
    Py_DECREF(result);
    return ret;
  } else if (result && !deriv && PyNumber_Check(result)) {
    *e_term = PyFloat_AsDouble(result);
    Py_DECREF(result);
    return 0;
  } else {
    if (result) {
      Py_DECREF(result);
      PyErr_SetString(PyExc_TypeError,
                      deriv ? "Callback function should return "
                              "(val, dvx, dvy, dvz)"
                            : "Callback function should return a number");
    }
    return 1;
  }
}

/* Return the user-defined local similarity. */
static int python_cb_simloc(void *data, const struct mod_alignment *alignment, 
                            int ires, int is1, int is2, float *val)
{
  PyObject *func, *arglist, *result, *int1, *int2, *int3;

  func = (PyObject *)data;
  int1 = PyInt_FromLong(ires);
  int2 = PyInt_FromLong(is1);
  int3 = PyInt_FromLong(is2);
  arglist = Py_BuildValue("(OOOO)", alignment->scriptobj, int1, int2, int3);
  Py_XDECREF(int1);
  Py_XDECREF(int2);
  Py_XDECREF(int3);
  if (!arglist) {
    return 1;
  }
  result = PyEval_CallObject(func, arglist);
  Py_DECREF(arglist);
  if (result && PyFloat_Check(result)) {
    *val = PyFloat_AsDouble(result);
    Py_DECREF(result);
    return 0;
  } else if (result) {
    Py_DECREF(result);
    PyErr_SetString(PyExc_TypeError,
                    "Callback function should return a floating-point number");
  }
  return 1;
}

/* Evaluate a user-defined functional form. */
static int python_cb_reval(void *data, float feat[], int iftyp[],
                           int modal[], int n_feat, float pcsr[],
                           int n_pcsr, gboolean deriv, float fderv[],
                           float *val)
{
  PyObject *func, *arglist, *result, *tuple1, *tuple2, *tuple3, *tuple4,
           *derobj;

  func = (PyObject *)data;
  tuple1 = python_from_float_array(feat, n_feat);
  tuple2 = python_from_int_array(iftyp, n_feat);
  tuple3 = python_from_int_array(modal, n_feat);
  tuple4 = python_from_float_array(pcsr, n_pcsr);
  derobj = PyBool_FromLong(deriv);
  arglist = Py_BuildValue("(OOOOO)", tuple1, tuple2, tuple3, tuple4, derobj);
  Py_XDECREF(tuple1);
  Py_XDECREF(tuple2);
  Py_XDECREF(tuple3);
  Py_XDECREF(tuple4);
  Py_XDECREF(derobj);
  if (!arglist) {
    return 1;
  }
  result = PyEval_CallObject(func, arglist);
  Py_DECREF(arglist);
  if (result && deriv && PySequence_Check(result)
      && PySequence_Size(result) == 2) {
    PyObject *o1, *o2;
    int ret;
    o1 = PySequence_GetItem(result, 0);
    o2 = PySequence_GetItem(result, 1);
    if (!PyNumber_Check(o1)) {
      PyErr_SetString(PyExc_ValueError, "Expected a number");
      ret = 1;
    } else {
      *val = PyFloat_AsDouble(o1);
      ret = !python_to_float_array(o2, n_feat, NULL, fderv, "fderv");
    }
    Py_DECREF(o1);
    Py_DECREF(o2);
    Py_DECREF(result);
    return ret;
  } else if (result && !deriv && PyNumber_Check(result)) {
    *val = PyFloat_AsDouble(result);
    Py_DECREF(result);
    return 0;
  } else {
    if (result) {
      Py_DECREF(result);
      PyErr_SetString(PyExc_TypeError,
                      deriv ? "Callback function should return (val, fderv)"
                            : "Callback function should return a number");
    }
    return 1;
  }
}

/* Evaluate the violation of a user-defined restraint form */
static int python_cb_violation(void *data, float feat[], int iftyp[],
                               int modal[], int n_feat, float pcsr[],
                               int n_pcsr, float *val)
{
  PyObject *func, *arglist, *result, *tuple1, *tuple2, *tuple3, *tuple4;

  func = (PyObject *)data;
  tuple1 = python_from_float_array(feat, n_feat);
  tuple2 = python_from_int_array(iftyp, n_feat);
  tuple3 = python_from_int_array(modal, n_feat);
  tuple4 = python_from_float_array(pcsr, n_pcsr);
  arglist = Py_BuildValue("(OOOO)", tuple1, tuple2, tuple3, tuple4);
  Py_XDECREF(tuple1);
  Py_XDECREF(tuple2);
  Py_XDECREF(tuple3);
  Py_XDECREF(tuple4);
  if (!arglist) {
    return 1;
  }
  result = PyEval_CallObject(func, arglist);
  Py_DECREF(arglist);
  if (result && PyNumber_Check(result)) {
    *val = PyFloat_AsDouble(result);
    Py_DECREF(result);
    return 0;
  } else {
    if (result) {
      Py_DECREF(result);
      PyErr_SetString(PyExc_TypeError,
                      "Callback function should return a number");
    }
    return 1;
  }
}

/* Evaluate the mean of a user-defined restraint form */
static int python_cb_mean(void *data, float feat[], int iftyp[],
                          int modal[], int n_feat, float pcsr[],
                          int n_pcsr, float *val)
{
  PyObject *func, *arglist, *result, *tuple1, *tuple2, *tuple3, *tuple4;

  func = (PyObject *)data;
  tuple1 = python_from_float_array(feat, n_feat);
  tuple2 = python_from_int_array(iftyp, n_feat);
  tuple3 = python_from_int_array(modal, n_feat);
  tuple4 = python_from_float_array(pcsr, n_pcsr);
  arglist = Py_BuildValue("(OOOO)", tuple1, tuple2, tuple3, tuple4);
  Py_XDECREF(tuple1);
  Py_XDECREF(tuple2);
  Py_XDECREF(tuple3);
  Py_XDECREF(tuple4);
  if (!arglist) {
    return 1;
  }
  result = PyEval_CallObject(func, arglist);
  Py_DECREF(arglist);
  if (result) {
    int ret = !python_to_float_array(result, n_feat, NULL, val, "val");
    Py_DECREF(result);
    return ret;
  } else {
    return 1;
  }
}

/* Carry out an optimizer action */
static int python_cb_action(void *data, const void *optobj)
{
  PyObject *func, *arglist, *result;

  func = (PyObject *)data;
  if (!(arglist = Py_BuildValue("(O)", optobj))) {
    return 1;
  }
  result = PyEval_CallObject(func, arglist);
  Py_DECREF(arglist);
  if (result && result != Py_None) {
    PyErr_SetString(PyExc_ValueError,
                    "Optimizer actions should return nothing");
    Py_DECREF(result);
    return 1;
  } else if (!result) {
    return 1;
  } else {
    return 0;
  }
}

/* Log an output message */
static int python_cb_log(void *data, const char *text)
{
  PyObject *func, *arglist, *result;

  func = (PyObject *)data;
  arglist = Py_BuildValue("(s)", text);
  if (!arglist) {
    return 1;
  }
  result = PyEval_CallObject(func, arglist);
  Py_DECREF(arglist);
  if (result && result != Py_None) {
    PyErr_SetString(PyExc_ValueError,
                    "Log functions should return nothing");
    Py_DECREF(result);
    return 1;
  } else if (!result) {
    return 1;
  } else {
    return 0;
  }
}

/* Checks that the provided Python object is a function, and sets
   parameters necessary to call it from a C callback function. The
   Python reference count is increased to prevent the function from
   being garbage collected. */
static int set_python_callback(void **funcpt, void **data, PyObject *pyfunc,
                               void *cfunc)
{
  *funcpt = cfunc;
  *data = pyfunc;
  if (!PyCallable_Check(pyfunc)) {
    PyErr_SetString(PyExc_TypeError, "Need a callable object!");
    return 1;
  }
  Py_INCREF(pyfunc);
  return 0;
}

/* Helper function to release the reference to a Python function; used
   by set_python_callback_free() */
static void python_callback_decref(void *data)
{
  PyObject *obj = data;
  Py_DECREF(obj);
}

/* Set up a Python callback, similar to set_python_callback(). Additionally,
   a 'free' function is registered to release the reference to the Python
   function when it is no longer required. */
static int set_python_callback_free(void **funcpt, cb_free *freept,
                                    void **data, PyObject *pyfunc, void *cfunc)
{
  if (set_python_callback(funcpt, data, pyfunc, cfunc) == 0) {
    *freept = python_callback_decref;
    return 0;
  } else {
    return 1;
  }
}
%}
