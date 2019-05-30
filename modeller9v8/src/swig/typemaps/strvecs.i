/* Convert Python lists to NULL-terminated C string vectors */

#ifdef SWIGPYTHON
%typemap(in) (char **STRVEC) {
    $1 = python_to_c_str_array($input, "$1_name");
    if (!$1) SWIG_fail;
}
/* Override the generic char ** argout typemap */
%typemap(argout) (char **STRVEC) {
}
%typemap(arginit) (char **STRVEC) {
  $1 = NULL;
}
%typemap(freearg) (char **STRVEC) {
  g_strfreev($1);
}
#endif

/* Apply typemaps to actual functions */
%apply char **STRVEC { char **exts };
%apply char **STRVEC { char **dirs };
