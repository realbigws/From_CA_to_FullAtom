/* Convert Python/Perl strings to/from C strings */


/* All non-const strings returned from Modeller functions are strdup'd, so
   we must free them. const strings point to static storage, so we must not.
   These work in conjunction with %newobject definitions, auto-generated for
   all char * functions and placed in swig-arrays.i */
%typemap(newfree) const char * {
}
%typemap(newfree) char * {
  if ($1) free($1);
}

/* Return cstrdup'd strings as Python or Perl strings. Do not free const
   strings. */
%typemap(in, numinputs=0) char ** (char *temp = 0) {
  $1 = &temp;
}
#ifdef SWIGPYTHON
%typemap(argout,fragment="t_output_helper") const char ** {
  if (*$1) {
    PyObject *o = PyString_FromString(*$1);
    $result = t_output_helper($result, o);
  }
}
%typemap(argout,fragment="t_output_helper") char ** {
  if (*$1) {
    PyObject *o = PyString_FromString(*$1);
    free(*$1);
    $result = t_output_helper($result, o);
  }
}
#endif
#ifdef SWIGPERL
%typemap(argout) const char ** {
  if (*$1) {
    if (argvi >= items) {
      EXTEND(sp, 1);
    }
    $result = sv_newmortal();
    sv_setpv($result, *$1);
    argvi++;
  }
}
%typemap(argout) char ** {
  if (*$1) {
    if (argvi >= items) {
      EXTEND(sp, 1);
    }
    $result = sv_newmortal();
    sv_setpv($result, *$1);
    argvi++;
    free(*$1);
  }
}
#endif
