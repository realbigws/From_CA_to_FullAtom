/* Custom typemaps to handle arrays, callback functions, exceptions, etc. */


/* Standard SWIG typemaps */
%include "typemaps.i"

/* Helper functions used by typemaps */
#ifdef SWIGPYTHON
%include "helperfuncs/python-init.i"
%include "helperfuncs/python-arrays.i"
%include "helperfuncs/python-callbacks.i"
%include "typemaps/python-dstr.i"
#endif

#ifdef SWIGPERL
%include "helperfuncs/perl-arrays.i"
#endif

%include "helperfuncs/exceptions.i"

/* Pass scripting language objects to new_foo() routines as opaque void * */
%typemap(in) const void * scriptobj {
  $1 = $input;
}

%include "typemaps/bools.i"
%include "typemaps/strings.i"
%include "typemaps/buffers.i"
%include "typemaps/strvecs.i"
%include "typemaps/fixed-size-arrays.i"
%include "typemaps/variable-size-arrays.i"

%include "typemaps/swig-arrays.i"

/* Get simple return values */
%apply float *OUTPUT { float * };
%apply int *OUTPUT { int * };

%include "typemaps/exceptions.i"

#ifdef SWIGPYTHON
%include "typemaps/python-callbacks.i"
#endif
