%module _cuser_optimizer

%{
#include "cuser_optimizer.h"
%}

%include "typemaps.i"
%apply int *OUTPUT { int * };

%include "cuser_optimizer.h"
