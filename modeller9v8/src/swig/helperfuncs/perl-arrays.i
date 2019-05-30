/* Convert a Perl list to a C array */
%define PERLTOARRAY(type, chkfunc, convfunc, errmsg)
%inline %{
static type *perl_to_ ## type ## _array(SV *input, int size, int *sizevar, char *varname) {
  int i;
  type *result;
  AV *tempav;
  if (!SvROK(input) || SvTYPE(SvRV(input)) != SVt_PVAV) {
    char errstr[160];
    sprintf(errstr, "Expected an array reference for %s", varname);
    croak(errstr);
  }
  tempav = (AV*)SvRV(input);
  if (sizevar) {
    *sizevar = size = av_len(tempav) + 1;
  } else if (av_len(tempav) + 1 != size) {
    char errstr[160];
    sprintf(errstr, "Expected a sequence of length %d for %s; got %d", size,
            varname, av_len(tempav) + 1);
    croak(errstr);
  }
  if (size == 0) {
    result = malloc(sizeof(type));
  } else {
    result = malloc(size * sizeof(type));
  }
  for (i = 0; i < size; i++) {
    SV **tv = av_fetch(tempav, i, 0);
    if (chkfunc(*tv)) {
      result[i] = (type)convfunc(*tv);
    } else {
      char errstr[160];
      sprintf(errstr, errmsg, varname, i);
      free(result);
      croak(errstr);
    }
  }
  return result;
}
%}
%enddef

PERLTOARRAY(int, SvIOK, SvIV, "%s[%d] must be an integer")
PERLTOARRAY(float, SvNIOK, SvNV, "%s[%d] must be a number")

%{
static char *perl_to_f_str_array(SV *input, int size, int *sizevar,
                                 int *strlenvar, char *varname) {
  int i, memsize;
  char *result;
  AV *tempav;
  *strlenvar = 1;
  if (!SvROK(input) || SvTYPE(SvRV(input)) != SVt_PVAV) {
    char errstr[160];
    sprintf(errstr, "Expected an array reference for %s", varname);
    croak(errstr);
  }
  tempav = (AV*)SvRV(input);
  if (sizevar) {
    *sizevar = size = av_len(tempav) + 1;
  } else if (av_len(tempav) + 1 != size) {
    char errstr[160];
    sprintf(errstr, "Expected a sequence of length %d for %s; got %d", size,
            varname, av_len(tempav) + 1);
    croak(errstr);
    return NULL;
  }
  for (i = 0; i < size; i++) {
    SV **tv = av_fetch(tempav, i, 0);
    if (SvPOK(*tv)) {
      int len;
      char *pt = SvPV(*tv, PL_na);
      len = strlen(pt);
      if (len > *strlenvar) *strlenvar = len;
    } else {
      char errstr[160];
      sprintf(errstr, "%s[%d] must be a string", varname, i);
      croak(errstr);
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
    SV **tv = av_fetch(tempav, i, 0);
    char *pt = SvPV(*tv, PL_na);
    memcpy(result + i * (*strlenvar), pt, strlen(pt));
  }
  return result;
}
%}
