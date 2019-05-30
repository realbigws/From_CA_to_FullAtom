import _modeller

def set_topvars(topvardict, vars, topini):
    for name in topvardict.keys():
        value = topvardict[name]
        if value is None:
            continue
        (typ, num) = topini[name]
        setfuncs = { 'S': [ str, [int, float] ],
                     'L': [ bool, [] ],
                     'R': [ float, [int] ],
                     'I': [ int, [] ] }
        (top_type, top_oktypes) = setfuncs[typ]
        if not isinstance(value, (list, tuple)):
            value = [value]
        check_top_dim(num, len(value))
        for val in value:
            val = check_top_type(name, val, top_type, top_oktypes)
        if typ == 'S':
            value = list(value)
            for indx in range(len(value)):
                value[indx] = str_macro(value[indx], vars)
        if num == 1:
            vars[name] = value[0]
        else:
            vars[name] = tuple(value)

def str_macro_sub(s, key, value):
    s = s.replace("${%s}" % key, value)
    return s.replace("$(%s)" % key, value)

def str_macro(value, vars):
    value = str(value)
    try:
        directory = vars['directory']
        root_name = vars['root_name']
        id1 = vars['id1']
        id2 = vars['id2']
        file_ext = vars['file_ext']
        file_id = vars['file_id']
    except KeyError:
        return value
    if value.find('${LIB}') >= 0 or value.find('$(LIB)') >= 0:
        value = str_macro_sub(value, 'LIB', _modeller.mod_libdir_get())
    if value.find('${DIR}') >= 0 or value.find('$(DIR)') >= 0:
        value = str_macro_sub(value, 'DIR', directory)
    if value.find('${JOB}') >= 0 or value.find('$(JOB)') >= 0:
        value = str_macro_sub(value, 'JOB', _modeller.mod_jobname_get())
    if value.find('${DEFAULT}') >= 0 or value.find('$(DEFAULT)') >= 0:
        value = str_macro_sub(value, 'DEFAULT',
                              "%s%s%04d%04d%s" % (root_name, file_id, id1,
                                                  id2, file_ext))
    return value

def check_top_type(name, value, top_type, top_oktypes):
    if not isinstance(value, top_type):
        if type(value) in top_oktypes:
            value = top_type(value)
        else:
            raise TypeError("Incorrect element type for variable %s"
                            % name.upper())
    return value

def check_top_dim(orig_dim, new_dim):
    if orig_dim != new_dim and orig_dim != 0:
        raise TypeError("Numbers of expected, actual arguments "
                        + "are different: %d  %d" % (orig_dim, new_dim))
