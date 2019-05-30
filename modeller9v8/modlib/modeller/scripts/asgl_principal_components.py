def asgl_principal_components(env, family):
    """Produce an ASGL plot for the principal components analysis;
       requires the .dat file produced by principal_components()"""

    table_file = family + '.dat'
    asgl_file = family + '.top'

    fh = open(asgl_file, 'w')

    fh.write('READ_TABLE FILE = \'' + table_file + '\'\n')
    fh.write('SET XY_COLUMNS = 3 4\n')
    fh.write('WORLD POSITION = 12 1\n')
    fh.write('AXES2D\n')
    fh.write('PLOT2D PLOT2D_LINE_TYPE = 0, PLOT2D_SYMBOL_TYPE = 4, ;\n')
    fh.write('       LABEL_COLUMN = 2, LABEL_LOCATION = 2, LABEL_FONT = 5\n')
    fh.write('RESET_CAPTIONS\n')
    fh.write('CAPTION CAPTION_POSITION = 1, ;\n')
    fh.write('        CAPTION_TEXT = \'FAMILY ' + family + '\'\n')
    fh.write('CAPTION CAPTION_POSITION = 1, ;\n')
    fh.write('        CAPTION_TEXT = \'PRINCIPAL COMPONENTS CLUSTERING\'\n')
    fh.write('CAPTION CAPTION_POSITION = 2, CAPTION_TEXT = \'p\'\n')
    fh.write('CAPTION CAPTION_POSITION = 3, CAPTION_TEXT = \'q\'\n')
    fh.close()

    env.system('asgl ' + family)
