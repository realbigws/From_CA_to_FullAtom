"""Classes to run legacy TOP scripts"""

import _modeller
from _modeller import mod_log_write
from modeller.top_interpreter.commands import Commands
from modeller.top_interpreter.variables import Variables

__docformat__ = "epytext en"

class TOPfile:
    """Class to run a legacy TOP script"""

    def __init__(self, topfile):
        """Prepare to run a TOP script from C{topfile}"""
        self.topcmds = Commands()
        self.variables = Variables(self.topcmds)
        self.topcmds.vars = self.variables
        # Read top.ini
        self._read_top_ini(self.variables)
        # Read libraries
        libs = self.topcmds.get_libs()
        _modeller.mod_libraries_read_libs(libs,
                                          self.variables['restyp_lib_file'])
        # Read default atom classes file
        _modeller.mod_atom_classes_read(self.topcmds.get_gprsr(),
                                        self.variables['atom_classes_file'])
        if isinstance(topfile, file):
            self.lines = self._read_top(topfile)
        else:
            # Set job name
            if topfile.endswith('.top'):
                _modeller.mod_jobname_set(topfile[:-4])
            else:
                _modeller.mod_jobname_set(topfile)
            fh = open(topfile, 'r')
            self.lines = self._read_top(fh)
            fh.close()

    def run(self):
        """Run the entire TOP script"""
        self._runlines(self.lines)

    def _read_top_ini(self, vars):
        """Parse the top.ini file containing variable definitions"""
        vars.topini = {}
        fh = open(_modeller.mod_inifil_get(), 'r')
        while True:
            line = fh.readline()
            if line.startswith('--- KEYWORDS:'):
                self._read_topvars(fh, vars)
            elif line == '':
                break
        fh.close()

    def _read_topvars(self, fh, vars):
        """Read all TOP variables from top.ini"""
        while True:
            line = fh.readline()
            if line == '' or line.startswith('--- END OF FILE'):
                break
            line = self._remove_comment(line.strip())
            if len(line) > 0:
                try:
                    self._parse_top_ini_line(line.split(None, 4), vars)
                except TypeError, detail:
                    mod_log_write("read_top__E> " + str(detail))
                    mod_log_write("             top.ini line: " + line)
                    raise

    def _parse_top_ini_line(self, words, vars):
        """Read a single line from top.ini"""
        typ = words[1][0]
        varname = words[2].lower()
        dim = int(words[3])
        vars.topini[varname] = (typ, dim)
        if len(words) > 4:
            try:
                vars[varname] = self._parse_var(words[4])
            except SyntaxError:    # Ignore 'read-only variable' errors
                pass
        elif dim == 0:
            vars[varname] = []
        else:
            raise TypeError("No default value given for variable %s!"
                            % varname.upper())

    def _read_top(self, fh):
        """Read and return all of the TOP command lines from the file C{fh}"""
        lines = []
        lastline = ''
        while True:
            line = fh.readline()
            if line == '': break
            line = line.expandtabs(1)
            line = self._remove_comment(line.strip())
            if len(line) > 0:
                if lastline.endswith(';'):
                    lines[-1] = lastline[:-1] + line
                else:
                    lines.append(line)
                lastline = lines[-1]
        return self._process_includes(lines)

    def _process_includes(self, lines):
        """Expand out any INCLUDE commands in C{lines}"""
        newlines = []
        for line in lines:
            (cmd, vars) = self._break_line(line)
            if cmd == 'INCLUDE':
                self._set_top_vars(vars)
                try:
                    fh = self._open_include_file(self.variables['include_file'])
                except IOError, detail:
                    mod_log_write("runlines__E> " + str(detail))
                    mod_log_write("             TOP Command line: " + line)
                    raise
                newlines.extend(self._read_top(fh))
                fh.close()
            else:
                newlines.append(line)
        return newlines

    def _open_include_file(self, fname):
        """Open and return an include file, looking in Modeller's include
           file search path and possibly adding a .top extension"""
        bindir = _modeller.mod_bindir_get()
        dirs = [ '.', bindir ]
        exts = [ '', '.top' ]
        for dir in dirs:
            for ext in exts:
                inc = dir + '/' + fname + ext
                try:
                    return open(inc, "r")
                except IOError:
                    pass
        raise IOError("Could not find include file: %s" % fname)

    def _runlines(self, lines):
        """Run a collection of TOP command lines"""
        subrout = self._process_subroutines(lines)
        callstack = []
        linenum = 0
        indxca = 0
        while indxca < len(lines):
            linenum = linenum + 1
            line = lines[indxca]
            loglevel = self.variables['output_control'][1]
            (cmd, vars) = self._break_line(line)
            try:
                if loglevel > 0 and cmd != 'SUBROUTINE':
                    self._write_act(line, linenum, indxca + 1)
                if cmd == 'DO':
                    indxca = self._handle_do(vars, indxca, lines, True)
                else:
                    self._set_top_vars(vars)
                    _modeller.mod_top_pre()
                    indxca = self._run_top_cmd(cmd, indxca, lines, callstack,
                                               subrout)
                    _modeller.mod_top_post()
            except (IndexError, SyntaxError, TypeError), detail:
                mod_log_write("runlines__E> " + str(detail))
                mod_log_write("             TOP Command line: " + line)
                raise

    def _write_act(self, line, linenum, indxca):
        """Write the current TOP command line out to the log"""
        lenwrap = 57
        line = line.strip()
        splstr = line.split('=')
        for i in range(len(splstr) - 1):
            if not splstr[i].endswith(' '):
                splstr[i] = splstr[i] + ' '
        line = '='.join(splstr)
        thisline = line[:lenwrap]
        line = line[lenwrap:]
        if len(line) > 0:
            thisline = thisline + ';'
        else:
            thisline = thisline + ' '
        mod_log_write('TOP_________>%6d%5d %s' % (linenum, indxca, thisline))
        while len(line) > 0:
            thisline = line[:lenwrap]
            line = line[lenwrap:]
            if len(line) > 0:
                thisline = thisline + ';'
            else:
                thisline = thisline + ' '
            mod_log_write('                      ' + thisline)
        mod_log_write('')

    def _handle_do(self, var, indxca, lines, firsttime):
        """Handle the DO command, returning the line number of the next
           command to execute (either right after the DO or after the END_DO)"""
        ind = var.find('=')
        varname = var[:ind].strip().lower()
        varval = var[ind+1:].strip()
        varval = self._parse_var(varval)
        if not isinstance(varval, list):
            varval = [varval]

        if len(varval) == 2:
            varval.append(1)
        if len(varval) != 3:
            raise SyntaxError("Wrong number of arguments to DO command")
        else:
            for var in varval:
                if not isinstance(var, (int, float)):
                    raise SyntaxError("DO commands require INTEGER or " +
                                      "REAL arguments")
            if firsttime:
                self.variables[varname] = varval[0]
            else:
                self.variables[varname] = self.variables[varname] + varval[2]
        if self.variables[varname] > varval[1]:
            return self._skip_to('END_DO', indxca, lines, incs=['DO'],
                                 decs=['END_DO']) + 1
        else:
            return indxca + 1

    def _process_subroutines(self, lines):
        """Look for any TOP subroutine definitions in C{lines}"""
        subrout = {}
        for num, line in enumerate(lines):
            (cmd, vars) = self._break_line(line)
            if cmd == 'SUBROUTINE':
                self._set_top_vars(vars)
                rname = self.variables['routine']
                if rname in subrout:
                    mod_log_write("Warning: subroutine %s redefined" % (rname))
                subrout[rname] = num
        return subrout

    def _split_quoted(self, line, splch):
        """Break up C{line} on C{splch}, except when it's within a
           quoted string"""
        start = 0
        end = 0
        numquote = 0
        last = ''
        splist = []
        for ch in line:
            if ch == "'" and last != '\\':
                numquote += 1
            if ch == splch and numquote % 2 == 0:
                splist.append(line[start:end])
                start = end + 1
            last = ch
            end += 1
        if end > start:
            splist.append(line[start:end])
        return splist

    def _break_line(self, line):
        """Break up a TOP command line into a command and variable
           assignments"""
        ind = line.find(' ')
        if ind >= 0:
            cmd = line[:ind].upper()
            vars = line[ind+1:].strip()
            if len(vars) > 0:
                if cmd == 'DO':
                    vars = vars.replace(',', ' ')
                else:
                    vars = self._split_quoted(vars, ',')
        else:
            cmd = line.upper()
            vars = ''
        return (cmd, vars)

    def _remove_comment(self, line):
        """Strip any trailing comment from a TOP command line"""
        ind = line.find('#')
        if ind >= 0:
            return line[:ind]
        else:
            return line

    def _set_top_vars(self, vars):
        """Set TOP variables by parsing any variable=value pairs"""
        for var in vars:
            if len(var) == 0: continue
            ind = var.find('=')
            if ind < 0:
                ind = var.find(' ')
            if ind < 0:
                raise IndexError("No '=' found in variable assignment")
            varname = var[:ind].strip().lower()
            varval = var[ind+1:].strip()
            varval = self._parse_var(varval,
                                     varname!='result' and varname!='variables')
            if varname in self.variables:
                self.variables[varname] = varval
            else:
                raise IndexError("Variable name not recognized: %s"
                                 % varname.upper())

    def _parse_var(self, cmdlin, subvar=True):
        """Convert a TOP value into a Python value or list of values"""
        vars = []
        curvar = ''
        delim = ''
        last = ''
        for ch in cmdlin:
            if ch == "'" and last != '\\':
                if delim == "'":
                    vars.append(curvar)
                    curvar = ''
                    delim = ''
                else:
                    delim = "'"
                    if len(curvar) > 0:
                        vars.append(self._get_top_rep(curvar))
                        curvar = ''
            elif ch == ' ' and delim != "'":
                if len(curvar) > 0:
                    vars.append(self._get_top_rep(curvar))
                    curvar = ''
            else:
                if ch == "'" and last == '\\':
                    curvar = curvar[:-1] + ch
                else:
                    curvar = curvar + ch
            last = ch
        if delim == "'":
            raise SyntaxError("Unmatched quote")
        if len(curvar) > 0:
            vars.append(self._get_top_rep(curvar))
        if subvar:
            # Allow for a TOP oddity: if a string value matches the name of a
            # TOP variable, substitute it in (as if the string were unquoted).
            newvars = []
            for var in vars:
                if isinstance(var, str) and var.upper() == var:
                    varlow = var.lower()
                    try:
                        getvar = self.variables[var.lower()]
                        if isinstance(getvar, (list, tuple)):
                            newvars.extend(getvar)
                        else:
                            newvars.append(getvar)
                    except IndexError:
                        newvars.append(var)
                else:
                    newvars.append(var)
            vars = newvars
        if len(vars) == 1:
            vars = vars[0]
        return vars

    def _get_top_rep(self, curvar):
        """Convert a string to an int/float/string variable using TOP rules"""
        try:
            return int(curvar)
        except ValueError:
            try:
                retval = float(curvar)
                if int(retval) == retval:
                    return int(retval)
                else:
                    return retval
            except ValueError:
                return curvar

    def __atom_files_get(self, aln, num):
        """Helper function to get atom_files from an alignment"""
        alnseq = _modeller.mod_alignment_alnsequence_get(aln, num)
        return _modeller.mod_alnsequence_atom_files_get(alnseq)

    def _sync_alignment(self, num, topname, getfunc):
        """Make sure that the TOP variable (ALIGN_CODES or ATOM_FILES) matches
           the actual content of the alignment"""
        aln = self.topcmds.get_aln(num)
        len_aln = _modeller.mod_alignment_nseq_get(aln)
        value = []
        for i in range(len_aln):
            value.append(getfunc(aln, i))
        self.variables[topname] = value

    def _run_top_cmd(self, cmd, indxca, lines, callstack, subrout):
        """Run a single TOP command, and return the number of the next line"""
        excludes = [ 'SET', 'INCLUDE', 'DEFINE_STRING', 'DEFINE_INTEGER',
                     'DEFINE_REAL', 'DEFINE_LOGICAL', 'SUBROUTINE', 'CALL',
                     'END_SUBROUTINE', 'RETURN', 'OPERATE', 'STRING_OPERATE',
                     'IF', 'END_IF', 'STRING_IF', 'ELSE', 'DO', 'END_DO',
                     'EXIT', 'CYCLE', 'WRITE_TOP', 'RESET' ]
        if cmd not in excludes:
            lowcmd = cmd.lower()
            try:
                eval('self.topcmds.' + lowcmd + '()')
            except AttributeError:
                raise SyntaxError("Invalid TOP command: " + cmd)
            except (_modeller.ModellerError, ZeroDivisionError, IOError,
                    MemoryError, EOFError):
                self.variables['error_status'] = 1
                _modeller.mod_top_error(1, self.variables['stop_on_error'])
            if lowcmd in ('read_alignment', 'expand_alignment',
                          'sequence_search', 'delete_alignment'):
                self._sync_alignment(1, 'align_codes',
                                     _modeller.mod_alignment_codes_get)
                self._sync_alignment(1, 'atom_files', self.__atom_files_get)
            elif lowcmd == 'read_alignment2':
                self._sync_alignment(2, 'align_codes2',
                                     _modeller.mod_alignment_codes_get)
                self._sync_alignment(2, 'atom_files2', self.__atom_files_get)
        elif cmd == 'WRITE_TOP':
            self._write_top(self.variables['file'], lines)
        elif cmd == 'RESET':
            self._reset()
        elif cmd == 'DEFINE_STRING':
            self._make_top_variables('')
        elif cmd == 'DEFINE_REAL':
            self._make_top_variables(0.0)
        elif cmd == 'DEFINE_INTEGER':
            self._make_top_variables(0)
        elif cmd == 'DEFINE_LOGICAL':
            self._make_top_variables(False)
        elif cmd == 'SUBROUTINE':
            return self._skip_to('END_SUBROUTINE', indxca, lines) + 1
        elif cmd == 'CYCLE':
            return self._skip_to('END_DO', indxca, lines, incs=['DO'],
                                 decs=['END_DO'])
        elif cmd == 'EXIT':
            return self._skip_to('END_DO', indxca, lines, incs=['DO'],
                                 decs=['END_DO']) + 1
        elif cmd == 'END_DO':
            indxca = self._skip_to('DO', indxca, lines, incs=['END_DO'],
                                   decs=['DO'], indinc=-1)
            line = lines[indxca]
            (cmd, vars) = self._break_line(line)
            return self._handle_do(vars, indxca, lines, False)
        elif cmd == 'IF' or cmd == 'STRING_IF':
            # Make sure there is a terminating END_IF
            self._skip_to('END_IF', indxca, lines, incs=['IF', 'STRING_IF'],
                          decs=['END_IF'])
            if cmd == 'IF':
                res = self._numeric_if(self.variables['operation'],
                                       self.variables['arguments'])
            else:
                res = self._string_if(self.variables['operation'],
                                      self.variables['string_arguments'])
            if not res:
                return self._skip_to(['ELSE', 'END_IF'], indxca, lines,
                                     incs=['IF', 'STRING_IF'],
                                     decs=['END_IF']) + 1
        elif cmd == 'ELSE':
            return self._skip_to('END_IF', indxca, lines,
                                 incs=['IF', 'STRING_IF'], decs=['END_IF']) + 1
        elif cmd == 'CALL':
            callstack.append(indxca)
            return self._jump_to_subroutine(self.variables['routine'],
                                            subrout) + 1
        elif cmd == 'END_SUBROUTINE' or cmd == 'RETURN':
            if len(callstack) > 0:
                return callstack.pop() + 1
            else:
                raise SyntaxError(cmd + " but not within a SUBROUTINE")
        elif cmd == 'STRING_OPERATE':
            res = self._string_operate(self.variables['operation'],
                                       self.variables['string_arguments'])
            name = self.variables['result'].lower()
            if name in self.variables:
                self.variables[name] = res
            else:
                raise IndexError("Bad variable name for RESULT: %s"
                                 % name.upper())
        elif cmd == 'OPERATE':
            res = self._operate(self.variables['operation'],
                                self.variables['arguments'])
            name = self.variables['result'].lower()
            if name in self.variables:
                self.variables[name] = res
            else:
                raise IndexError("Bad variable name for RESULT: %s"
                                 % name.upper())
        return indxca + 1

    def _write_top(self, filename, lines):
        """Write the entire TOP script out to a file"""
        fh = open(filename, "w")
        for line in lines:
            fh.write(line + '\n')
        fh.close()

    def _reset(self):
        """Reset the state of the TOP interpreter"""
        self.variables = Variables(self.topcmds)
        self.topcmds.vars = self.variables
        # Read top.ini
        self._read_top_ini(self.variables)

    def _string_if(self, operation, arguments):
        """Process the STRING_IF TOP command"""
        operation = operation.upper()
        if operation == 'EQ':
            return arguments[0] == arguments[1]
        elif operation == 'NE':
            return arguments[0] != arguments[1]
        elif operation == 'INDEX':
            return arguments[0].find(arguments[1]) >= 0
        else:
            raise SyntaxError("Invalid OPERATION for STRING_IF: %s"
                              % operation)

    def _numeric_if(self, operation, orig_arguments):
        """Process the (numeric) IF TOP command"""
        operation = operation.upper()
        arguments = []
        for arg in orig_arguments:
            try:
                arguments.append(float(arg))
            except ValueError:
                raise SyntaxError("Argument for numeric IF must be REAL " +
                                  "or INTEGER: " + str(arg))
        if operation == 'EQ':
            return abs(arguments[0] - arguments[1]) < 1e-10
        elif operation == 'GT':
            return arguments[0] > arguments[1]
        elif operation == 'LT':
            return arguments[0] < arguments[1]
        elif operation == 'GE':
            return arguments[0] >= arguments[1] - 1e-10
        elif operation == 'LE':
            return arguments[0] <= arguments[1] + 1e-10
        elif operation == 'NE':
            return abs(arguments[0] - arguments[1]) > 1e-10
        else:
            raise SyntaxError("Invalid OPERATION for IF: %s" % operation)

    def _string_operate(self, operation, arguments):
        """Process the STRING_OPERATE TOP command"""
        if operation.upper() == 'CONCATENATE':
            val = ''
            for arg in arguments: val += str(arg)
            # The old Fortran code did this, since strings can't end
            # in whitespace
            return val.replace('@', ' ')
        else:
            raise SyntaxError("Invalid OPERATION for STRING_OPERATE: "
                              + operation)

    def _operate(self, operation, arguments):
        """Process the OPERATE TOP command"""
        operation = operation.upper()
        if operation == 'SUM':
            return sum(arguments)
        elif operation == 'MULTIPLY':
            val = 0
            for arg in arguments: val = val * arg
            return val
        elif operation == 'DIVIDE':
            return arguments[0] / arguments[1]
        elif operation == 'POWER':
            return arguments[0] ** arguments[1]
        elif operation == 'MOD':
            return arguments[0] % arguments[1]
        else:
            raise SyntaxError("Invalid OPERATION for OPERATE: %s" % operation)

    def _skip_to(self, gotocmds, indxca, lines, incs=[], decs=[], indinc = 1):
        """Return the line number of the next command in C{gotocmds}"""
        if not isinstance(gotocmds, list):
            gotocmds = [gotocmds]
        inclevel = 0
        indxca = indxca + indinc
        while indxca < len(lines) and indxca >= 0:
            (cmd, vars) = self._break_line(lines[indxca])
            if cmd in gotocmds and inclevel == 0:
                return indxca
            elif cmd in incs:
                inclevel = inclevel + 1
            elif cmd in decs:
                inclevel = inclevel - 1
                if inclevel < 0:
                    raise SyntaxError("Mismatched " + str(incs)
                                      + " and " + str(decs))
            indxca = indxca + indinc
        raise SyntaxError("Cannot find matching " + str(gotocmds))

    def _jump_to_subroutine(self, rname, subrout):
        """Get the line number of the given subroutine"""
        if rname in subrout:
            return subrout[rname]
        else:
            raise SyntaxError("No such subroutine: " + rname)

    def _make_top_variables(self, initval):
        """Make new user-defined TOP variables"""
        vars = self.variables['variables']
        if not isinstance(vars, list):
            vars = [ vars ]
        for var in vars:
            var = var.lower()
            self.variables[var] = initval
