import _modeller
from modeller.error import ModellerError

class Logger(object):
    """Allow messages to be printed to the Modeller logfile"""

    def __init__(self):
        self.header_written = False
        # Use our own write function for all log output
        _modeller.mod_log_function_set(self._write_with_header)

    def level(self, output=1, notes=0, warnings=0, errors=1, memory=0):
        """Set the logging level - i.e. which messages to output"""
        for (indx, val) in enumerate((output, notes, warnings, errors, memory)):
            _modeller.mod_log_set(indx, val)

    def none(self):
        """Don't print any log messages"""
        self.level(0, 0, 0, 0, 0)

    def minimal(self):
        """Only print important logging messages"""
        self.level(output=1, notes=0, warnings=0, errors=1, memory=0)

    def verbose(self):
        """Print verbose logging messages"""
        self.level(output=1, notes=1, warnings=1, errors=1, memory=1)

    def very_verbose(self):
        """Print verbose logging messages, with dynamic memory breakdown"""
        self.level(output=1, notes=1, warnings=1, errors=1, memory=2)

    def message(self, routine, message):
        """Print an informative message to the log"""
        if self.output:
            self.__wrterr(routine, '_', message)

    def warning(self, routine, message):
        """Print a warning message to the log"""
        if self.warnings:
            self.__wrterr(routine, 'W', message)

    def error(self, routine, message):
        """Print an error message to the log"""
        if self.errors:
            self.__wrterr(routine, 'E', message)
        raise ModellerError(message)

    def write_header_once(self):
        """Write the Modeller header to the log, but only once"""
        if not self.header_written:
            if self.output:
                self.write(_modeller.mod_header_get())
            self.header_written = True
        # Use standard method for further logging
        _modeller.mod_log_function_set(self.write)

    def _write_with_header(self, text):
        """Write to the logfile, with header prepended"""
        self.write_header_once()
        self.write(text)

    def write(self, text):
        """Write directly to the logfile"""
        print text

    def __wrterr(self, routine, code, message):
        if len(routine) > 10:
            routine = routine[:10]
        elif len(routine) < 10:
            routine = routine + '_' * (10 - len(routine))
        _modeller.mod_log_write("%s_%s> %s" % (routine, code, message))

    def __set_output(self, val):
        _modeller.mod_log_set(0, val)
    def __get_output(self):
        return _modeller.mod_log_get(0)
    def __set_notes(self, val):
        _modeller.mod_log_set(1, val)
    def __get_notes(self):
        return _modeller.mod_log_get(1)
    def __set_warnings(self, val):
        _modeller.mod_log_set(2, val)
    def __get_warnings(self):
        return _modeller.mod_log_get(2)
    def __set_errors(self, val):
        _modeller.mod_log_set(3, val)
    def __get_errors(self):
        return _modeller.mod_log_get(3)
    def __set_memory(self, val):
        _modeller.mod_log_set(4, val)
    def __get_memory(self):
        return _modeller.mod_log_get(4)
    output = property(__get_output, __set_output)
    notes = property(__get_notes, __set_notes)
    warnings = property(__get_warnings, __set_warnings)
    errors = property(__get_errors, __set_errors)
    memory = property(__get_memory, __set_memory)

# Global log object
log = Logger()
