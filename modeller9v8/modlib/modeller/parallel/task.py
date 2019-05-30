import data_types
import os, tempfile, shutil
from communicator import TaskSetupError

class task(object):
    _results = None
    run_in_tempdir = False

    def __init__(self, *args, **vars):
        self.input_files = []
        self.output_files = []
        self._args = args
        self._vars = vars

    def __str__(self):
        return "<task>"

    def _setup(self):
        """Do pre-run setup, e.g. making a temporary run directory"""
        try:
            if self.run_in_tempdir:
                self._cwd = os.getcwd()
                self._tmpdir = tempfile.mkdtemp()
                os.chdir(self._tmpdir)
        except Exception, detail:
            # Wrap errors that occur in the setup phase
            raise TaskSetupError(detail)

    def _do_run(self, master):
        """Actually run the task"""
        try:
            ret = self.run(*self._args, **self._vars)
            for transfer in self.output_files:
                master.send_data(data_types.TransferFile(transfer))
        finally:
            if self.run_in_tempdir:
                os.chdir(self._cwd)
                shutil.rmtree(self._tmpdir, ignore_errors=True)
        return ret

    def run(self):
        raise NotImplementedError
