import os
import job, sge_qsub_array

class sge_qsub_job(job.job):
    """A parallel job which automatically starts Sun Grid Engine processes
       using 'qsub'."""

    def __init__(self, options, maxslave, seq=(), modeller_path=None,
                 host=None):
        job.job.__init__(self, seq, modeller_path, host)
        self.options = options
        self.maxslave = maxslave
        self.pending_arrays = []

    def expand_for_tasks(self):
        numslaves = len(self)
        numdesired = len(self.tasks)
        if self.maxslave is not None and self.maxslave < numdesired:
            numdesired = self.maxslave
        if numslaves < numdesired:
            array = sge_qsub_array.sge_qsub_array(self.options,
                                                  numdesired - numslaves)
            self.extend(array)
            self.pending_arrays.append(array)

    def start_processes(self, port):
        """Start all new slaves"""
        job.job.start_processes(self, port)
        jobname = self.get_name()
        for array in self.pending_arrays:
            array.start(jobname)
        self.pending_arrays = []
