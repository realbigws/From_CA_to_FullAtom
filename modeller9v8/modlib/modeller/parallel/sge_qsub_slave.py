import os
import re
import weakref
import slave
from mypopen4 import MyPopen4

class sge_qsub_slave(slave.slave):
    """A slave started with Sun Grid Engine 'qsub'"""
    standard_options = '-j y -cwd -r n -V'

    def __init__(self, options, array=None):
        self._jobid = "(unknown)"
        self._options = options
        slave.slave.__init__(self)
        if array is not None:
            self._array = weakref.ref(array)

    def __del__(self):
        if self._jobid != "(unknown)":
            os.system("qdel %s &" % self._jobid)

    def _get_array(self):
        if hasattr(self, "_array"):
            return self._array()
        else:
            return None

    def start(self, path, id, output):
        slave.slave.start(self, path, id, output)
        array = self._get_array()
        if array is None:
            self._start_single(path, id, output)
        else:
            self._start_array(path, id, output, array)

    def _start_single(self, path, id, output):
        name = os.path.basename(output)
        qsub = "qsub -S /bin/sh -o '%s' -N '%s' %s %s" % \
              (output, name, self.standard_options, self._options)
        cmd = "%s -slave %s" % (path, id)
        print "%s | %s" % (cmd, qsub)
        a = MyPopen4(qsub)
        (input, output) = (a.stdin, a.stdout)
        print >> input, cmd
        input.close()
        outlines = output.readlines()
        output.close()
        for line in outlines:
            print line.rstrip('\r\n')
        a.require_clean_exit()
        self._set_jobid(outlines)

    def _start_array(self, path, id, output, array):
        array.start_slave(path, id, output, self.standard_options)

    def _set_jobid(self, outlines):
        """Try to figure out the job ID from the SGE qsub output"""
        if len(outlines) > 0:
            m = re.compile(r"\d+").search(outlines[0])
            if m:
                self._jobid = int(m.group())

    def __repr__(self):
        return "<SGE qsub slave, ID %s>" % self._jobid
