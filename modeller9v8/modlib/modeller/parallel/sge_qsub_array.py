import re
from modeller.parallel.sge_qsub_slave import sge_qsub_slave
from modeller.parallel.mypopen4 import MyPopen4

class sge_qsub_array(list):
    """An array of slaves started with Sun Grid Engine 'qsub'"""

    def __init__(self, options, numslave, seq=()):
        list.__init__(self, seq)
        self._options = options
        self._pending_slaves = []
        for i in range(numslave):
            self.append(sge_qsub_slave(options, array=self))

    def start_slave(self, path, id, output, standard_options):
        self._path = path
        self._standard_options = standard_options
        self._pending_slaves.append((id, output))

    def start(self, jobname):
        qsub = "qsub -S /bin/sh -N '%s' -o sge-errors %s %s -t 1-%d" % \
               (jobname, self._options, self._standard_options,
                len(self._pending_slaves))
        print qsub
        a = MyPopen4(qsub)
        (input, output) = (a.stdin, a.stdout)
        id = " ".join([repr(s[0]) for s in self._pending_slaves])
        out = " ".join([repr(s[1]) for s in self._pending_slaves])
        print >> input, "#!/bin/sh"
        print >> input, "id=( '' %s )" % id
        print >> input, "out=( '' %s )" % out
        print >> input, "myid=${id[$SGE_TASK_ID]}"
        print >> input, "myout=${out[$SGE_TASK_ID]}"
        print >> input, "%s -slave $myid > $myout 2>&1" % self._path
        input.close()
        outlines = output.readlines()
        output.close()
        for line in outlines:
            print line.rstrip('\r\n')
        a.require_clean_exit()
        self._set_jobid(outlines)

    def _set_jobid(self, outlines):
        """Try to figure out the job ID from the SGE qsub output"""
        if len(outlines) > 0:
            m = re.compile(r"\d+").search(outlines[0])
            if m:
                self._jobid = int(m.group())
                for (num, slave) in enumerate(self):
                    slave._jobid = "%d.%d" % (self._jobid, num+1)
