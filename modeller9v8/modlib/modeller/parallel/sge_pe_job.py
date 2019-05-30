import os
import job, sge_pe_slave, local_slave

class sge_pe_job(job.job):
    """A parallel job containing processes on all Sun Grid Engine slave nodes"""

    def __init__(self, seq=(), modeller_path=None, host=None):
        job.job.__init__(self, seq, modeller_path, host)
        pe = os.environ['PE_HOSTFILE']
        fh = open(pe, "r")
        while True:
            line = fh.readline()
            if line == '':
                break
            (node, num, queue) = line.split(None, 2)
            for i in range(int(num)):
                self.append(sge_pe_slave.sge_pe_slave(node))
        # Replace first slave with a local_slave, as this is ourself, and SGE
        # won't let us start this process with qrsh (as we are already
        # occupying the slot)
        if len(self) > 0:
            self[0] = local_slave.local_slave()
