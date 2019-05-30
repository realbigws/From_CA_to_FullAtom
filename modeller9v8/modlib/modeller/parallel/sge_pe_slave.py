import slave
from modeller.parallel.myspawn import myspawn

class sge_pe_slave(slave.slave):
    """A slave within a Sun Grid Engine parallel environment"""

    def __init__(self, nodename):
        self._nodename = nodename
        slave.slave.__init__(self)

    def start(self, path, id, output):
        slave.slave.start(self, path, id, output)
        myspawn("qrsh -inherit -V %s %s -slave %s" % (self._nodename, path, id),
                output)

    def __repr__(self):
        return "<SGE PE slave on %s>" % self._nodename
