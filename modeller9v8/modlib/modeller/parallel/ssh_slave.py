import slave
from modeller.parallel.myspawn import myspawn

class ssh_slave(slave.slave):
    """A slave on a remote host accessed via ssh or rsh"""

    def __init__(self, nodename, ssh_command='ssh'):
        self._nodename = nodename
        self._ssh = ssh_command
        slave.slave.__init__(self)

    def start(self, path, id, output):
        slave.slave.start(self, path, id, output)
        myspawn("%s %s %s -slave %s" % (self._ssh, self._nodename, path, id),
                output)

    def __repr__(self):
        return "<ssh slave on %s>" % self._nodename
