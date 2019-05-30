import slave
from modeller.parallel.myspawn import myspawn

class local_slave(slave.slave):
    """A slave running on the same machine"""

    def start(self, path, id, output):
        slave.slave.start(self, path, id, output)
        myspawn("%s -slave %s" % (path, id), output)

    def __repr__(self):
        return "<Slave on localhost>"
