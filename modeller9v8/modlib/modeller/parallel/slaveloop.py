import socket
import communicator, data_types
from slave_communicator import SlaveCommunicator
import threading
import time

class HeartBeat(threading.Thread):
    """Periodically send a 'heartbeat' back to the master, so that it can
       distinguish between failed nodes and long calculations"""
    timeout = 300

    def __init__(self, master):
        threading.Thread.__init__(self)
        self.master = master
        self.event = threading.Event()

    def cancel(self):
        """Stop the heartbeat"""
        self.event.set()

    def run(self):
        while True:
            self.event.wait(self.timeout)
            if self.event.isSet():
                break
            else:
                self.master.send_data(data_types.heartbeat())

def handle_slave_io(master, slavedict):
    """Handle all messages from master"""
    while True:
        try:
            cmdstr = master.get_command()
        except communicator.NetworkError:
            # Connection broken - shutdown slave
            break
        try:
            exec(cmdstr, slavedict)
        except Exception, detail:
            # Propagate errors to master, and reraise them here (but don't
            # send back erorrs we got from the master!)
            if not isinstance(detail, communicator.RemoteError):
                try:
                    master.send_data(communicator.ErrorWrapper(detail))
                except socket.error, detail2:
                    print "Warning: ignored exception " + str(detail2) \
                          + " when trying to send error state " + str(detail) \
                          + " back to master"
            raise

def slaveloop(addr):
    print "Slave startup: connect to master at %s" % addr
    (host, port, identifier) = addr.split(":", 2)
    port = int(port)
    lock = threading.Lock()
    master = SlaveCommunicator(lock, reconnect=(host, port, identifier))
    master.connect_back(host, port, identifier)
    th = HeartBeat(master)
    th.start()
    slavedict = {'master':master}
    exec('from modeller import *', slavedict)
    try:
        handle_slave_io(master, slavedict)
    finally:
        th.cancel()
    print "Slave shutdown"
