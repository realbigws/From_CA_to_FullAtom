from communicator import Communicator
import socket

class SlaveCommunicator(Communicator):
    connect_timeout = 600

    def __init__(self, lock=None, reconnect=None):
        Communicator.__init__(self, lock)
        self.reconnect = reconnect

    def connect_back(self, host, port, identifier):
        """Establish a TCP connection with the master"""
        s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        s.settimeout(self.connect_timeout)
        s.connect((host, port))
        s.sendall(identifier)
        s.settimeout(None)
        self.accept_connection(s)

    def _send(self, data):
        """Try to reopen the connection to the master if we got a broken pipe"""
        try:
            self.socket.sendall(data)
        except socket.error:
            if self.reconnect:
                self.connect_back(*self.reconnect)
                self.socket.sendall(data)
            else:
                raise
