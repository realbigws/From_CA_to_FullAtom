import xdrlib
import os

try:
    import cPickle as pickle
except ImportError:
    import pickle

class TransferFile(object):
    """Object used to transfer a file between Communicators"""
    def __init__(self, filename):
        self._filename = filename

class heartbeat:
    pass

class data_type:
    desc = "generic data"
    def recv(self, buffer):
        pass
    def send(self, obj):
        pass
    def __str__(self):
        return self.desc

class netstring(data_type):
    code = 'S'
    obj = str
    desc = "string"
    def recv(self, buffer):
        p = xdrlib.Unpacker(buffer[1:])
        obj = p.unpack_string()
        return (obj, buffer[1+p.get_position():])
    def send(self, obj):
        p = xdrlib.Packer()
        p.pack_string(obj)
        return self.code + p.get_buffer()

class netinteger(data_type):
    code = 'I'
    obj = int
    desc = "integer"
    def recv(self, buffer):
        p = xdrlib.Unpacker(buffer[1:])
        obj = p.unpack_int()
        return (obj, buffer[1+p.get_position():])
    def send(self, obj):
        p = xdrlib.Packer()
        p.pack_int(obj)
        return self.code + p.get_buffer()

class netfloat(data_type):
    code = 'F'
    obj = float
    desc = "floating point"
    def recv(self, buffer):
        p = xdrlib.Unpacker(buffer[1:])
        obj = p.unpack_float()
        return (obj, buffer[1+p.get_position():])
    def send(self, obj):
        try:
            p = xdrlib.Packer()
            p.pack_float(obj)
            return self.code + p.get_buffer()
        # Python < 2.5 can fail trying to send Inf or NaN floats, so fall back
        # to pickling in this case:
        except SystemError:
            return netpickle.send(obj)

class netpickle(data_type):
    code = 'P'
    obj = None
    desc = "Python pickled object"
    def recv(self, buffer):
        p = xdrlib.Unpacker(buffer[1:])
        obj = p.unpack_string()
        return (pickle.loads(obj), buffer[1+p.get_position():])
    def send(self, obj):
        p = xdrlib.Packer()
        try:
            p.pack_string(pickle.dumps(obj, -1))
        # Python < 2.5 can fail trying to send Inf or NaN floats in binary
        # mode, so fall back to the old protocol in this case:
        except SystemError:
            p.pack_string(pickle.dumps(obj, 0))
        return self.code + p.get_buffer()

class NetFile(data_type):
    code = 'f'
    obj = TransferFile
    desc = "Transferred file"
    def recv(self, buffer):
        filename, buffer = netstring.recv(buffer[1:])
        filelen, buffer = netinteger.recv(buffer)
        if len(buffer) < filelen:
            raise IndexError("File not completely read")
        filename = os.path.basename(filename)
        f = file(filename, 'wb')
        f.write(buffer[:filelen])
        return TransferFile(filename), buffer[filelen:]
    def send(self, obj):
        filename = obj._filename
        f = file(filename, 'rb')
        buf = f.read()
        return self.code + netstring.send(filename) \
               + netinteger.send(len(buf)) + buf

class netcommand(netstring):
    code = 'C'
    desc = "command"

typemap = {}
cmdmap = {}
for name in dir():
    var = eval(name)
    try:
        if issubclass(var, data_type) and var is not data_type:
            exec("%s = %s()" % (name, name))
            var = eval(name)
            cmdmap[var.code] = var
            typemap[var.obj] = var
    except TypeError:
        pass
