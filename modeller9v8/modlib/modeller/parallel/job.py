import socket, random, os, select
import sys
import modeller
import _modeller
from communicator import NetworkError, TaskSetupError
import slavestate
import data_types

class job(list):
    connect_timeout = 7200
    heartbeat_timeout = 7200

    def __init__(self, seq=(), modeller_path=None, host=None):
        list.__init__(self, seq)
        self.slave_startup_commands = []
        self.tasks = []
        if modeller_path is not None:
            self.modeller_path = modeller_path
        else:
            self.modeller_path = self.get_default_modeller_path()
        if host:
            self.host = host
        else:
            # Get primary IP address of this machine
            self.host = socket.gethostbyname_ex(socket.gethostname())[-1][0]
        self.listensock = self.__listen_on_random_port(self.host,
                                                       self.connect_timeout)
        self.pending_slaves = {}
        self.connected_slaves = {}
        self.cwd = os.getcwd()

    def get_default_modeller_path(self):
        bindir = modeller.info.bindir
        modpy = os.path.join(bindir, 'modpy.sh')
        modslave = os.path.join(bindir, 'modslave.py')
        if os.path.exists(modslave) and not _modeller.mod_embedded_get():
            path = sys.executable + " " + modslave
            if os.path.exists(modpy):
                path = modpy + " " + path
            return path
        else:
            return "mod" + _modeller.mod_short_version_get()

    def get_name(self):
        job = _modeller.mod_jobname_get()
        if job == '(stdin)':
            job = 'stdout'
        else:
            job = os.path.basename(job)
        return job

    def start_processes(self, port):
        job = self.get_name()
        for (num, slave) in enumerate(self):
            if slave.get_state() == slavestate.init:
                id = self.__get_id(num)
                self.pending_slaves[id] = slave
                addr = "%s:%d" % (self.host, port)
                output = job + ".slave%d" % num
                slave.start(self.modeller_path, "%s:%s" % (addr, id), output)

    def expand_for_tasks(self):
        pass

    def accept_slave(self, sock, pending_slaves, connected_slaves):
        # Make sure the new socket is blocking (on some platforms this socket
        # inherits non-blocking status from the listening socket)
        sock.setblocking(True)
        id = sock.recv(1024)
        if id and id in pending_slaves:
            slave = pending_slaves.pop(id)
            connected_slaves[id] = slave
            print "Identified slave %s " % str(slave)
            slave.accept_connection(sock)
            slave.set_directory(self.cwd)
            if sys.path[0] != '':
                slave.set_python_search_path(sys.path[0])
            for cmd in self.slave_startup_commands:
                slave.run_cmd(cmd)
            slave.set_log_level(modeller.log)
            return slave
        elif id and id in connected_slaves:
            slave = connected_slaves[id]
            print "Reconnect from slave %s " % str(slave)
            slave.accept_connection(sock)
        else:
            print "Ignoring request from unknown slave"

    def start(self):
        """Start all non-running slaves"""
        (s, port) = self.listensock
        self.start_processes(port)
        while len(self.pending_slaves) > 0:
            (conn, addr) = s.accept()
            self.accept_slave(conn, self.pending_slaves, self.connected_slaves)
        print "All slaves connected OK"

    def __listen_on_random_port(self, host, timeout):
        """Open a listening socket on a random high-numbered port"""
        s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        s.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
        tries = 0
        while True:
            port = random.randint(10000, 60000)
            try:
                s.bind((host, port))
            # gaierror is a subclass of error, so catch it separately
            except socket.gaierror:
                raise
            except socket.error:
                tries += 1
                if tries > 10: raise
            else:
                break
        s.listen(15)
        s.settimeout(timeout)
        return (s, port)

    def queue_task(self, taskobj):
        self.tasks.append(taskobj)

    def run_all_tasks(self):
        """Run all tasks and return all the results, in the same order that they
           were submitted, when all the jobs have completed."""
        tasks = self.tasks[:]
        print "Running %d tasks on %d slaves" % (len(tasks), len(self))
        self.push_tasks_to_slaves()
        while True:
            try:
                for task in self._finish_all_tasks():
                    pass
            except IndexError:
                break
        if len(self.tasks) > 0:
            raise ValueError("Ran out of slaves to run tasks")
        return [task._results for task in tasks]

    def yield_tasks_unordered(self):
        """Run all tasks and return their results (as a generator), in
           whatever order they complete."""
        print "Running %d tasks on %d slaves" % (len(self.tasks), len(self))
        self.push_tasks_to_slaves()
        while True:
            try:
                for task in self._finish_all_tasks():
                    yield task._results
            except IndexError:
                break
        if len(self.tasks) > 0:
            raise ValueError("Ran out of slaves to run tasks")

    def _finish_all_tasks(self):
        """Waits for tasks to finish from any running slave; if any pending
           slaves try to connect, start them up. Return each finished task,
           as a generator."""
        (s, port) = self.listensock
        while True:
            events = self.get_next_events()
            if events is None:
                self.kill_all_running_slaves()
            else:
                for obj in events:
                    task = self._process_event(obj, s)
                    if task:
                        yield task

    def _process_event(self, obj, listensock):
        """Handle a single event returned from a slave"""
        if obj == listensock:
            # New slave just connected to the listening socket
            (conn, addr) = listensock.accept()
            slave = self.accept_slave(conn, self.pending_slaves,
                                      self.connected_slaves)
            if slave and len(self.tasks) > 0:
                slave.run_task(self.tasks.pop(0))
        elif obj.running_task():
            # A slave returned data
            try:
                task = obj.task_results()
                if task:
                    # The slave completed its task
                    print "%s on %s completed" % (str(task), str(obj))
                    if len(self.tasks) > 0:
                        obj.run_task(self.tasks.pop(0))
                    return task
                else:
                    # The slave sent back a heartbeat; check for any dead slaves
                    self.kill_timed_out_slaves()
            except (NetworkError, TaskSetupError), detail:
                self.kill_slaves((obj,), detail)
        else:
            print "Warning: slave %s reports data, but is not running a task" \
                  % str(obj)

    def kill_slaves(self, slaves, err=""):
        if err != "":
            err = "(%s) " % err
        for s in slaves:
            print "%s failed %s- removing from %s" % (s, err, self)
            task = s.kill()
            if task:
                self.tasks.append(task)
        self.push_tasks_to_slaves()

    def kill_all_running_slaves(self):
        running = [a for a in self if a.running_task()]
        self.kill_slaves(running)
        raise NetworkError("Did not hear from any running slave in %d seconds"
                           % self.heartbeat_timeout)

    def kill_timed_out_slaves(self):
        timedout = [a for a in self if a.running_task() and \
                                    a.contact_timeout(self.heartbeat_timeout)]
        if len(timedout) > 0:
            print "Did not hear from slaves %s in %d seconds" % \
                  (str(timedout), self.heartbeat_timeout)
            self.kill_slaves(timedout)

    def push_tasks_to_slaves(self):
        (s, port) = self.listensock
        self.start_processes(port)
        for slave in [a for a in self if a.ready_for_task()]:
            try:
                t = self.tasks.pop(0)
            except IndexError:
                break
            try:
                slave.run_task(t)
            # If a network error occurred, kill the slave and requeue the task
            except socket.error, detail:
                print "slave %s failed on run task with %s; removing from job" \
                      % (slave, detail)
                slave.kill()
                self.tasks.insert(0, t)
        self.expand_for_tasks()
        self.start_processes(port)

    def get_next_events(self):
        (s, port) = self.listensock
        running = [a for a in self if a.running_task()]
        slavemap = {}
        fileno = s.fileno()
        slavemap[fileno] = s
        if len(running) == 0:
            if len(self.tasks) == 0:
                raise IndexError("No more tasks")
            if len(self.pending_slaves) == 0:
                raise IndexError("No more slaves")
        # poll() for new events, or fall back to select() on platforms
        # which don't have poll().
        try:
            poll = select.poll()
        except AttributeError:
            poll = None
        if poll:
            return self.__next_events_poll(poll, running, slavemap, fileno)
        else:
            return self.__next_events_select(running, slavemap, fileno)

    def __next_events_poll(self, poll, running, slavemap, fileno):
        poll.register(fileno, select.POLLIN)
        for slave in running:
            fileno = slave.socket.fileno()
            slavemap[fileno] = slave
            poll.register(fileno, select.POLLIN)
        ready = poll.poll(self.heartbeat_timeout * 1000)
        if len(ready) == 0:
            return None
        else:
            return [slavemap[fd[0]] for fd in ready]

    def __next_events_select(self, running, slavemap, fileno):
        waitin = [ fileno ]
        for slave in running:
            fileno = slave.socket.fileno()
            slavemap[fileno] = slave
            waitin.append(fileno)
        (ready,rout,rerr) = select.select(waitin, [], [],
                                          self.heartbeat_timeout)
        if len(ready) == 0:
            return None
        else:
            return [slavemap[fd] for fd in ready]

    def __get_id(self, num):
        """Return a random identifier, used to make sure the right slaves
           connect back to us."""
        id = "%d:" % num
        for i in range(0, 8):
            id += chr(random.randint(0, 25) + ord('A'))
        return id

    def __repr__(self):
        return "<Parallel job [" + \
               ", ".join([str(node) for node in self]) + "]>"
