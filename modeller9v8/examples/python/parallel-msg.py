from modeller import *
from modeller.parallel import *

# Create an empty parallel job, and then add a single slave process running
# on the local machine
j = job()
j.append(local_slave())

# Start all slave processes (note: this will only work if 'modxxx' - where
# xxx is the Modeller version - is in the PATH; if not, use modeller_path
# to specify an alternate location)
j.start()

# Have each slave read in a PDB file (provided by us, the master) and
# return the PDB resolution back to us
for slave in j:
    slave.run_cmd('''
env = environ()
env.io.atom_files_directory = ["../atom_files"]
log.verbose()
code = master.get_data()
mdl = model(env, file=code)
master.send_data(mdl.resolution)
''')
    slave.send_data('1fdn')
    data = slave.get_data()
    print slave, "returned model resolution: ", data
