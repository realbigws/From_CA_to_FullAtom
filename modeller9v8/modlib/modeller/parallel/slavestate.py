class SlaveState:
    def __init__(self, desc):
        self.desc = desc

init = SlaveState("Initial state")
pending = SlaveState("Pending startup")
connected = SlaveState("Connected, with no active task")
running_task = SlaveState("Connected, and running a task")
dead = SlaveState("Disconnected due to failure")
