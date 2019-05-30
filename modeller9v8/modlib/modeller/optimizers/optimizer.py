import _modeller

class optimizer(object):
    atmsel = None

    def optimize(self, atmsel):
        raise NotImplementedError

    def get_selection(self):
        return self.atmsel

    def _update_params(self, params, ok_keys, vars):
        for key in vars.keys():
            if key in ok_keys:
                params[key] = vars[key]
            else:
                raise ValueError("Unrecognized parameter: %s" % key)
