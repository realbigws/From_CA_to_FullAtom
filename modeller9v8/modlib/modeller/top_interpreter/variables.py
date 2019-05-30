"""Classes to get/set legacy TOP variables"""

import _modeller
from modeller.top_interpreter import util

__docformat__ = "epytext en"

class Variables:
    """Class to get/set legacy TOP variables"""

    _deftops = { 'include_file':'__mod', 'variables':[''], 'arguments':[''],
                 'string_arguments':[''], 'operation':'', 'result':'',
                 'routine':'', 'numb_of_sequences':None,
                 'numb_of_sequences2':None, 'n_schedule':None,
                 'schedule_step':0, 'topology_model':3, 'rand_seed':-8123,
                 'seq_id':-999.0, 'output_control':[1, 0, 1, 1, 0],
                 'sel1':[1], 'sel1out':[1], 'sel2':[1], 'sel3':[1],
                 'sel':[1], 'selout':[1] }
    _ourvars = {}
    _aln_nseq = [ 0, 0 ]
    _edat_members = { 'lennard_jones_switch':[0], 'coulomb_switch':[0],
                      'contact_shell':0.0, 'relative_dielectric':0.0,
                      'radii_factor':0.0, 'update_dynamic':0.0,
                      'nonbonded_sel_atoms':0, 'covalent_cys':False,
                      'excl_local':[False], 'nlogn_use':0, 'sphere_stdv':0.0,
                      'dynamic_pairs':False, 'dynamic_sphere':False,
                      'dynamic_coulomb':False, 'dynamic_lennard':False,
                      'dynamic_modeller':False, 'dynamic_access':False }
    _io_members = { 'hydrogen_io':False, 'hetatm_io':False, 'water_io':False,
                    'atom_files_directory':'' }

    def __init__(self, topcmds):
        self._topcmds = topcmds
        self._topvars = {}
        for var in self._deftops.keys():
            self._ourvars[var] = self._deftops[var]
        for var in self._edat_members.keys():
            self._ourvars[var] = self._edat_members[var]
        for var in self._io_members.keys():
            self._ourvars[var] = self._io_members[var]

    def __setitem__(self, key, value):
        if key in self.topini and key not in self._deftops \
                 and key not in self._edat_members \
                 and key not in self._io_members:
            (typ, num) = self.topini[key]
            if typ == 'L':
                value = self._parse_bool(value)
            if num != 1 and not isinstance(value, (list, tuple)):
                value = [ value ]
            util.set_topvars({key:value}, self._topvars, self.topini)
            if key == 'align_codes':
                self._to_alignment(1, value, _modeller.mod_alignment_codes_set)
            elif key == 'align_codes2':
                self._to_alignment(2, value, _modeller.mod_alignment_codes_set)
            elif key == 'atom_files':
                self._to_alignment(1, value, self.__atom_files_set)
            elif key == 'atom_files2':
                self._to_alignment(2, value, self.__atom_files_set)
        else:
            if key in self._ourvars and self._ourvars[key] == None:
                raise SyntaxError(key.upper() + " is a read-only variable")
            if key in self._ourvars and (isinstance(self._ourvars[key], bool) or
               (isinstance(self._ourvars[key], list) and
                isinstance(self._ourvars[key][0], bool))):
                value = self._parse_bool(value)
            if key in self._edat_members:
                edat = self._topcmds.get_edat()
                exec("edat." + key + " = value")
            elif key == 'atom_files_directory':
                io = self._topcmds.get_io()
                io.atom_files_directory = value.split(':')
            elif key in self._io_members:
                io = self._topcmds.get_io()
                if key[-3:] == '_io':
                    key = key[:-3]
                exec("io." + key + " = value")
            elif key == 'schedule_step':
                self._set_sched_step(value)
            elif key == 'topology_model':
                self._set_topmod(value)
            elif key == 'rand_seed':
                self._set_rand_seed(value)
            elif key == 'seq_id':
                self._set_seq_id(value)
            elif key == 'output_control':
                self._set_output_control(value)
            elif key == 'sel1out':
                self['sel1'] = value
            elif key == 'selout':
                num = self['pick_atoms_set']
                num = min(3, max(1, num))
                names = { 1:'sel1', 2:'sel2', 3:'sel3' }
                self[names[num]] = value
            else:
                if isinstance(value, str):
                    value = util.str_macro(value, self._topvars)
                self._ourvars[key] = value

    def __getitem__(self, key):
        if key in self._ourvars:
            if key in self._edat_members:
                edat = self._topcmds.get_edat()
                return eval("edat." + key)
            elif key == 'atom_files_directory':
                io = self._topcmds.get_io()
                return ':'.join(io.atom_files_directory)
            elif key in self._io_members:
                io = self._topcmds.get_io()
                if key[-3:] == '_io':
                    key = key[:-3]
                return eval("io." + key)
            elif key == 'numb_of_sequences':
                return self._get_aln_nseq(1)
            elif key == 'numb_of_sequences2':
                return self._get_aln_nseq(2)
            elif key == 'n_schedule':
                return self._get_n_schedule()
            elif key == 'schedule_step':
                return self._get_sched_step()
            elif key == 'topology_model':
                return self._get_topmod()
            elif key == 'rand_seed':
                return self._get_rand_seed()
            elif key == 'seq_id':
                return self._get_seq_id()
            elif key == 'output_control':
                return self._get_output_control()
            elif key == 'sel':
                num = self['pick_atoms_set']
                num = min(3, max(1, num))
                names = { 1:'sel1', 2:'sel2', 3:'sel3' }
                return self[names[num]]
            else:
                return self._ourvars[key]
        elif key in self.topini:
            value = self._topvars[key]
            return value
        else:
            raise IndexError("%s not a valid TOP variable" % key)

    def __contains__(self, item):
        return (item in self._ourvars or item in self.topini)

    def _get_aln_nseq(self, num):
        """Get the number of sequences in the given alignment"""
        aln = self._topcmds.get_aln(num)
        return _modeller.mod_alignment_nseq_get(aln)

    def _get_sched(self):
        """Get the schedule object for MODEL"""
        mdl = self._topcmds.get_mdl(1)
        return _modeller.mod_model_sched_get(mdl)

    def _get_topmod(self):
        """Get the current topology submodel (TOPOLOGY_MODEL)"""
        tpl = self._topcmds.get_tpl()
        return _modeller.mod_topology_submodel_get(tpl)

    def _set_topmod(self, val):
        """Set the current topology submodel (TOPOLOGY_MODEL)"""
        tpl = self._topcmds.get_tpl()
        _modeller.mod_topology_submodel_set(tpl, val)

    def _get_rand_seed(self):
        """Get the current random seed"""
        libs = self._topcmds.get_libs()
        return _modeller.mod_libraries_rand_seed_get(libs)

    def _set_rand_seed(self, val):
        """Set the current random seed"""
        libs = self._topcmds.get_libs()
        _modeller.mod_libraries_rand_seed_set(libs, val)

    def _get_seq_id(self):
        """Get the model's seqid to best template (SEQ_ID)"""
        mdl = self._topcmds.get_mdl(1)
        return _modeller.mod_model_seq_id_get(mdl)

    def _set_seq_id(self, val):
        """Set the model's seqid to best template (SEQ_ID)"""
        mdl = self._topcmds.get_mdl(1)
        _modeller.mod_model_seq_id_set(mdl, val)

    def _get_output_control(self):
        """Get the current set of log levels (OUTPUT_CONTROL)"""
        vals = []
        for indx in range(5):
            vals.append(_modeller.mod_log_get(indx))
        return vals

    def _set_output_control(self, vals):
        """Set the current set of log levels (OUTPUT_CONTROL)"""
        for (indx, val) in enumerate(vals):
            _modeller.mod_log_set(indx, val)

    def _get_n_schedule(self):
        """Get the number of schedule steps (N_SCHEDULE)"""
        return _modeller.mod_schedule_nsch_get(self._get_sched())

    def _get_sched_step(self):
        """Get the current schedule step (SCHEDULE_STEP)"""
        return _modeller.mod_schedule_step_get(self._get_sched())

    def _set_sched_step(self, val):
        """Set the current schedule step (SCHEDULE_STEP)"""
        return _modeller.mod_schedule_step_set(self._get_sched(), val)

    def __atom_files_set(self, aln, num, val):
        """Helper function to set atom_files in an alignment"""
        alnseq = _modeller.mod_alignment_alnsequence_get(aln, num)
        _modeller.mod_alnsequence_atom_files_set(alnseq, val)

    def _to_alignment(self, num, value, setfunc):
        """Update the alignment's internal data as a result of setting
           ATOM_FILES or ALIGN_CODES"""
        aln = self._topcmds.get_aln(num)
        len_aln = _modeller.mod_alignment_nseq_get(aln)
        self._aln_nseq[num - 1] = len_aln
        numal = min(len(value), len_aln)
        for num in range(numal):
            setfunc(aln, num, value[num])

    def _parse_bool(self, varval):
        """Convert any TOP bools ('on' and 'off') to Python (True/False)"""
        notlist = False
        if not isinstance(varval, list):
            varval = [ varval ]
            notlist = True
        newvarval = []
        for var in varval:
            if isinstance(var, str):
                var = var.lower()
                if var == 'on':
                    newvarval.append(True)
                elif var == 'off':
                    newvarval.append(False)
                else:
                    newvarval.append(var)
            else:
                newvarval.append(var)
        if notlist:
            return newvarval[0]
        else:
            return newvarval
