"""Classes to optionally build comparative models, and then refine the loops"""

from automodel import automodel
from modeller.util.modobject import modobject
from modeller import *
from modeller.optimizers import conjugate_gradients
import refine
import autosched

__docformat__ = "epytext en"

class loopmodel(automodel):
    """Optionally build comparative models, and then refine the loops"""
    loop = None
    inimodel = None
    loop_potential = None

    def __init__(self, env, sequence, alnfile=None, knowns=[], inimodel=None,
                 deviation=None, library_schedule=None, csrfile=None,
                 inifile=None, assess_methods=None, loop_assess_methods=None):
        automodel.__init__(self, env, alnfile, knowns, sequence, deviation,
                           library_schedule, csrfile, inifile, assess_methods)
        self.inimodel = inimodel
        self.loop = LoopData(env)
        self.loop.assess_methods = loop_assess_methods

    def make(self, exit_stage=0):
        """Build all models"""
        if self.loop.write_selection_only and self.final_malign3d:
            raise ValueError("Cannot turn on both loop.write_selection_only "
                             "and final_malign3d; please turn off one or both")
        if self.inimodel:
            self.env = self.loop.env
            self.build_seq(self.inimodel, 1)
        else:
            automodel.make(self, exit_stage)
        self.write_summary(self.loop.outputs, 'loop models')

    def get_loop_model_filename(self, sequence, id1, id2, file_ext):
        """Returns the loop model PDB name - usually of the form
           foo.BL000X000Y.pdb"""
        return modfile.default(file_id='.BL', file_ext=file_ext,
                               root_name=sequence, id1=id1, id2=id2)

    def multiple_models(self, atmsel):
        """Build all models, given the restraints and atom selection"""
        automodel.multiple_models(self, atmsel)
        envcopy = self.env
        self.env = self.loop.env
        try:
            for models in [a for a in self.outputs if a['failure'] is None]:
                num = models['num']
                filename = models['name']
                self.build_seq(filename, num)
        finally:
            self.env = envcopy

    def fit_models_on_template(self):
        """Superpose each of the generated models on the templates"""
        automodel.fit_models_on_template(self)
        for models in [a for a in self.outputs if a['failure'] is None]:
            num = models['num']
            filename = models['name']
            mdl = model(self.env, file=filename)
            aln = alignment(self.env)
            aln.append_model(mdl, align_codes=self.sequence)
            aln.append_model(mdl, align_codes=self.sequence)
            for loopmodels in [a for a in self.loop.outputs \
                               if a['failure'] is None and a['num'] == num]:
                filename = loopmodels['name']
                id1 = loopmodels['loopnum']
                mdl2 = model(self.env, file=filename)
                atmsel = selection(mdl).only_atom_types('CA')
                atmsel.superpose(mdl2, aln, fit=True)
                filename = self.get_loop_model_filename(self.sequence, id1,
                                                        num, '_fit.pdb')
                mdl2.write(file=filename)

    def read_potential(self):
        """Reads in the group_restraints statistical potential used for
           loop modeling. Redefine if you want to use a different potential"""
        return group_restraints(self.env, classes='$(LIB)/atmcls-melo.lib',
                                parameters='$(LIB)/melo1-dist.lib')

    def read_loop_alignment(self, filename):
        """Read in the alignment used to build loop models"""
        aln = alignment(self.env)
        for i in range(2):
            aln.append_model(self, align_codes=self.sequence,
                             atom_files=filename)
        return aln

    def create_loop_topology(self, aln, filename):
        """Build the initial topology of the loop"""
        self.clear_topology()
        self.generate_topology(aln[-1])

        # Save self.seq_id, since otherwise transfer_xyz will set it to 100%
        seq_id = self.seq_id
        self.transfer_xyz(aln)
        self.seq_id = seq_id

        # Fill in any missing atoms using internal coordinates
        self.build(build_method='INTERNAL_COORDINATES', initialize_xyz=False)

        self.res_num_from(model(self.env, file=filename), aln)
        self.special_patches(aln)

    def build_seq(self, filename, num):
        """Build all loop models for a given starting model"""
        self.csrfile = self.sequence + '.lrsr'
        self.read_top_par()
        if self.loop_potential is None:
            self.loop_potential = self.read_potential()
        oldgprsr = self.group_restraints
        self.group_restraints = None
        self.read(file=filename)
        self.group_restraints = self.loop_potential

        aln = self.read_loop_alignment(filename)
        self.create_loop_topology(aln, filename)
        atmsel = self._check_select_loop_atoms()
        self.loop_restraints(atmsel, aln)

        # Select corresponding restraints only:
        # only necessary to eliminate inefficiencies in 'special_restraints'
        # because MAKE_RSRS works with selected atoms now:
        self.restraints.unpick_all()
        self.restraints.pick(atmsel)
        self.restraints.condense()
        self.restraints.write(file=self.csrfile)

        # Calculate energy for the original (raw) loop:
        self.env.edat.nonbonded_sel_atoms = 1
        atmsel.energy()

        # Prepare the starting structure (comment it out if
        # the input PDB file is a better initial structure):
        self.build_ini_loop(atmsel)

        ini_model = "%s.IL%04d%04d%s" % (self.sequence, 0, num, self.pdb_ext)
        self.write(file=ini_model)

        sched = self.get_loop_schedule()
        if self.parallel_job is not None:
            self.parallel_loop_models(atmsel, ini_model, num, sched)
        else:
            self.multiple_loop_models(atmsel, ini_model, num, sched)
        self.group_restraints = oldgprsr

    def multiple_loop_models(self, atmsel, ini_model, num, sched):
        """Build all loop models for a given initial model"""
        for id1 in range(self.loop.starting_model, self.loop.ending_model + 1):
            r = self.single_loop_model(atmsel, ini_model, num, id1, sched)
            self.loop.outputs.append(r)

    def parallel_loop_models(self, atmsel, ini_model, num, sched):
        """Build all loop models for a given model using a parallel job"""
        from modeller.automodel.parallel import looptask

        job = self.parallel_job
        for id1 in range(self.loop.starting_model, self.loop.ending_model + 1):
            job.queue_task(looptask(self, atmsel, ini_model, num, id1, sched))
        self.loop.outputs.extend(job.run_all_tasks())

    def get_loop_schedule(self):
        """Get loop optimization schedule"""
        libsched = self.loop.library_schedule
        return libsched.make_for_model(self) * self.env.schedule_scale

    def new_loop_trace_file(self, id1, id2):
        """Open a new loop optimization trace file"""
        filename = modfile.default(file_ext='', file_id='.DL',
                                   root_name=self.sequence, id1=id1, id2=id2)
        return file(filename, 'w')

    def get_loop_actions(self):
        """Get actions to carry out during loop optimization"""
        return self.get_optimize_actions()

    def read_ini_loop_model(self, ini_model):
        """Read in the initial loop model"""
        # Make sure we read every atom, since when we write out the model, it
        # writes out every atom, not just non-HET/non-hydrogen/non-water
        io = io_data(copy=self.env.io)
        io.hetatm = io.water = io.hydrogen = True
        self.read(file=ini_model, io=io)

    def single_loop_model(self, atmsel, ini_model, num, id1, sched,
                          parallel=False):
        """Build a single loop model"""
        self.tracefile = self.new_loop_trace_file(id1, num)

        if parallel:
            self.read_top_par()
            self.read_ini_loop_model(ini_model)
            aln = self.read_loop_alignment(ini_model)
            self.create_loop_topology(aln, ini_model)
        else:
            self.read_ini_loop_model(ini_model)

        atmsel.randomize_xyz(deviation=5.0)

        if parallel:
            self.group_restraints = self.read_potential()
            self.rd_restraints()

        filename = self.get_loop_model_filename(self.sequence, id1, num,
                                                self.pdb_ext)
        out = {'name':filename, 'loopnum':id1, 'num':num, 'failure':None}

        actions = self.get_loop_actions()
        try:
            # Refine without the rest of the protein:
            self.env.edat.nonbonded_sel_atoms = 2
            self.optimize_loop(atmsel, sched, actions)
            # Refine in the context of the rest of the protein:
            self.env.edat.nonbonded_sel_atoms = 1
            self.optimize_loop(atmsel, sched, actions)

            (out['molpdf'], out['pdfterms']) = atmsel.energy()

            self.to_iupac()
        except (ModellerError, OverflowError), detail:
            if len(str(detail)) > 0:
                out['failure'] = detail
            else:
                out['failure'] = 'Optimization failed'
        else:
            self.loop_model_analysis(atmsel, ini_model, filename, out, id1, num)
        del self.tracefile
        return out

    def loop_model_analysis(self, atmsel, ini_model, filename, out, id1, num):
        """Energy evaluation and assessment, and write out the loop model"""
        self.user_after_single_loop_model()

        if self.loop.write_selection_only:
            self.select_loop_atoms().write(file=filename)
        else:
            self.write(file=filename)

        if self.accelrys:
            # Accelrys wants their analysis *after* the model is written, so
            # that the written-out model keeps the original template-derived
            # Biso, rather than getting the energy profile Biso:
            self.loop_model_analysis_accelrys(ini_model, id1, num)

        # Do model assessment if requested
        self.assess(atmsel, self.loop.assess_methods, out)

    def loop_model_analysis_accelrys(self, ini_model, id1, num):
        """Additional model analysis for Accelrys build"""
        # select all atoms
        allat = selection(self)
        # write total violations
        for (id, norm) in (('.EL', False), ('.NEL', True)):
            allat.energy(output='LONG ENERGY_PROFILE', normalize_profile=norm,
                         file=modfile.default(file_id=id, file_ext='',
                                              root_name=self.sequence,
                                              id1=id1, id2=num))

    def user_after_single_loop_model(self):
        """Used for any user analysis after building each loop model. Redefine
           as you see fit."""
        pass

    def build_ini_loop(self, atmsel):
        """Create the initial conformation of the loop. By default we place
           all atoms on a line between the loop termini, but you may want
           to use a different conformation, in which case you should redefine
           this routine. For example, if you want to leave the initial PDB
           file untouched, use a one-line 'pass' routine."""
        atmsel.unbuild()
        self.build(build_method='3D_INTERPOLATION', initialize_xyz=False)

    def build_charmm_loop_restraints(self, atmsel, rsr, aln):
        """Build loop restraints from CHARMM libraries"""
        dih_lib_only = True
        mnch_lib = 1
        res_at = 1
        for typ in ('bond', 'angle', 'improper', 'dihedral',
                    'phi-psi_binormal'):
            rsr.make(atmsel, aln=aln, restraint_type=typ,
                     spline_on_site=self.spline_on_site,
                     dih_lib_only=dih_lib_only,
                     mnch_lib=mnch_lib, restraint_sel_atoms=res_at)
        for typ in ('omega', 'chi1', 'chi2', 'chi3', 'chi4'):
            rsr.make(atmsel, aln=aln, restraint_type=typ+'_dihedral',
                     spline_on_site=self.spline_on_site,
                     dih_lib_only=dih_lib_only, mnch_lib=mnch_lib,
                     restraint_sel_atoms=res_at, spline_range=4.0,
                     spline_dx=0.3, spline_min_points=5)

    def loop_restraints(self, atmsel, aln):
        """Construct restraints for loop modeling"""
        rsr = self.restraints
        rsr.clear()
        if self.loop.library_restraints is not None:
            self.build_library_restraints(atmsel, rsr,
                                          self.loop.library_restraints)
        else:
            self.build_charmm_loop_restraints(atmsel, rsr, aln)
        self.special_restraints(aln)

    def _check_select_loop_atoms(self):
        """Get loop atoms to be refined, and check for sanity"""
        atmsel = self.select_loop_atoms()
        if not hasattr(atmsel, "get_atom_indices"):
            raise ModellerError("you must return a selection object " +
                                "from select_loop_atoms")
        elif len(atmsel) == 0:
            raise ModellerError("no atoms selected for loop refinement")
        elif atmsel.get_model() is not self:
            raise ModellerError("selection is defined on the wrong model")
        elif len(atmsel.only_no_topology()) > 0:
            raise ModellerError("some selected residues have no topology")
        else:
            print "%d atoms selected for loop refinement" % len(atmsel)
        return atmsel

    def optimize_loop(self, atmsel, sched, actions):
        """Optimize the geometry of a single loop model"""
        for step in sched:
            step.optimize(atmsel, max_iterations=self.loop.max_var_iterations,
                          output=self.optimize_output, min_atom_shift=0.001,
                          actions=actions)
        if self.loop.md_level:
            self.loop.md_level(atmsel, actions)
        cg = conjugate_gradients()
        cg.optimize(atmsel, max_iterations=1000, output=self.optimize_output,
                    min_atom_shift=0.00001, actions=actions)

    def select_loop_atoms(self):
        """The default loop atom selection routine. This selects all atoms near
           gaps in the alignment. You can redefine this routine to select a
           different region, and in fact this is necessary if you are refining a
           PDB file, as no alignment is available in this case."""
        if len(self.knowns) == 0:
            raise ModellerError("No alignment: you must redefine " +
                                "select_loop_atoms")
        aln = self.read_alignment()
        loops = self.loops(aln, minlength=5, maxlength=15, insertion_ext=2,
                           deletion_ext=1)
        sel = selection(loops).only_std_residues()
        if len(sel) == 0:
            raise ModellerError("No loops detected for refinement: " +
                                "you must redefine select_loop_atoms")
        return sel


class LoopData(modobject):
    """Methods and data unique to loop modeling"""

    starting_model = 1
    ending_model = 1
    write_selection_only = False
    md_level = None
    env = None
    outputs = None
    assess_methods = None
    library_restraints = None
    max_var_iterations = 200
    library_schedule = autosched.loop

    def __init__(self, env):
        self.md_level = refine.slow
        self.env = env.copy()
        self.env.edat.contact_shell = 7.0
        self.env.edat.dynamic_modeller = True
        self.env.edat.dynamic_sphere = True
        self.outputs = []

    def use_library_restraints(self, librestraints):
        """Use a restraints library"""
        self.library_restraints = librestraints
