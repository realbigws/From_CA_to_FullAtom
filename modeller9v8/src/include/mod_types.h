/** \file mod_types.h      Data types used by Modeller.
 *
 *             Part of MODELLER, Copyright(c) 1989-2008 Andrej Sali
 */

#ifndef MOD_TYPES_H
#define MOD_TYPES_H

#include <stdio.h>
#include <glib.h>

#include "fortran-pointer-types.h"
#include "mod_file.h"

/** The structures and defines below are automatically generated from their
    Fortran equivalents in the Modeller source code.
    Do not edit them manually! */

#define LENSPRM 150
#define MAXMOD 11
#define MAXSCH 50
#define MCNSRES 200
#define MAPRES 100
#define MDIHTYP 9
#define MRESCLS 22
#define MRESGRP 5
#define MRESINCLS 64
#define MAXDATATYP 14
#define NITER 10
#define NMEMTYP 6
#define LENF 600
#define LEN_ATMNAM 4
#define LEN_CHAIN 1
#define LEN_RESCHR 1
#define LEN_RESCRM 4
#define LEN_RESNUM 5
#define LEN_RESSTR 3
#define MAXCOD 40
#define MAXFLD 80
#define NRTYP 21

/** Classification of atoms into classes */
struct mod_atom_classes {
  /** Number of residue/atom types */
  int nattmod;
  /** Number of atom classes */
  int ngrpatm;
  /** Maximum number of residue/atom types */
  int maxattmod;
  /** Maximum number of atom classes */
  int maxgrpatm;
  /** Atom class for each residue/atom type; dimension(maxattmod) */
  struct mod_int1_array iatmgrp;
  /** Beginning of each class in atmmod and iatmgrp; dimension(maxgrpatm) */
  struct mod_int1_array iattmod;
  /** Names of residue/atom types belonging to a class; dim(2,maxattmod) */
  struct mod_char2_array attmod;
  /** Name of each atom class: dimension(maxgrpatm) */
  struct mod_char1_array grpatm;
};


/** Atomic coordinates and other data shared by template structures and models */
struct mod_coordinates {
  /** Number of real atoms */
  int natm;
  /** Number of real and pseudo atoms */
  int natmp;
  /** Maximum number of atoms (real and pseudo) */
  int maxatm;
  /** Residue -1 for each atom; dimension(maxatm) */
  struct mod_int1_array iresatm;
  /** Beginning atom -1 for each residue; dimension(maxres) */
  struct mod_int1_array iatmr1;
  /** Occupancy; dimension(maxatm) */
  struct mod_float1_array occ;
  /** Average sidechain Biso; dimension(maxatm) */
  struct mod_float1_array biso;
  /** Cartesian coordinates for real and pseudo atoms; dimension(maxatm) */
  struct mod_float1_array x, y, z;
  /** Atomic accessibility; dimension(maxatm) */
  struct mod_float1_array atmacc;
  /** IUPAC atom names; dimension(maxatm) */
  struct mod_char1_array atmnam;
  /** Atom element names; dimension(maxatm) */
  struct mod_char1_array element;
  /** Residue numbers (including insertion codes); dimension(maxres) */
  struct mod_char1_array resnum;
};


/** Residue sequence data shared by alignments and models */
struct mod_sequence {
  /** Number of residues in the sequence */
  int nres;
  /** Maximum number of residues in the sequence */
  int maxres;
  /** Number of segments in the sequence */
  int nseg;
  /** Maximum number of segments in the sequence */
  int maxseg;
  /** Segment names (chain id's where possible); dimension(maxseg) */
  struct mod_char1_array segid;
  /** Residue type; dimension(maxres) */
  struct mod_int1_array irestyp;
  /** First residue of segment i; dimension(maxseg) */
  struct mod_int1_array iress1;
  /** Last residue of segment i; dimension(maxseg) */
  struct mod_int1_array iress2;
  /** Resolution, -1 if it does not exist */
  float resol;
  /** R-factor, -1 if it does not exist */
  float rfactr;
  /** Protein name */
  char *name;
  /** Source organism */
  char *source;
  /** Type of entry (sequence, structure, structureN, structureX, structureN) */
  char *prottyp;
  /** First and last residue chain ID */
  char chnrng[2][LEN_RESCHR];
  /** First and last residue number */
  char rng[2][LEN_RESNUM];
};


/** Residue sequence data specific to alignments */
struct mod_alnsequence {
  /** TRUE if the code is read in */
  gboolean codein;
  /** Full atom filename obtained from atom_files or codes */
  char *atmfull;
  /** Filename to read PDB information from */
  char *atom_files;
  char *mematmfull[NMEMTYP];
  char memchnrng[2][NMEMTYP][LEN_CHAIN];
  char memrng[2][NMEMTYP][LEN_RESNUM];
  /** PDB code for the sequence */
  char codes[MAXCOD];
};


/** A type to hold all protein structure information */
struct mod_structure {
  /** TRUE if the protein has all data defined */
  gboolean accepts;
  /** TRUE if the structural information was read from an atom file */
  gboolean read_from_file;
  /** Coordinates */
  struct mod_coordinates cd;
  /** Bin -1 for dihedral angle classes; dimension(maxres,mdihtyp) */
  struct mod_int2_array idihc;
  /** Bin -1 for Wilmot mainchain conformation; dimension(maxres) */
  struct mod_int1_array imnchw;
  /** Atom indices of selected atom type 1; dimension(maxres) */
  struct mod_int1_array idsta1;
  /** Atom indices of selected atom type 2; dimension(maxres) */
  struct mod_int1_array idsta2;
  /** Number of neighbor residues of each residue; dimension(maxres) */
  struct mod_int1_array neigh;
  /** Residue indices of neighbors; dimension(maxngh,maxres) */
  struct mod_int2_array ineigh;
  /** Bin -1 for MODELLER secondary structure; dimension(maxres) */
  struct mod_int1_array isstruc;
  /** All dihedral angles for each residue; dimension(maxres,mdihtyp) */
  struct mod_float2_array dih;
  /** Residue accessibilities; dimension(maxres) */
  struct mod_float1_array acc;
  /** Mainchain curvature; dimension(maxres) */
  struct mod_float1_array curvn;
};


/** An alignment of protein sequences */
struct mod_alignment {
  /** Number of sequences in the alignment */
  int nseq;
  /** Maximum number of sequences in the alignment */
  int maxseq;
  /** Length of alignment */
  int naln;
  /** Maximum length of alignment */
  int maxaln;
  /** Number of comments in the aln file */
  int ncomment;
  /** Maximum number of comments in the aln file */
  int maxcomment;
  /** Residue index alignment; dimension(maxaln,maxseq) */
  struct mod_int2_array ialn;
  /** Biases of the comparison matrix at each alignment position; dim(maxaln) */
  struct mod_uchar1_array fix_position;
  /** Residue SS pred confidence index alignment; dimension(maxaln,maxseq) */
  struct mod_int2_array iconf;
  /** predicted secondary structure alignment; dimension(maxaln,maxseq) */
  struct mod_char2_array pred_ss;
  /** Inverted residue index alignment, gaps removed; dimension(maxaln,maxseq) */
  struct mod_int2_array invaln;
  /** Residue type alignment; dimension(maxaln,maxseq) */
  struct mod_char2_array caln;
  /** Comments; dimension(maxcomment) */
  struct mod_char1_array comment;
  /** Alignment accuracy; dimension(maxaln) */
  struct mod_float1_array alnacc;
  /** Alignment profile; dimension(maxaln,nprof) */
  struct mod_float2_array prof;
  /** Structural data for relevant sequences; dimension(maxseq) */
  struct mod_derv1_array struc;
  /** Residue sequence information; dimension(maxseq) */
  struct mod_derv1_array seq;
  /** Alignment-specific residue sequence information; dimension(maxseq) */
  struct mod_derv1_array alnseq;
  /** Fractional global similarity; dimension(nseq,nseq) when not null */
  struct mod_float2_array fractglsim;
  /** Scripting language object attached to this type (if any) */
  void *scriptobj;
  /** Pointers for maintaining a doubly-linked list */
  struct mod_derv_pt prev, next;
};


/** Information about pseudo atoms */
struct mod_pseudo_atoms {
  /** Number of real atoms defining all pseudo atoms */
  int npsdatm;
  /** Maximum number of real atoms defining all pseudo atoms */
  int maxpsdatm;
  /** Number of scalar intermediate results saved for derivatives calculation */
  int npsdsv;
  /** Maximum number of scalar intermediate results */
  int maxpsdsv;
  /** Number of all pseudo atoms (= natmp-natm) */
  int npseudo;
  /** Maximum number of all pseudo atoms */
  int maxpseudo;
  /** Intermediate results for derivatives; dimension(maxpsdsv) */
  struct mod_float1_array psdsv;
  /** Atom -1 of pseudo atom i; dimension(maxpseudo) */
  struct mod_int1_array iatpsd;
  /** Indices of real atoms defining a pseudoatm; dimension(maxpsdatm) */
  struct mod_int1_array indpsd;
  /** Starting -1 in vector of atom indices defining p. atom; dim(maxpseudo) */
  struct mod_int1_array ipsd;
  /** -1 in psdsv of the intermediate result for each p. atom; dim(maxpseudo) */
  struct mod_int1_array ipsdsv;
  /** Type of the pseudo atom i; dimension(maxpseudo) */
  struct mod_int1_array ipstyp;
  /** Number of real atoms defining the pseudo atom; dimension(maxpseudo) */
  struct mod_int1_array npsdef;
};


/** All EM density information */
struct mod_density {
  /** The resolution used with the function that calculates atomic density */
  float resolution;
  /** The voxel size of the density grid, also used for rendering the density */
  float vox_size;
  /** The density grid dimensions */
  int nx, ny, nz;
  /** The position on the map where to place the protein center of mass */
  float px, py, pz;
  /** Sum of squares over all grid points */
  float grid_sqr_sum;
  /** Filter type */
  int filter_type;
  /** The function to calculate the data (hard sphere, Gaussian...) */
  /** change the name of top_density_type for function_type!!!!! */
  int function_type;
  /** The cross-correlation coefficient between the map and the protein */
  float cc;
  /** Normalization factor */
  double norm_factor;
  /** The type of correlation coefficient between the map and the protein */
  int cc_func_type;
  /** Scripting language object attached to this type (if any) */
  void *scriptobj;
  /** Filter values (if needed) */
  struct mod_float1_array filter_values;
  /** The density grid; dimension(nz,ny,nx) */
  struct mod_float3_array grid;
  /** Pointers for maintaining a doubly-linked list */
  struct mod_derv_pt prev, next;
};


/**  */
/** SAXSDATA data type. This objects contains all relevant information */
/** for incorporation of SAXS constraints to modeling. In particular */
/** experimental SAXS data, formfactors, q, ... */
/** 12/08/05 FF - added w_s, wswitch, s_hybrid, normsq_exp */
/** last change 06/02/06 FF - added use_offset */
/** 02/26/07 FF - added scorespace, pr_smooth, autocorr */
/**  */
struct mod_saxsdata {
  /** ----- sampling of reciprocal space ----------------------------------- */
  /** number of sampling points */
  int ns;
  /** maximum number of sampling point */
  int maxs;
  /** number of different atoms (correlates to file!) */
  int natomtyp;
  /** number of formfactors per scattering center to be used. */
  /** 2 or larger if isomorphic replacement or Cysteine labeling */
  /** works... */
  int nscatts;
  /** maximum frequency in space (4\pi\sin(theta)/\lambda) in \AA */
  double s_max;
  /** minimum frequency in space */
  double s_min;
  /** sampling density in reciprocal space (s_max/(ns-1)) */
  double s_mesh;
  /** bandpass to exclude frequency below / above */
  double s_low, s_hi;
  /** frequency above which s^2 weighting is applied if 'hybrid' */
  /** weighting is specified */
  double s_hybrid;
  /** norm^2 of I(s) given the weighting scheme wswitch and the error */
  /** sigma_exp */
  double normsq_exp;
  /** scaling factor of model spectrum to match experimental one */
  double c;
  /** offset of experimental data (optional) */
  double offset;
  /** b-factor for Gaussian rolloff (optional) */
  double bfac;
  /** magnitude of Gaussian rolloff (optional) */
  double rolloff;
  /** chi_square of experimental and model saxs data */
  double chi_sq;
  /** electron density of solvent - default=0.334 e/A^3 (H2O) */
  double rho_solv;
  /** density in r */
  double dr;
  /** density in r of experimental data */
  double dr_exp;
  /** no of r-samples of experimental data p_r_exp */
  int nr_exp;
  /** no of sampling points for sinc function per 1 unit */
  int mesh_sinc;
  /** density of sampling for sinc function */
  double sinc_dens;
  /** array containing atom index for each atom (~formfactor file) */
  struct mod_int1_array atmindx;
  /** mesh in reciprocal space s(maxs) */
  struct mod_double1_array s;
  /** weighting factors in frequency space */
  struct mod_double1_array w_s;
  /** ----- intensity and radial distribution function --------------------- */
  /** from model - FF: later extend to multiple scattering factors - */
  /** i.e. isomorphous replacement, gold labeling */
  /** calculated SAXS profile of model intensity(maxs) */
  struct mod_double1_array intensity;
  /** calculated radial distribution function p(maxs) of model */
  struct mod_double1_array p_r;
  /** formfactors of scattering centers formfactor(natomtyp, maxs) */
  struct mod_double2_array formfactor;
  /** ----- intensity and radial distribution function - EXPERIMENT -------- */
  /** FF: later: extend to multiple f's */
  /** atom indices saxs restraint is working on */
  struct mod_int1_array isaxsatm;
  /** length of isaxsatm array */
  int n_isaxsatm;
  /** no of r-samples */
  int nr;
  /** measured saxsdata int_exp(maxs) */
  struct mod_double1_array int_exp;
  /** corresponding error */
  struct mod_double1_array sigma_exp;
  /** measured p(r) - derived from int_exp */
  struct mod_double1_array p_r_exp;
  /** radius mesh corresponding to p_r_exp */
  struct mod_double1_array r_exp;
  /** approximated error of p_r_exp - derived from int_exp */
  struct mod_double1_array p_r_sig;
  /** model p(r) on mesh of p_r_exp */
  struct mod_double1_array p_r_resamp;
  /** switch for using lookup tables for sinc and cos functions */
  gboolean use_lookup;
  /** switch for using additive constant in exp. data */
  gboolean use_offset;
  /** switch for using Gaussian rolloff on model */
  gboolean use_rolloff;
  /** switch for using nitrogen formfactor for convolution */
  gboolean use_conv;
  /** switch for mixed conformations */
  gboolean mixflag;
  /** switch smoothing of p_r */
  gboolean pr_smooth;
  /** include autocorrelation term in P(r)? */
  gboolean autocorr;
  /** lookup table for sinc function */
  struct mod_double1_array sinc_lookup;
  /** lookup table for cos function */
  struct mod_double1_array cos_lookup;
  /** Pointers for maintaining a doubly-linked list */
  struct mod_derv_pt prev, next;
  /** ----- scoring in reciprocal or real space 'reciprocal' 'real' */
  char scorespace[LENF];
  /** how should I(s) be computed? real space via P(r) or reciprocal? */
  char spaceflag[LENF];
  /** type of representation - 'allatm' or 'CA' */
  char represtyp[LENF];
  /** weighting scheme in frequency space */
  char wswitch[LENF];
  /** filename of used formfactors */
  char formfac_file[LENF];
};


/** The variables for dealing with symmetry restraints */
struct mod_symmetry {
  /** Number of atoms involved in all symmetry restraints */
  int natmsym;
  /** Maximum number of atoms involved in all symmetry restraints */
  int maxatmsym;
  /** Number of symmetry restraints (segments) */
  int nsegsym;
  /** Maximum number of symmetry restraints */
  int maxsegsym;
  /** Weight of the symmetry objective function term; dimension(maxatmsym) */
  struct mod_float1_array wghsym;
  /** Index into iatmsym1,2 of the first atom pair in each rsr; dim(maxsegsym) */
  struct mod_int1_array isegsym;
  /** First atom index; dimension(maxatmsym) */
  struct mod_int1_array iatsym1;
  /** Second atom index; dimension(maxatmsym) */
  struct mod_int1_array iatsym2;
};


/** Hash data for keeping track of excluded contacts */
struct mod_hash_contacts {
  /** The number of hash buckets */
  int nbuckets;
  /** Maximum number of hash buckets */
  int maxbuckets;
  /** The number of excluded pairs in each bucket; dimension(maxbuckets) */
  struct mod_int1_array nentries;
  /** Maximum number of excluded pairs per bucket */
  int maxentries;
  /** Atom indices of excluded pairs; dimension(maxentries,maxbuckets) */
  struct mod_int2_array ientry, jentry;
};


/** All schedule information (for variable target function optimization) */
struct mod_schedule {
  /** Physical restraint type scaling factors; dimension(maxphycns,maxsch) */
  struct mod_float2_array scaln;
  /** Default optimization method */
  int imethn[MAXSCH];
  /** The min and max residue separations for which restraints are calculated */
  int nrangn[2][MAXSCH];
  /** The number of steps in the schedule */
  int nsch;
  /** The current optimization step in the schedule */
  int step;
};


/** Topology information for a model */
struct mod_model_topology {
  /** Number of angles */
  int nang;
  /** Number of bonds */
  int nbnd;
  /** Number of dihedrals */
  int ndih;
  /** Number of impropers */
  int nimp;
  /** Maximum number of angles */
  int maxang;
  /** Maximum number of bonds */
  int maxbnd;
  /** Maximum number of dihedrals */
  int maxdih;
  /** Maximum number of impropers */
  int maximp;
  /** Angles in a model - dimension(3,maxang) */
  struct mod_int2_array iata;
  /** Bonds in a model - dimension(2,maxbnd) */
  struct mod_int2_array iatb;
  /** Dihedrals in a model - dimension(4,maxdih) */
  struct mod_int2_array iatd;
  /** Impropers in a model - dimension(4,maximp) */
  struct mod_int2_array iati;
  /** Number of atoms bonded to each atom - dimension(mdl%maxatm) */
  struct mod_int1_array natc;
  /** Atoms bonded to each atom - dimension(maxbat,mdl%maxatm) */
  struct mod_int2_array iatc;
};


/** A rigid body */
struct mod_rigid_body {
  /** Number of atoms */
  int natm;
  /** Maximum number of atoms in the body */
  int maxatm;
  /** Atoms in the body - dimension(maxatm) */
  struct mod_int1_array indatm;
  /** Position of center of mass (in world space) */
  float com[3];
  /** Current rotation, as a set of Euler angles to transform to world space */
  /** (CG) or a quaternion to transform to body space (MD) */
  float rot[4];
  /** Scaling factor for CG state and derivatives */
  float scalefact;
  /** Inverse scaling factor for CG state and derivatives */
  float invscalefact;
  /** Total body mass */
  float mass;
  /** Moment of inertia */
  float minert[3];
  /** Linear velocity */
  float linvel[3];
  /** Rotational (quaternion) velocity */
  float rotvel[4];
  /** Force on the center of mass */
  float force[3];
  /** Torque about the center of mass, in world space */
  float torque_world[3];
  /** Torque about the center of mass, in body space */
  float torque_body[3];
  /** Positions of atoms (in body space, rel. to COM) - dimension(3,maxatm) */
  struct mod_float2_array bodypos;
  /** Whether this body is picked for optimization */
  gboolean picked;
  /** Pointers to maintain a doubly-linked list */
  struct mod_derv_pt prev, next;
};


/** A list of Modeller restraints */
struct mod_restraints {
  /** Number of atoms in indatm() array */
  int natmcns;
  /** Maximum number of atoms in indatm() array */
  int maxatmcns;
  /** Current number of all static restraints */
  int ncsr;
  /** Maximum number of all static restraints */
  int maxcsr;
  /** Current number of restraint specs */
  int nspec;
  /** Maximum number of restraint specs */
  int maxspec;
  /** No. of excluded atom pairs (in bonds, angles, dih, & explicitly excl.) */
  int nexcl;
  /** Maximum no. of excluded atom pairs */
  int maxexcl;
  /** Current number of all selected static restraints */
  int nicsr;
  /** Current number of nonbond pairs */
  int npairs;
  /** Maximum number of nonbond pairs */
  int maxpairs;
  /** Number of parameters in pcsr() array */
  int nprmcns;
  /** Maximum number of parameters in pcsr() array */
  int maxprmcns;
  /** Minimum allowed residue span for calculation of dynamic non-bonded rsrs */
  int nrang1;
  /** Maximum allowed residue span for calculation of dynamic non-bonded rsrs */
  int nrang2;
  /** Number of global parameter arrays */
  int narray;
  /** Maximum number of global parameter arrays */
  int maxarray;
  /** Global parameter arrays; dimension(maxarray) */
  struct mod_derv1_array arrays;
  /** Parameters for all restraints as read from a file; dimension(maxprmcns) */
  struct mod_float1_array pcsr;
  /** Nonbond pair list; dimension(2,mpairs) */
  struct mod_int2_array iapairs;
  /** Starting position in indatm() for each restraint; dimension(maxcsr) */
  struct mod_int1_array iatm;
  /** Atom indices for each restraint; dimension(maxatmcns) */
  struct mod_int1_array indatm;
  /** Indices for all selected static rsrs; dimension(maxcsr) */
  struct mod_int1_array indcsr;
  /** Starting position in pcsr() for each restraint; dimension(maxcsr) */
  struct mod_int1_array ipcsr;
  /** Various restraint specs for each restraint; dimension(maxspec) */
  /**  */
  /** (1) math form (lower,upper,mono,poly,vdw,electr); */
  /** (2) polymodality of poly-modal Gaussians; */
  /** (3) feature type (dist, ang, dih); */
  /** (4) physical feature type (bond, angle, CA-CA dist...) */
  /** (5) number of atoms in restraint */
  /** (6) number of parameters in restraint */
  /** (7) number of atoms in first group */
  /** (8) future use */
  struct mod_int1_array itcsr;
  /** Starting position in itcsr() for each restraint; dimension(maxcsr) */
  struct mod_int1_array ispec;
  /** Excluded atom pairs; dimension(2,maxexcl) */
  struct mod_int2_array iexcl;
};


/** Parameters and atom classes/groups for dynamic_modeller restraints */
struct mod_group_restraints {
  /** Atom classes */
  struct mod_atom_classes atclass;
  /** Current number of restraint specs */
  int nspec;
  /** Maximum number of restraint specs */
  int maxspec;
  int nacsr;
  int nattacns;
  int maxattacns;
  int maxacsr;
  int nprmacns;
  int maxprmacns;
  /** dimension(maxprmacns) */
  struct mod_float1_array pacsr;
  struct mod_int1_array iarsr;
  struct mod_int2_array iarsr2;
  /** dimension(maxacsr) */
  struct mod_int1_array iatt;
  /** dimension(maxattacns) */
  struct mod_int1_array indatt;
  /** dimension(maxacsr) */
  struct mod_int1_array ipacsr;
  /** dimension(maxspec) */
  struct mod_int1_array itacsr;
  /** Starting position in itacsr() for each restraint; dimension(maxacsr) */
  struct mod_int1_array ispec;
  /** Scripting language object attached to this type (if any) */
  void *scriptobj;
  /** Pointers for maintaining a doubly-linked list */
  struct mod_derv_pt prev, next;
};


/** An internal coordinate, for building the model */
struct mod_intcoor {
  /** Bond length between first two atoms */
  float b1ic;
  /** Bond length between last two atoms */
  float b2ic;
  /** Torsion angle made by the four atoms */
  float pic;
  /** Bond angle between first three atoms */
  float t1ic;
  /** Bond angle between last three atoms */
  float t2ic;
  /** Atom indices defining the internal coordinate */
  int iar[4];
  /** Flag indicating that this is an improper torsion */
  gboolean tar;
};


/** All internal coordinates, as a vector */
struct mod_intcoor_vec {
  /** The internal coordinates themselves */
  struct mod_derv1_array d;
  /** The number of ICs in the vector */
  int dim;
};


/** \brief A protein model */
/**  */
/** Contains the model's topology, Cartesian coordinates, and properties */
/** of its atoms and residues. The model can consist of several segments, */
/** each of them built with a different residue topology library. Segments */
/** may or may not be connected. */
struct mod_model {
  /** Number of disulfides in model */
  int nss;
  /** Maximum number of disulfides in model */
  int maxss;
  /** File handle for debug output */
  struct mod_file *fh_deb;
  /** Sequence identity of closest template */
  float seq_id;
  /** Energy from last energy or optimization */
  float last_energy;
  /** PDB REMARK line(s) */
  GString *remark;
  /** Coordinates */
  struct mod_coordinates cd;
  /** Residue sequence information */
  struct mod_sequence seq;
  /** PDB file header */
  GString *header;
  /** Atomic charges; dimension(maxatm) */
  struct mod_float1_array charge;
  /** Velocities in dx/dt for real atoms; dimension(maxatm) */
  struct mod_float1_array vx, vy, vz;
  /** Derivatives for the selected atoms; dimension(maxatm) */
  struct mod_float1_array dvx, dvy, dvz;
  /** Position of atom type in ATGRPT classification; dimension(maxatm) */
  struct mod_int1_array iatta;
  /** Position of atom type in CHARMM topology file; dimension(maxatm) */
  struct mod_int1_array iattyp;
  /** For each residue, for each dihedral angle type, the atom indices of */
  /** a dihedral (idihres(1,i,j)=nundf, if undefined); */
  /** dimension(4,mdihtyp,maxres) */
  struct mod_int3_array idihres;
  /** Number of atomic neighbors within CONTACT_SHELL; dimension(maxatm) */
  struct mod_int1_array natngh;
  /** Residue indices of S-S bonds; dimension(2,maxss) */
  struct mod_int2_array iss;
  /** Pseudo atom information */
  struct mod_pseudo_atoms psd;
  /** Excluded pair information */
  struct mod_hash_contacts hashc;
  /** Internal coordinates, for building the model */
  struct mod_intcoor_vec intcr;
  /** All current schedule information */
  struct mod_schedule sched;
  /** All symmetry restraints for this model */
  struct mod_symmetry sym;
  /** All restraints for this model */
  struct mod_restraints rsr;
  /** All group restraints for this model */
  struct mod_derv_pt gprsr;
  /** All topology information for this model */
  struct mod_model_topology mtp;
  /** Rigid body definitions (if any) */
  struct mod_derv_pt rigbod;
  /** Scripting language object attached to this type (if any) */
  void *scriptobj;
  /** Pointers for maintaining a doubly-linked list */
  struct mod_derv_pt prev, next;
};


/** A database of sequences */
struct mod_sequence_db {
  /** Total number of amino acid residues in the database */
  gint64 nresdb;
  /** Maximum total number of amino acid residues in the database */
  gint64 maxresdb;
  /** Number of chains in the database */
  int nchn;
  /** Maximum number of chains in the database */
  int maxchn;
  /** Index of first residue in each chain; dimension(maxchn) */
  struct mod_int641_array iseqdb;
  /** Length of sequence in each chain; dimension(maxchn) */
  struct mod_int1_array nseqdb;
  /** Codes of sequences in the database; dimension(maxchn) */
  struct mod_char1_array codep;
  /** Type of sequence - only structure 'X' or sequence 'S'; dimension(maxchn) */
  struct mod_char1_array prottypdb;
  /** Residue types of the sequences in the database; dimension(maxresdb) */
  struct mod_uchar1_array seqdb;
  /** Resolutions of sequence (PIR format only); dimension(maxchn) */
  struct mod_float1_array resoldb;
  /** Scripting language object attached to this type (if any) */
  void *scriptobj;
  /** Pointers for maintaining a doubly-linked list */
  struct mod_derv_pt prev, next;
};


/** A sequence profile */
struct mod_profile {
  /** Number of sequences in the profile */
  int nseq;
  /** Number of positions in the profile */
  int npos;
  /** Maximum number of sequences in the profile */
  int maxseq;
  /** Maximum number of positions in the profile */
  int maxpos;
  /** Index number of the database sequence; dimension(maxseq) */
  struct mod_int1_array iprofile;
  /** The sequences in the profile; dimension(maxpos, maxseq) */
  struct mod_int2_array sprofile;
  /** Code of the database sequence; dimension(maxseq) */
  struct mod_char1_array code;
  /** Protein type of the sequence (X:structure,S:sequence); dimension(maxseq) */
  struct mod_char1_array prottyp;
  /** Iteration from which the alignment was derived; dimension(maxseq) */
  struct mod_int1_array iter;
  /** Beginning position of target in alignment; dimension(maxseq) */
  struct mod_int1_array i1res;
  /** End position of target in alignment; dimension(maxseq) */
  struct mod_int1_array i2res;
  /** Beginning position of dbseq in alignment; dimension(maxseq) */
  struct mod_int1_array j1res;
  /** End position of dbseq in alignment; dimension(maxseq) */
  struct mod_int1_array j2res;
  /** Number of equivalent positions; dimension(maxseq) */
  struct mod_int1_array neqv;
  /** Length of db sequence; dimension(maxseq) */
  struct mod_int1_array nres;
  /** Percentage sequence identity from alignment; dimension(maxseq) */
  struct mod_float1_array fid;
  /** E-value of the alignment; dimension(maxseq) */
  struct mod_double1_array evalue;
  /** Substitution matrix offset for local alignment */
  float matrix_offset;
  /** Number of iterations from profile building */
  int n_prof_iterations;
  /** Gap creation and extension penalties */
  float gap_penalties_1d[2];
  /** Residue-residue scoring file */
  char *rr_file;
  /** The file the profile was read from */
  char *filename;
  /** Scripting language object attached to this type (if any) */
  void *scriptobj;
  /** Pointers for maintaining a doubly-linked list */
  struct mod_derv_pt prev, next;
};


/** CHARMM residue type information */
struct mod_residue_typ {
  /** Number of atoms */
  int natm;
  /** Charges for all atoms */
  float charge[MAPRES];
  /** Total charge */
  float tchgrt;
  /** Distance, angle, dihedral, angle, distance for each IC entry */
  float ricrt[5][MCNSRES];
  /** Number of angles */
  int nang;
  /** Number of bonds */
  int nbnd;
  /** Number of dihedrals */
  int ndih;
  /** Number of IC records */
  int nic;
  /** Number of impropers */
  int nimp;
  /** CHARMM names */
  int resrt;
  /** True if residue read in from TOPOLOGY file */
  gboolean inrt;
  /** (I,+,-,#) IUPAC atom names for all bonds */
  char bnd[2][MCNSRES][LEN_ATMNAM+3];
  /** (I,+,-,#) IUPAC atom names for all angles */
  char ang[3][MCNSRES][LEN_ATMNAM+3];
  /** (I,+,-,#) IUPAC atom names for all dihedrals */
  char dih[4][MCNSRES][LEN_ATMNAM+3];
  /** (I,+,-,#) IUPAC atom names for all impropers */
  char imp[4][MCNSRES][LEN_ATMNAM+3];
  /** (I,+,-,#) IUPAC atom names for all IC entries */
  char icart[4][MCNSRES][LEN_ATMNAM+3];
  /** IUPAC atom names (2 extra chars needed for patch 1: or 2: prefix) */
  char atiuprt[MAPRES][LEN_ATMNAM+2];
  /** CHARMM atom names */
  char atresrt[MAPRES][LEN_ATMNAM];
  /** CHARMM atom names for all topology models */
  char atmmod[MAPRES][MAXMOD][LEN_ATMNAM];
  /** N and C patch residues */
  char ptchrt[2][LEN_RESCRM];
};


/** A patching residue */
struct mod_patch_type {
  /** Residue type information */
  struct mod_residue_typ restyp;
  /** Number of atoms deleted from a patched residue */
  int ndelart;
  /** Number of bonds deleted from a patched residue */
  int ndelbrt;
  /** Number of dihedrals deleted from a patched residue */
  int ndeldrt;
  /** Number of angles deleted from a patched residue */
  int ndelgrt;
  /** Number of impropers deleted from a patched residue */
  int ndelirt;
  /** Atoms to be deleted from the patched residue */
  int delart[MAPRES];
  /** Bonds to be deleted from the patched residue */
  char delbrt[2][MCNSRES][LEN_ATMNAM+3];
  /** Dihedrals to be deleted from the patched residue */
  char deldrt[4][MCNSRES][LEN_ATMNAM+3];
  /** Angles to be deleted from the patched residue */
  char delgrt[3][MCNSRES][LEN_ATMNAM+3];
  /** Impropers to be deleted from the patched residue */
  char delirt[4][MCNSRES][LEN_ATMNAM+3];
};


/** Information from a residue topology file */
struct mod_topology {
  /** The submodel (old TOPOLOGY_MODEL) of this topology (for radii etc.) */
  int submodel;
  /** Number of CHARMM residue types */
  int nresrt;
  /** Number of  residue topology libs in memory */
  int nrtl;
  /** Number of sub-topology models */
  int ntopmod;
  /** Whether to auto-generate angles */
  gboolean autang;
  /** Whether to auto-generate dihedrals */
  gboolean autdih;
  /** Number of CHARMM patches */
  int npatch;
  /** Maximum number of CHARMM patches */
  int maxpatch;
  /** CHARMM residue types; dimension(maxrestyp) */
  struct mod_derv1_array restyp;
  /** CHARMM patch types; dimension(maxpatch) */
  struct mod_derv1_array patchtyp;
  /** Mass (in atomic units) of each CHARMM atom type; dimension(maxatmtyp) */
  struct mod_float1_array amrt;
  /** Element name of each CHARMM atom type; dimension(maxatmtyp) */
  struct mod_char1_array element;
  /** Mapping from patch names to indices */
  void *patch_names;
};


/** Information from a parameter file */
struct mod_parameters {
  /** Number of bonds in the file */
  int nbnd;
  /** Maximum number of bonds in the file */
  int maxbnd;
  /** Number of angles in the file */
  int nang;
  /** Maximum number of angles in the file */
  int maxang;
  /** Number of dihedrals in the file */
  int ndih;
  /** Maximum number of dihedrals in the file */
  int maxdih;
  /** Number of impropers in the file */
  int nimp;
  /** Maximum number of impropers in the file */
  int maximp;
  /** IUPAC atom names of the bond; dimension(2,maxbnd) */
  struct mod_char2_array bnd;
  /** IUPAC atom names of the angle; dimension(3,maxang) */
  struct mod_char2_array ang;
  /** IUPAC atom names of the dihedral; dimension(4,maxdih) */
  struct mod_char2_array dih;
  /** IUPAC atom names of the improper; dimension(4,maximp) */
  struct mod_char2_array imp;
  /** Force constant for each bond; dimension(maxbnd) */
  struct mod_float1_array fbnd;
  /** Force constant for each angle; dimension(maxang) */
  struct mod_float1_array fang;
  /** Force constant for each dihedral; dimension(maxdih) */
  struct mod_float1_array fdih;
  /** Force constant for each improper; dimension(maximp) */
  struct mod_float1_array fimp;
  /** Mean value for each bond; dimension(maxbnd) */
  struct mod_float1_array abnd;
  /** Mean value for each angle; dimension(maxang) */
  struct mod_float1_array aang;
  /** Mean value for each dihedral; dimension(maxdih) */
  struct mod_float1_array adih;
  /** Mean value for each improper; dimension(maximp) */
  struct mod_float1_array aimp;
  /** Multiplicity for each dihedral; dimension(maxdih) */
  struct mod_int1_array idih;
  /** Multiplicity for each improper; dimension(maximp) */
  struct mod_int1_array iimp;
  /** Energy well depth for Lennard-Jones interaction; dimension(maxatmtyp) */
  struct mod_float1_array eminat;
  /** Energy well depth for Lennard-Jones interaction, 2nd set; dim(maxatmtyp) */
  struct mod_float1_array eminat2;
  /** r/2 for Lennard-Jones interaction; dimension(maxatmtyp) */
  struct mod_float1_array rminat;
  /** r/2 for Lennard-Jones interaction, 2nd set; dimension(maxatmtyp) */
  struct mod_float1_array rminat2;
};


/** Random number generator state */
struct mod_randstate {
  /** Random seed */
  int seed;
  /** 'ran3' state */
  int ran3_inext, ran3_inextp, ran3_iff, ran3_ma[55];
  /** 'iran1' state */
  int iran1_iy, iran1_iff, iran1_ir[97];
};


/** Wilmot phi/psi mainchain conformation classes */
struct mod_wilmot_data {
  /** Bin indices for Wilmot conform for Phi,Psi areas 10x10^o */
  int iwilmot[36][36][NRTYP];
  /** Mainchain class characters in IWILMOT() */
  char mnchcls[64];
};


/** Residue dihedral angle type definitions */
struct mod_resdih_data {
  /** Number of dihedral angle types */
  int ndihc;
  /** Maximal number of dihedral angle optima */
  int ndiht;
  /** Dihedral angle names */
  char dihnam[MDIHTYP][LEN_RESSTR+1];
};


/** Information from residue, dihedral, topology, parameter libraries */
struct mod_libraries {
  /** Number of residue types */
  int nrestyp;
  /** Maximum number of residue types */
  int maxrestyp;
  /** Number of atom types */
  int natmtyp;
  /** Maximum number of atom types */
  int maxatmtyp;
  /** Maximum number of residue dihedral angle optima */
  int maxdihopt;
  /** Number of secondary structure templates */
  int nsstyp;
  /** Maximum number of secondary structure templates */
  int maxsstyp;
  /** Maximum number of residues in a secondary structure template */
  int maxresss;
  /** Number of residues in a secondary structure template; dim(maxsstyp) */
  struct mod_int1_array nresss;
  /** Secondary structure type of a SS template; dimension(maxsstyp) */
  struct mod_int1_array itypss;
  /** Distances in sec. structure template; dim(maxresss,maxresss,maxsstyp) */
  struct mod_float3_array dstss;
  /** DRMS and MAX(delta distance) cutoffs for sec. structure template; */
  /** dimension(6,maxsstyp) */
  struct mod_float2_array sscut;
  /** BLK residue type */
  int iblktyp;
  /** Gap residue type */
  int igaptyp;
  /** Water residue type */
  int iwattyp;
  /** Random number state */
  struct mod_randstate rstat;
  /** Three letter (PDB/IUPAC) residue codes; dimension(maxrestyp) */
  struct mod_char1_array brkres;
  /** Four letter (CHARMM) residue codes; dimension(maxrestyp) */
  struct mod_char1_array chmres;
  /** One letter residue codes; dimension(maxrestyp) */
  struct mod_char1_array modres;
  /** HETATM indicator; dimension(maxrestyp) */
  struct mod_bool1_array hetatm;
  /** Mapping from all residue types to standard ones; dimension(maxrestyp) */
  struct mod_uchar1_array stdres;
  /** Default residue type to map non-standard types to */
  int default_stdres;
  /** Mapping table for one letter codes */
  int modmap[256];
  /** Number of mainchain conformation classes */
  int nmnch;
  /** Wilmot phi/psi mainchain conformation classes */
  struct mod_wilmot_data wilmot;
  /** Residue dihedral angle type definitions */
  struct mod_resdih_data resdih;
  /** Mainchain conformation classes; phi/psi means; dim(maxmnch,2,maxrestyp) */
  struct mod_float3_array amnch;
  /** Mainchain conformation classes; phi/psi stdevs; dim(maxmnch,2,maxrestyp) */
  struct mod_float3_array smnch;
  /** Mainchain conf. classes; cross-corr. coeff.; dim(maxmnch,maxrestyp) */
  struct mod_float2_array rmnch;
  /** Mainchain conformation weights; dimension(maxmnch,maxrestyp) */
  struct mod_float2_array wmnch;
  /** Residue-residue comparison matrix; dimension(maxrestyp,maxrestyp) */
  struct mod_float2_array rrwght;
  /** Dihedral library optimum mean; dimension(maxdihopt,mdihtyp,maxrestyp) */
  struct mod_float3_array adihlib;
  /** Dihedral library optimum angle range; dim(2,maxdihopt,mdihtyp,maxrestyp) */
  struct mod_float4_array dihcls;
  /** Dihedral library optimum stdev; dim(maxdihopt,mdihtyp,maxrestyp) */
  struct mod_float3_array sdihlib;
  /** Dihedral library optimum weight; dim(maxdihopt,mdihtyp,maxrestyp) */
  struct mod_float3_array wdihlib;
  /** Number of dihedrals in the dihedral library; dim(maxrestyp) */
  struct mod_int1_array ndihlib;
  /** Number of optima for each dihedral library entry; dim(mdihtyp,maxrestyp) */
  struct mod_int2_array ndihopt;
  /** Atom names for each dihedral library entry; dim(4,mdihtyp,maxrestyp) */
  struct mod_char3_array dihlib;
  /** van der Waals radii for each atom type in each topology submodel; */
  /** dimension(maxatmtyp,maxmod) */
  struct mod_float2_array vdwcnt;
  /** 1-4 van der Waals radii for each atom type in each topology submodel; */
  /** dimension(maxatmtyp,maxmod) */
  struct mod_float2_array vdwcnt14;
  /** CHARMM atom names for each atom type; dimension(maxatmtyp) */
  struct mod_char1_array chmatm;
  /** Number of topology submodels defined for VDW radii */
  int nvdwmod;
  /** Number of residue groupings that can be selected in MDT */
  int nresgrp;
  /** Number of classes within each residue group */
  int nrescls[MRESGRP];
  /** The topology library */
  struct mod_topology tpl;
  /** The parameter library */
  struct mod_parameters prm;
  /** Scripting language object attached to this type (if any) */
  void *scriptobj;
  /** Pointers for maintaining a doubly-linked list */
  struct mod_derv_pt prev, next;
  /** One-letter residue types for each residue group/class */
  char resgrp[MRESCLS][MRESGRP][MRESINCLS];
};


/** Controls reading from atom files */
struct mod_io_data {
  /** Whether to read HETATM records */
  gboolean hetatm;
  /** Whether to read hydrogen atoms */
  gboolean hydrogen;
  /** Whether to read water residues */
  gboolean water;
  /** Directories in which to look for atom files */
  char **atom_files_directory;
  /** Pointers for maintaining a doubly-linked list */
  struct mod_derv_pt prev, next;
};


/** Linked list of pointers to saxsdata structures */
struct mod_saxslist {
  /** SAXS structure pointer */
  struct mod_derv_pt d;
  /** Next structure in list */
  struct mod_derv_pt next;
};


/** Variables used during evaluation of the objective function */
struct mod_energy_data {
  /** User-defined energy term callbacks */
  void *user_terms;
  /** User-defined local similarity function */
  void *user_simloc;
  /** User-defined distance similarity function */
  void *user_distsimloc;
  /** The range for Lennard-Jones interaction smoothing to 0 */
  float lennard_jones_switch[2];
  /** The range for Coulomb interaction smoothing to 0 */
  float coulomb_switch[2];
  /** Distance cutoff for calculation of the non-bonded pairs list */
  float contact_shell;
  /** Relative dielectric */
  float relative_dielectric;
  /** Factor for van der Waals radii */
  float radii_factor;
  /** Whether to do dynamic pairs irrespective of anything */
  gboolean dynamic_pairs;
  /** Whether to use dynamic soft-sphere repulsion terms */
  gboolean dynamic_sphere;
  /** Whether to use dynamic Coulomb energy terms */
  gboolean dynamic_coulomb;
  /** Whether to use dynamic Lennard-Jones energy terms */
  gboolean dynamic_lennard;
  /** Whether to use dynamic MODELLER (statistical) non-bonded restraints */
  gboolean dynamic_modeller;
  /** Whether to use dynamic accessibility energy terms */
  gboolean dynamic_access;
  /** When to update the pairs list */
  float update_dynamic;
  /** Only include nonbond pairs that have at least this many selected atoms */
  int nonbonded_sel_atoms;
  /** Whether to consider SG-SG bonds similarly to the polypeptide chain */
  gboolean covalent_cys;
  /** Whether to exclude bonds/angles/diheds/excl pairs from distance rsrs */
  gboolean excl_local[4];
  /** Number of residues at which to begin using NlogN nonbond pairs routine */
  int nlogn_use;
  /** Standard deviation of soft-sphere repulsion */
  float sphere_stdev;
  /** Coulombic interaction conversion factor */
  float coulcnv;
  /** Electrostatics switching parameters */
  float relct3, relct6;
  /** Lennard-Jones switching parameters */
  float rlnrd3, rlnrd6;
  /** Accumulative change in positions that triggers dynamical violations */
  /** update (in Angstroms) */
  float shifttot;
  /** Maximal range for statistical pairwise potential */
  float statpotmax;
  /** Counters for various operations during optimization */
  int iter[NITER];
  /** Density to fit the model against (if not null) */
  struct mod_derv_pt den;
  /** List of SAXS data to fit (if not null) */
  struct mod_derv_pt saxslst;
  /** Lennard-Jones Aij and Bij parameters for each atom type pair; */
  /** dim(natmtyp,natmtyp) */
  struct mod_float2_array aij, bij;
  /** Scaled van der Waals radii for each atom type; dim(natmtyp) */
  struct mod_float1_array svdwcnt;
  /** Pointers for maintaining a doubly-linked list */
  struct mod_derv_pt prev, next;
};


/** Optimizer state information */
struct mod_optimizer {
  /** Scripting language object attached to this type (if any) */
  void *scriptobj;
  /** Any periodic actions to carry out */
  void *actions;
  /** Current step */
  int step;
  /** Number of objective function evaluations */
  int funcs;
  /** Initial energy */
  float init_e;
  /** Current energy */
  float current_e;
  /** Average shift of atoms in angstroms */
  float shiftavr;
  /** Maximal shift of atoms in angstroms */
  float shiftmax;
};


/** All selected real and pseudo atoms, rigid or non-rigid */
struct mod_energy_selection {
  /** Total number of selected atoms (real and pseudo) */
  int n_inds;
  /** Number of selected pseudo atoms (always at the end of the selection) */
  int n_selpsd;
  /** Selected atom indices */
  struct mod_int1_array inds;
  /** Indicator for whether each real or pseudo atom in the model is selected */
  struct mod_bool1_array picatm;
};


/** All selected real (non-pseudo) atoms and rigid bodies, for moving */
struct mod_optimize_selection {
  /** Number of selected real non-rigid atoms */
  int n_inds;
  /** Number of selected rigid bodies */
  int n_selbody;
  /** Number of degrees of freedom */
  int ndegf;
  /** Selected real non-rigid atom indices */
  struct mod_int1_array inds;
};


/** CG optimizer state information */
struct mod_cg_optimizer {
  /** Base optimizer type */
  struct mod_optimizer opt;
  /** Current gradient */
  float gsq;
  /** Pointers for maintaining a doubly-linked list */
  struct mod_derv_pt prev, next;
};


/** MD optimizer state information */
struct mod_md_optimizer {
  /** Base optimizer type */
  struct mod_optimizer opt;
  /** Average energy */
  float aver_e;
  /** Highest energy */
  float max_e;
  /** Lowest energy */
  float min_e;
  /** Standard deviation of the energy */
  float stdev_e;
  /** Average atomic shift (along 1 axis) */
  float smdavr;
  /** Kinetic energy */
  float kinetic_e;
  /** Kinetic temperature */
  float kinetic_temp;
  /** Set temperature */
  float temperature;
  /** Highest energy step during MD */
  int imdmax;
  /** Lowest energy step during MD */
  int imdmin;
  /** Guiding factor for self-guided MD */
  float guide_factor;
  /** Guiding time for self-guided MD */
  float guide_time;
  /** Friction factor for Langevin dynamics */
  float friction;
  /** Model coordinates with the lowest energy during MD; dimension(natm) */
  struct mod_float1_array xmin, ymin, zmin;
  /** Guiding forces for self-guided dynamics; dimension(natm) */
  struct mod_float1_array guide_x, guide_y, guide_z;
  /** Temporary storage for random forces during SGLD; dimension(3,natm) */
  struct mod_float2_array randf;
  /** Pointers for maintaining a doubly-linked list */
  struct mod_derv_pt prev, next;
};


/** Quasi Newton optimizer state information */
struct mod_qn_optimizer {
  /** Base optimizer type */
  struct mod_optimizer opt;
  /** Pointers for maintaining a doubly-linked list */
  struct mod_derv_pt prev, next;
};


/** State optimizer state information */
struct mod_state_optimizer {
  /** Base optimizer type */
  struct mod_optimizer opt;
  /** Model to be optimized */
  struct mod_derv_pt mdl;
  /** Nonbond pair function */
  int ipair;
  /** Number of physical restraint types */
  int nphycns;
  /** Selected atoms for optimization */
  struct mod_optimize_selection opsel;
  /** Selected atoms for energy calculation */
  struct mod_energy_selection ensel;
  /** Current state vector; dimension(nvar) */
  struct mod_float1_array state;
  /** Energy function scaling terms; dimension(nphycns) */
  struct mod_float1_array scaln;
  /** Pointers for maintaining a doubly-linked list */
  struct mod_derv_pt prev, next;
};


/** A type to hold a selection of one or more word options in a TOP variable */
struct mod_optvar {
  /** Whether 'given' can contain only one option */
  gboolean single;
  /** \internal The index of 'given' within 'possibles' */
  int numopt;
  /** The selection of options to choose from (blank separated) */
  char possibles[LENF];
  /** Option(s) selected by the user (blank separated) */
  char given[LENF];
};


/** PSSM data type. This object holds a the number positions in the PSSM */
/** (npos), the filename of the actual sequence profile the PSSM was */
/** derived from (fname), the actual position-specific scoring matrix */
/** itself (wmat). */
struct mod_pssmobj {
  /** Number of positions in the pssm */
  int npos;
  /** The pssm itself */
  struct mod_short2_array wmat;
  /** Name of the file that corresponds to this pssm */
  char fname[LENF];
};


/** PSSM Database data type. */
/** The database contains an array of pointers to pssm objects. */
/** The ipssm(:) refers to the ith pssm in the db object. */
struct mod_pssmdbobj {
  /** Number of pssms in the object. */
  int npssm;
  struct mod_derv1_array ipssm;
  /** Scripting language object attached to this type (if any) */
  void *scriptobj;
  /** Pointers for maintaining a doubly-linked list */
  struct mod_derv_pt prev, next;
};


/** Storage for the neighbor list used in GB/SA. */
struct mod_gbsa_neighbor {
  /** Atom indices of all close contacts; dim(natm,maxcontact) */
  struct mod_int2_array contact;
  /** Number of close contacts for each atom; dim(natm) */
  struct mod_int1_array ncontact;
  /** Number of atoms in the system */
  int natm;
  /** Maximum number of contacts for any atom */
  int maxcontact;
  /** Self-pointer to allow us to free this object from a C pointer */
  struct mod_derv_pt self;
};


/** Library of feature data used by MDTs */
struct mod_mdt_library {
  /** Number of library features */
  int nfeat;
  /** Allocated number of library features */
  int maxfeat;
  /** MDT library features */
  struct mod_mdt_libfeature *features;
  /** Currently unused */
  int unused1;
  /** Currently unused */
  int unused2;
  /** Currently unused */
  char unused3[2][LEN_ATMNAM];
};


/** Most MDT data */
struct mod_mdt {
  /** Number of elements in the array */
  int nelems;
  /** Maximum number of elements in the array */
  int maxelems;
  /** Used for handling specialties of the MDT in mdt.F, plot.F, and MODELLER */
  int nbinx;
  /** Number of features in this MDT (dimensionality) */
  int nfeat;
  /** Maximum number of features in this MDT */
  int maxfeat;
  /** MDT features */
  struct mod_mdt_feature *features;
  /** Used for handling specialties of the MDT */
  int nfeat1;
  /** TRUE if data file/action i be done in precalc() */
  gboolean readin[MAXDATATYP];
  /** Number of proteins in comparison (1,2,3) */
  int nprotcmp;
  /** Type of storage for the MDT bins (float, double, etc.) */
  int bin_type;
  /** The MDT array itself (access with mod_mdt_bin_get/set); */
  /** dimension(maxelems) */
  char *bindata;
};


#endif /* MOD_TYPES_H */
