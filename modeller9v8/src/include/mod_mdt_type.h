/** \file mod_mdt_type.h   MDT types.
 *
 *             Part of MODELLER, Copyright(c) 1989-2010 Andrej Sali
 */

#ifndef MOD_MDTTYPE_H
#define MOD_MDTTYPE_H

#include <glib.h>
#include "mod_types.h"

#ifdef __cplusplus
extern "C" {
#endif

/** Data file types needed to be read in before evaluating a feature */
typedef enum {
  MOD_MDTF_STRUCTURE = 1,
  MOD_MDTF_SECSTRUC,
  MOD_MDTF_PSA,
  MOD_MDTF_NEIGHBORS,
  MOD_MDTF_CONTACT,
  MOD_MDTF_PHI,
  MOD_MDTF_PSI,
  MOD_MDTF_OMEGA,
  MOD_MDTF_CHI1,
  MOD_MDTF_CHI2,
  MOD_MDTF_CHI3,
  MOD_MDTF_CHI4,
  MOD_MDTF_DIHEDRALS,
  MOD_MDTF_CHI5
} mod_mdt_file;

/** Property to be calculated before evaluating a feature */
typedef enum {
  /* obsolete properties removed */
  MOD_MDTC_NONE = 9,
  /* obsolete properties removed */
  MOD_MDTC_MNRAMA = 12,    /* mainchain conformation (Ramachandran) */
  /* obsolete properties removed */
  MOD_MDTC_CHI1CL = 18,    /* chi1 class */
  MOD_MDTC_CHI2CL,
  MOD_MDTC_CHI3CL,
  MOD_MDTC_CHI4CL,
  MOD_MDTC_OMEGA,     /* obsolete */
  MOD_MDTC_PHICL,
  MOD_MDTC_PSICL,
  MOD_MDTC_OMEGACL,
  MOD_MDTC_CHI5,      /* obsolete */
  MOD_MDTC_CHI5CL
} mod_mdt_calc;

/** Type of scan over available data for a feature */
typedef enum {
  MOD_MDTS_PROTEIN = 0,
  MOD_MDTS_RESIDUE,
  MOD_MDTS_RESIDUE_PAIR,
  MOD_MDTS_ATOM,           /* any atom (first entry in alignment only) */
  MOD_MDTS_ATOM_PAIR,
  MOD_MDTS_TUPLE,          /* atom tuple (first entry in alignment only) */
  MOD_MDTS_TUPLE_PAIR,
  MOD_MDTS_BOND,
  MOD_MDTS_ANGLE,
  MOD_MDTS_DIHEDRAL
} mod_mdt_scan;

/** Protein(s) on which a given feature acts */
typedef enum {
  MOD_MDTP_A = 1,
  MOD_MDTP_B,
  MOD_MDTP_AB,
  MOD_MDTP_C,
  MOD_MDTP_AC
} mod_mdt_protein;

/** Type of storage for the MDT bins */
typedef enum {
  MOD_MDTB_FLOAT = 1,
  MOD_MDTB_DOUBLE,
  MOD_MDTB_INT32,
  MOD_MDTB_UINT32,
  MOD_MDTB_INT16,
  MOD_MDTB_UINT16,
  MOD_MDTB_INT8,
  MOD_MDTB_UINT8
} mod_mdt_bin_type;

/** MDT bin type */
struct mod_mdt_bin {
  /** Start and end of bin range */
  float rang1, rang2;
  /** Human-readable symbol */
  char *symbol;
};

/** MDT library feature */
struct mod_mdt_libfeature {
  /** Indices of data files/actions to be done for each feature in precalc() */
  int *idatfil;
  /** Proteins to read data from (1=protein 1, 2=protein 2, 3=proteins 1 and 2,
      4=protein 3, 5=proteins 1 and 3) */
  int iknown;
  /** Type of scan */
  mod_mdt_scan iresfeat;
  /** Type of precalculation required */
  int idatfeat;
  /** >1 if asymmetric (then N^2 loop, not (N*(N-1)/2)) */
  int isymm;
  /** The number of actions in precalc */
  int ndatfil;
  /** Feature name */
  char *name;
  /** The type of the feature in mdt.bin:
      0 ... no symbols, no ranges specified (hardwired);
      1 ... no ranges specified, only symbols;
      2 ... real ranges and symbols specified explicitly;
      3 ... real ranges specified with beginning and width (no symbols) */
  int itsymb;
  /** Number of bins (including undefined) */
  int nbins;
  /** Bin data; dimension(nbins) */
  struct mod_mdt_bin *bins;
};

/** MDT feature */
struct mod_mdt_feature {
  /** Starting bin -1 for each selected feature */
  int istart;
  /** Ending bin -1 for each selected feature */
  int iend;
  /** Number of used bins for each feature */
  int nbins;
  /** Integer type of each feature */
  int ifeat;
  /** Number of positions to skip for lookup of bin indices */
  int stride;
};

struct mod_mdt;
struct mod_mdt_library;

/** Register a new MDT library feature */
void mod_mdt_libfeature_register(struct mod_mdt_library *mlib, int ifeat,
                                 const char *name, mod_mdt_calc precalc_type,
                                 mod_mdt_protein protein_type,
                                 mod_mdt_scan scan_type, gboolean asymmetric,
                                 ...);

/** Redimension a single MDT library feature. */
void mod_mdt_libfeature_nbins_set(struct mod_mdt_libfeature *feat, int nbins);

/** Set the number of bin elements in an MDT. */
void mod_mdt_nelems_set(struct mod_mdt *mdt, int nelems);

/** Set the number of features in an MDT. */
void mod_mdt_nfeat_set(struct mod_mdt *mdt, int nfeat);

/** Size in bytes of a single bin value */
size_t mod_mdt_bin_get_size(const struct mod_mdt *mdt);

/** Get a bin value */
static inline double mod_mdt_bin_get(const struct mod_mdt *mdt, int offset)
{
  switch (mdt->bin_type) {
  case MOD_MDTB_FLOAT:
    return (double)(((float *)mdt->bindata)[offset]);
  case MOD_MDTB_DOUBLE:
    return ((double *)mdt->bindata)[offset];
  case MOD_MDTB_INT32:
    return (double)(((gint32 *)mdt->bindata)[offset]);
  case MOD_MDTB_UINT32:
    return (double)(((guint32 *)mdt->bindata)[offset]);
  case MOD_MDTB_INT16:
    return (double)(((gint16 *)mdt->bindata)[offset]);
  case MOD_MDTB_UINT16:
    return (double)(((guint16 *)mdt->bindata)[offset]);
  case MOD_MDTB_INT8:
    return (double)(((gint8 *)mdt->bindata)[offset]);
  case MOD_MDTB_UINT8:
    return (double)(((guint8 *)mdt->bindata)[offset]);
  }
  return 0.0;
}

/** Set a bin value */
static inline void mod_mdt_bin_set(struct mod_mdt *mdt, int offset, double val)
{
  switch (mdt->bin_type) {
  case MOD_MDTB_FLOAT:
    ((float *)mdt->bindata)[offset] = (float)val;
    break;
  case MOD_MDTB_DOUBLE:
    ((double *)mdt->bindata)[offset] = val;
    break;
  case MOD_MDTB_INT32:
    ((gint32 *)mdt->bindata)[offset] = (gint32)val;
    break;
  case MOD_MDTB_UINT32:
    ((guint32 *)mdt->bindata)[offset] = (guint32)val;
    break;
  case MOD_MDTB_INT16:
    ((gint16 *)mdt->bindata)[offset] = (gint16)val;
    break;
  case MOD_MDTB_UINT16:
    ((guint16 *)mdt->bindata)[offset] = (guint16)val;
    break;
  case MOD_MDTB_INT8:
    ((gint8 *)mdt->bindata)[offset] = (gint8)val;
    break;
  case MOD_MDTB_UINT8:
    ((guint8 *)mdt->bindata)[offset] = (guint8)val;
    break;
  }
}

#ifdef __cplusplus
}
#endif

#endif   /* MOD_MDTTYPE_H */
