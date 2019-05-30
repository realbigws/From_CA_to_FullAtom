/** \file mod_hdf5.h       Utility functions for accessing the HDF5 library.
 *
 *             Part of MODELLER, Copyright(c) 1989-2010 Andrej Sali
 */

#ifndef MOD_HDF5_H
#define MOD_HDF5_H

#include "hdf5.h"
#include "hdf5_hl.h"
#include "mod_file.h"

#define HDF_CHECK(x) if ((x) < 0) { mod_hdf_handle_error(); *ierr = 1; return; }
#define HDF_CHECK_NORET(x) if ((x) < 0) { mod_hdf_handle_error(); *ierr = 1; }

#ifdef __cplusplus
extern "C" {
#endif

/** Convert the topmost HDF5 error into a Modeller error */
void mod_hdf_handle_error(void);

/** Create an HDF5 file, with Modeller-style variable expansion */
hid_t mod_hdf_create(const char *filename, unsigned flags, hid_t create_plist,
                     hid_t access_plist, struct mod_file *file_info);

/** Open an HDF5 file, with Modeller-style compression and variable expansion */
hid_t mod_hdf_open(const char *filename, unsigned flags, hid_t access_plist,
                   struct mod_file *file_info);

/** Close an HDF5 file, and clean up after any decompression */
herr_t mod_hdf_close(hid_t file_id, struct mod_file *file_info);

/** Get the size of an HDF5 dataset - it must be nD */
herr_t mod_dataset_get_ndsize(hid_t loc_id, const char *dset_name, int n,
                              hsize_t *dims);

/** Read an HDF5 dataset, but make sure the dimensions are correct first */
herr_t mod_dataset_read(hid_t loc_id, const char *dset_name, int rank,
                        const hsize_t *dims, hid_t type_id, void *buffer);

/** Read an int HDF5 dataset, but make sure the dimensions are correct first */
herr_t mod_dataset_read_int(hid_t loc_id, const char *dset_name, int rank,
                            const hsize_t *dims, int *buffer);

/** Read an int64 HDF5 dataset, but make sure the dimensions are correct
    first */
herr_t mod_dataset_read_int64(hid_t loc_id, const char *dset_name, int rank,
                              const hsize_t *dims, gint64 *buffer);

/** Read a float HDF5 dataset, but make sure the dimensions are correct first */
herr_t mod_dataset_read_float(hid_t loc_id, const char *dset_name, int rank,
                              const hsize_t *dims, float *buffer);

/** Read a double HDF5 dataset, but make sure the dimensions are correct
    first */
herr_t mod_dataset_read_double(hid_t loc_id, const char *dset_name, int rank,
                               const hsize_t *dims, double *buffer);

/** Read a char HDF5 dataset, but make sure the dimensions are correct first */
herr_t mod_dataset_read_char(hid_t loc_id, const char *dset_name, int rank,
                             const hsize_t *dims, char *buffer);

/** Append to an extensible 1D HDF5 dataset */
herr_t mod_dataset_append(hid_t dataset, hsize_t *dataset_size, hsize_t bufsize,
                          hid_t type_id, const void *buffer);

/** Append to an extensible uchar 1D HDF5 dataset */
herr_t mod_dataset_append_uchar(hid_t dataset, hsize_t *dataset_size,
                                hsize_t bufsize, const unsigned char *buffer);

/** Write part of a sequence database in HDF5 format */
herr_t mod_hdf_write_seqdb(hid_t file_id, int nchn, gint64 iseqdb[],
                           int nseqdb[], char codep[], float resoldb[],
                           char prottypdb[], int codep_len);

#ifdef __cplusplus
}
#endif
#endif  /* MOD_HDF5_H */
