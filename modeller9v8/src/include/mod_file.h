/** \file mod_file.h       File-handling routines.
 *
 *             Part of MODELLER, Copyright(c) 1989-2010 Andrej Sali
 */

#ifndef MOD_FILE_H
#define MOD_FILE_H

#include <stdio.h>
#include <glib.h>

#ifdef __cplusplus
extern "C" {
#endif

/** Prototype for a function to free callback data */
typedef void (*cb_free)(void *data);

/** Callback function to read data from a file */
typedef size_t (*cb_readfunc)(void *data, char *buffer, size_t bufsize,
                              int *ierr);

/** Callback function to write data to a file */
typedef void (*cb_writefunc)(void *data, const char *buffer, size_t bufsize,
                             int *ierr);

/** A structure to keep track of a file */
struct mod_file {
  /** The FILE pointer used by the low-level C routines */
  FILE *filept;
  /** Callback function to read data from the file */
  cb_readfunc readfunc;
  /** Callback function to write data to the file */
  cb_writefunc writefunc;
  /** Callback function to free data */
  cb_free freefunc;
  /** Data to pass to callback functions */
  void *data;
  /** The fpos_t pointer used by the low-level C routines.
      Should be treated as opaque, since if Large File Support is active,
      this is actually an fpos64_t pointer. */
  fpos_t *fpos;
  /** The name of the file, before variable expansion and/or uncompression */
  char *name;
  /** The real name of the file */
  char *realname;
  /** Access mode ('r', 'w', 'rb', 'wb') */
  char *action;
  /** Error indicator (0=no error, otherwise system errno) */
  int error_num;
  /** The line number within the file */
  int linenum;
  /** The number of line-truncated warnings printed so far */
  int truncwarns;
  /** Set TRUE after readlinef iff the line didn't end in EOL */
  gboolean noend;
  /** Set TRUE after readlinef iff some of the line was discarded */
  gboolean discarded;
  /** TRUE iff the file is compressed */
  gboolean compressed;
};

/** Open a file, and uncompress it if necessary. */
struct mod_file *mod_file_open(const char *path, const char *mode);

/** Create a new mod_file using the given callbacks for stream operations. */
struct mod_file *mod_file_open_stream(cb_readfunc readfunc,
                                      cb_writefunc writefunc,
                                      cb_free freefunc, void *data);

/** Close an open file, and do any other necessary tidy-up if it was
    compressed. The initial value of ierr_inout is used (if an error was already
    set, it is not modified, but emergency cleanup is done here). */
void mod_file_close(struct mod_file *fh, int *ierr_inout);

/** Set the error indicator if any write errors were encountered on the file. */
void mod_file_check_errors(struct mod_file *fh, int *ierr);

/** Write a binary buffer to a file. Return TRUE on success. */
gboolean mod_file_write_buffer(struct mod_file *fh, const char *buffer,
                               size_t bufsize, int *ierr);

/** Read a binary buffer from a file. Return TRUE on success. */
gboolean mod_file_read_buffer(struct mod_file *fh, char *buffer,
                              size_t bufsize, int *ierr);

/** Read a single line from a file */
void mod_file_read_line(FILE *fp, GString *str, int *ierr, int *eof);

/** Read all data from a file */
void mod_file_read_contents(FILE *fp, char **text, unsigned *filelen,
                            int *ierr);

/** Expand out all environment variables in a filename. Additionally, the
    Modeller-specific ${LIB} and ${JOB} variables are expanded. The
    expanded string is returned, and should be freed when no longer needed. */
char *mod_file_expand(const char *filename);

/** Delete a file. */
void mod_file_delete(const char *file, int *ierr);

/** Return TRUE iff the given file exists. */
gboolean mod_file_exists(const char *file);

/** Deprecated interface to mod_file_exists() for TOP interpreter. */
void mod_inquire(const char *file, int *file_exists);

/** Uncompress the file, if already compressed */
char *mod_file_prep(const char *infil, const char *action,
                    gboolean *compressed, int *ierr);

/** Puts the file in the original state. */
void mod_file_unprep(const char *fname, const char *realname,
                     const char *action, gboolean compressed, int *ierr);

/** Return TRUE iff the given file, or a compressed version, exists */
gboolean mod_file_exists_z(const char *filename);

/** Look for the given file, trying the directories and extensions given.
    The full path of the file is returned, if found, otherwise NULL. */
char *mod_file_exists_z_path(const char *filename, char **dirs, char **exts);

/** Deprecated interface to mod_file_exists_z_path for TOP interpreter */
char *mod_fullfn(const char *filename, char **dirs, char **exts, int *ierr);

#ifdef __cplusplus
}
#endif
#endif  /* MOD_FILE_H */
