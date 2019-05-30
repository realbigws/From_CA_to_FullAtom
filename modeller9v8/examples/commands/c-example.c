#include <glib.h>
#include <stdio.h>
#include <stdlib.h>

#include <modeller.h>

/* Example of using Modeller from a C program. This simply reads in a PDB
 * file, prints out some data from that file, and then writes out a new
 * file in MMCIF format.
 *
 * To compile, use (where XXX is your Modeller version):
 * gcc -Wall -o c-example c-example.c `modXXX --cflags --libs` \
 *     `pkg-config --cflags --libs glib-2.0`
 * (If you use a compiler other than gcc, or a non-Unix system, you may need
 * to run 'modXXX --cflags --libs' manually and construct suitable compiler
 * options by hand.)
 *
 * To run, you must ensure that the Modeller dynamic libraries are in your
 * search path. This can be done on most systems by adding the directory
 * reported by 'modXXX --libs' to the LD_LIBRARY_PATH environment variable.
 * (On Mac, set DYLD_LIBRARY_PATH instead. On Windows, PATH. On AIX, LIBPATH.)
 *
 * You must also ensure that Modeller knows where it was installed,
 * and what the license key is. You can either do this by setting the
 * MODINSTALLXXX and KEY_MODELLERXXX environment variables accordingly, or
 * by calling the mod_install_dir_set() and mod_license_key_set() functions
 * before you call mod_start(). For example, if Modeller is installed in
 * /lib/modeller on a 32-bit Linux system, the following would work from the
 * command line (all on one line), where KEY is your license key:
 *     KEY_MODELLERXXX=KEY MODINSTALLXXX=/lib/modeller/
 *     LD_LIBRARY_PATH=/lib/modeller/lib/i386-intel8 ./c-example
 */


/* Exit, reporting the Modeller error, iff one occurred. */
void handle_error(int ierr)
{
  if (ierr != 0) {
    GError *err = mod_error_get();
    fprintf(stderr, "Modeller error: %s\n", err->message);
    g_error_free(err);
    exit(1);
  }
}

int main(void)
{
  struct mod_libraries *libs;
  struct mod_model *mdl;
  struct mod_io_data *io;
  struct mod_file *fh;
  int ierr, *sel1, nsel1;

  /* Uncomment these lines to hard code install location and license key,
     rather than setting MODINSTALLXXX and KEY_MODELLERXXX environment
     variables (see above) */
  /* mod_install_dir_set("/lib/modeller"); */
  /* mod_license_key_set("KEY"); */

  mod_start(&ierr);
  handle_error(ierr);
  mod_header_write();

  mod_log_set(2, 1);
  libs = mod_libraries_new(NULL);
  mod_libraries_read_libs(libs, "${LIB}/restyp.lib", &ierr);
  handle_error(ierr);
  mod_libraries_rand_seed_set(libs, -8123);

  mdl = mod_model_new(NULL);
  io = mod_io_data_new();
  mod_model_read(mdl, io, libs, "../atom_files/2nbt.pdb", "PDB",
                 "FIRST:@LAST:  ", 7, &ierr);
  handle_error(ierr);
  printf("Model of %s solved at resolution %f, rfactor %f\n", mdl->seq.name,
         mdl->seq.resol, mdl->seq.rfactr);
  fh = mod_file_open("new.cif", "w");
  if (fh) {
    mod_selection_all(mdl, &sel1, &nsel1);
    mod_model_write(mdl, libs, sel1, nsel1, fh, "MMCIF", 0, 1, &ierr);
    g_free(sel1);
    mod_file_close(fh, &ierr);
  }
  handle_error(ierr);
  mod_libraries_free(libs);
  mod_model_free(mdl);
  mod_io_data_free(io);

  mod_end();
  return 0;
}
