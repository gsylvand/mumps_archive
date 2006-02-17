/*

   THIS FILE IS PART OF MUMPS VERSION 4.6.1
   This Version was built on Fri Feb 17 14:27:51 2006


  This version of MUMPS is provided to you free of charge. It is public
  domain, based on public domain software developed during the Esprit IV
  European project PARASOL (1996-1999) by CERFACS, ENSEEIHT-IRIT and RAL. 
  Since this first public domain version in 1999, the developments are
  supported by the following institutions: CERFACS, ENSEEIHT-IRIT, and
  INRIA.

  Main contributors are Patrick Amestoy, Iain Duff, Abdou Guermouche,
  Jacko Koster, Jean-Yves L'Excellent, and Stephane Pralet.

  Up-to-date copies of the MUMPS package can be obtained
  from the Web pages http://www.enseeiht.fr/apo/MUMPS/
  or http://graal.ens-lyon.fr/MUMPS


   THIS MATERIAL IS PROVIDED AS IS, WITH ABSOLUTELY NO WARRANTY
   EXPRESSED OR IMPLIED.  ANY USE IS AT YOUR OWN RISK.


  User documentation of any code that uses this software can
  include this complete notice. You can acknowledge (using
  references [1], [2], and [3] the contribution of this package
  in any scientific publication dependent upon the use of the
  package. You shall use reasonable endeavours to notify
  the authors of the package of this publication.

   [1] P. R. Amestoy, I. S. Duff and  J.-Y. L'Excellent (1998),
   Multifrontal parallel distributed symmetric and unsymmetric solvers,
   in Comput. Methods in Appl. Mech. Eng., 184,  501-520 (2000).

   [2] P. R. Amestoy, I. S. Duff, J. Koster and  J.-Y. L'Excellent,
   A fully asynchronous multifrontal solver using distributed dynamic
   scheduling, SIAM Journal of Matrix Analysis and Applications,
   Vol 23, No 1, pp 15-41 (2001).

   [3] P. R. Amestoy and A. Guermouche and J.-Y. L'Excellent and
   S. Pralet (2005), Hybrid scheduling for the parallel solution
   of linear systems. Accepted to Parallel Computing.

*/
/*    $Id: smumps_io_basic.h,v 1.28 2006/01/19 17:15:23 aguermou Exp $  */

#define MAX_FILE_SIZE 1879048192 /* (2^31)-1-(2^27) */
/*                                                      */
/* Important Note :                                     */
/* ================                                     */
/* On GNU systems, __USE_GNU must be defined to have    */
/* access to the O_DIRECT I/O flag.                     */
/*                                                      */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#ifndef _WIN32
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <sys/time.h>
#include <time.h>

#endif

#ifndef _WIN32
#if defined(SP_) || defined(__USE_GNU)
#define SMUMPS_IO_FLAGS O_RDWR|O_DIRECT
#else
#define SMUMPS_IO_FLAGS O_RDWR
#endif
#endif

#define IO_SYNC      0

#define IO_ASYNC_TH  1

#define IO_ASYNC_AIO 2

#define IO_READ 1
#define IO_WRITE 0

#define UNITIALIZED "NAME_NOT_INITIALIZED"

#if defined (CINES_)
#define SMUMPS_OOC_DEFAULT_DIR "/tmpp/eagullo"
#else
#define SMUMPS_OOC_DEFAULT_DIR "/tmp"
#endif

#define SEPARATOR "/"

int _smumps_next_file();
void _smumps_update_current_file_position();
int _smumps_compute_where_to_write(const double to_be_written);
int _smumps_prepare_pointers_for_write(double to_be_written,int * pos_in_file, int * file_number);
int smumps_io_do_write_block(void * address_block,int * block_size,int * pos_in_file,int * file_number,int * ierr);
int smumps_io_do_read_block(void * address_block,int * block_size,int * from_where,int * file_number,int * ierr);
int _smumps_compute_nb_concerned_files(int * block_size, int * nb_concerned_files);
int smumps_free_file_pointers();
int smumps_init_file_structure(int* _myid, int* total_size_io,int* size_element);
int smumps_init_file_name(char* smumps_dir,char* smumps_file,int* smumps_dim_dir,int* smumps_dim_file,int* _myid);
int smumps_io_alloc_file_struct(int* nb);
int smumps_io_get_nb_files(int* nb_files);
int smumps_io_get_file_name(int* indice,char* name,int* length);
int smumps_io_alloc_pointers(int * dim);
int smumps_io_init_vars(int* myid_arg, int* nb_file_arg,int* size_element,int* async_arg);
int smumps_io_set_file_name(int* indice,char* name,int* length);
int smumps_io_open_files_for_read();
int smumps_io_set_last_file(int* dim);

#ifndef _WIN32  
#ifndef WITHOUT_PFUNC

int smumps_io_protect_pointers();
int smumps_io_unprotect_pointers();
int smumps_io_init_pointers_lock();
int smumps_io_destroy_pointers_lock();

#endif /*WITHOUT_PFUNC*/
#endif /*_WIN32*/


