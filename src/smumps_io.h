/*

   THIS FILE IS PART OF MUMPS VERSION 4.6
   This Version was built on Tue Jan 24 09:35:01 2006


  This version of MUMPS is provided to you free of charge. It is public
  domain, based on public domain software developed during the Esprit IV
  European project PARASOL (1996-1999) by CERFACS, ENSEEIHT-IRIT and RAL. 
  Since this first public domain version in 1999, the developments are
  supported by the following institutions: CERFACS, ENSEEIHT-IRIT, and
  INRIA Rhone-Alpes.

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
/*    $Id: smumps_io.h,v 1.12 2005/12/14 10:23:36 jylexcel Exp $ */

#if defined(Add_)
#define smumps_is_there_finished_request smumps_is_there_finished_request_
#define smumps_clean_request smumps_clean_request_
#define smumps_test_request smumps_test_request_
#define smumps_wait_request smumps_wait_request_
#define smumps_wait_all_requests smumps_wait_all_requests_
#define smumps_low_level_init_ooc_c smumps_low_level_init_ooc_c_
#define smumps_low_level_write_ooc_c smumps_low_level_write_ooc_c_
#define smumps_low_level_read_ooc_c smumps_low_level_read_ooc_c_
#define smumps_low_level_direct_read smumps_low_level_direct_read_
#define smumps_clean_io_data_c smumps_clean_io_data_c_
#define smumps_get_max_nb_req smumps_get_max_nb_req_
#define smumps_get_max_file_size smumps_get_max_file_size_ 
#define smumps_ooc_get_nb_files smumps_ooc_get_nb_files_
#define smumps_ooc_get_file_name smumps_ooc_get_file_name_
#define smumps_ooc_set_file_name smumps_ooc_set_file_name_
#define smumps_ooc_init_vars smumps_ooc_init_vars_
#define smumps_ooc_start_low_level smumps_ooc_start_low_level_
#define smumps_ooc_alloc_pointers smumps_ooc_alloc_pointers_
#define smumps_ooc_print_stats smumps_ooc_print_stats_
#define smumps_ooc_remove_file smumps_ooc_remove_file_
#define smumps_ooc_init_vars smumps_ooc_init_vars_
#elif defined(Add__)
#define smumps_is_there_finished_request smumps_is_there_finished_request__
#define smumps_clean_request smumps_clean_request__
#define smumps_test_request smumps_test_request__
#define smumps_wait_request smumps_wait_request__
#define smumps_wait_all_requests smumps_wait_all_requests__
#define smumps_low_level_init_ooc_c smumps_low_level_init_ooc_c__
#define smumps_low_level_write_ooc_c smumps_low_level_write_ooc_c__
#define smumps_low_level_read_ooc_c smumps_low_level_read_ooc_c__
#define smumps_low_level_direct_read smumps_low_level_direct_read__
#define smumps_clean_io_data_c smumps_clean_io_data_c__
#define smumps_get_max_nb_req smumps_get_max_nb_req__
#define smumps_get_max_file_size smumps_get_max_file_size__ 
#define smumps_ooc_get_nb_files smumps_ooc_get_nb_files__
#define smumps_ooc_get_file_name smumps_ooc_get_file_name__
#define smumps_ooc_set_file_name smumps_ooc_set_file_name__
#define smumps_ooc_init_vars smumps_ooc_init_vars__
#define smumps_ooc_start_low_level smumps_ooc_start_low_level__
#define smumps_ooc_alloc_pointers smumps_ooc_alloc_pointers__
#define smumps_ooc_print_stats smumps_ooc_print_stats__
#define smumps_ooc_remove_file smumps_ooc_remove_file__
#define smumps_ooc_init_vars smumps_ooc_init_vars__
#elif defined(UPPER)
#define smumps_is_there_finished_request SMUMPS_IS_THERE_FINISHED_REQUEST
#define smumps_clean_request SMUMPS_CLEAN_REQUEST
#define smumps_test_request SMUMPS_TEST_REQUEST
#define smumps_wait_request SMUMPS_WAIT_REQUEST
#define smumps_wait_all_requests SMUMPS_WAIT_ALL_REQUESTS
#define smumps_low_level_init_ooc_c SMUMPS_LOW_LEVEL_INIT_OOC_C
#define smumps_low_level_write_ooc_c SMUMPS_LOW_LEVEL_WRITE_OOC_C
#define smumps_low_level_read_ooc_c SMUMPS_LOW_LEVEL_READ_OOC_C
#define smumps_low_level_direct_read SMUMPS_LOW_LEVEL_DIRECT_READ
#define smumps_clean_io_data_c SMUMPS_CLEAN_IO_DATA_C
#define smumps_get_max_nb_req SMUMPS_GET_MAX_NB_REQ
#define smumps_get_max_file_size SMUMPS_GET_MAX_FILE_SIZE
#define smumps_ooc_get_nb_files SMUMPS_OOC_GET_NB_FILES
#define smumps_ooc_get_file_name SMUMPS_OOC_GET_FILE_NAME
#define smumps_ooc_set_file_name SMUMPS_OOC_SET_FILE_NAME
#define smumps_ooc_init_vars SMUMPS_OOC_INIT_VARS
#define smumps_ooc_start_low_level SMUMPS_OOC_START_LOW_LEVEL
#define smumps_ooc_alloc_pointers SMUMPS_OOC_ALLOC_POINTERS
#define smumps_ooc_print_stats SMUMPS_OOC_PRINT_STATS
#define smumps_ooc_remove_file SMUMPS_OOC_REMOVE_FILE
#define smumps_ooc_init_vars SMUMPS_OOC_INIT_VARS
#endif

int smumps_is_there_finished_request(int* flag,int* ierr);

int smumps_clean_request(int* request_id,int* ierr);

int smumps_test_request(int* request_id,int* flag,int* ierr);

int smumps_wait_request(int* request_id,int* ierr);

int smumps_wait_all_requests();

int smumps_low_level_init_ooc_c(int* _myid, int* total_size_io,int* size_element,
			       int* async,char* smumps_dir, char* smumps_file,
			       int* smumps_dim_dir, int* smumps_dim_file,
			       int* ierr);

int smumps_low_level_write_ooc_c( const int * strat_IO, 
				 void * address_block,
				 int * block_size,
				 int * pos_in_file,
				 int * file_number,
				 int * inode,
				 int * request_arg,
				 int * ierr);

int smumps_low_level_read_ooc_c( const int * strat_IO, 
			        void * address_block,
				int * block_size,
				int * from_where,
				int * file_number,
				int * inode,
				int * request_arg,
				int * ierr);

int smumps_low_level_direct_read(void * address_block,
				   int * block_size,
				   int * from_where,
				   int * file_number,
				   int * ierr);

int smumps_clean_io_data_c(int* myid,int* ierr);

int smumps_get_max_nb_req(int *max,int* ierr);

int smumps_get_max_file_size(double * max_ooc_file_size);

int smumps_ooc_get_nb_files(int* nb_files);

int smumps_ooc_get_file_name(int* indice,char* name,int* length);

int smumps_ooc_set_file_name(int* indice,char* name,int* length,int* ierr);

int smumps_ooc_alloc_pointers(int* dim,int* ierr);

int smumps_ooc_init_vars(int* myid_arg, int* nb_file_arg,
			int* size_element,int* async,
			char* smumps_dir, char* smumps_file,
			int* smumps_dim_dir, int* smumps_dim_file,
			int *ierr);

int smumps_ooc_start_low_level(int* ierr);

int smumps_ooc_print_stats();

int smumps_ooc_remove_file(char *name,int* ierr);

