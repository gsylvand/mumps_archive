/*

   THIS FILE IS PART OF MUMPS VERSION 4.3.2
   This Version was built on Wed Nov 12 16:57:09 2003

 COPYRIGHT (c) 1996-2003 P. R. Amestoy, I. S. Duff, J. Koster,
                       J.-Y. L'Excellent

  CERFACS      , Toulouse    (France)  (http://www.cerfacs.fr)
  ENSEEIHT-IRIT, Toulouse    (France)  (http://www.enseeiht.fr)
  INRIA                      (France)  (http://www.inria.fr)
  PARALLAB     , Bergen      (Norway)  (http://www.parallab.uib.no)

 All rights reserved.

  Your use or distribution of the package implies that you agree
  with this License. Up-to-date copies of the MUMPS package can be
  obtained from the Web page http://www.enseeiht.fr/apo/MUMPS/

  This package is provided to you free of charge. It was
  initially based on public domain software developed during
  the European Esprit IV project PARASOL (1996-1999).

  THIS MATERIAL IS PROVIDED AS IS, WITH ABSOLUTELY NO WARRANTY
  EXPRESSED OR IMPLIED.  ANY USE IS AT YOUR OWN RISK.

  Permission is hereby granted to use or copy this
  package provided that the Copyright and this License is
  retained on all copies and that the package is used
  under the same terms and conditions. User documentation
  of any code that uses this software should include this
  complete Copyright notice and this License.

  You can modify this code but, at no time shall the right
  or title to all or any part of this package pass to you.
  All information relating to any alteration or addition
  made to this package for the purposes of extending the
  capabilities or enhancing the performance of this package
  shall be made available free of charge to the authors for
  any purpose.

  You shall acknowledge (using references [1] and [2])
  the contribution of this package in any publication
  of material dependent upon the use of the package.
  You shall use reasonable endeavours to notify
  the authors of the package of this publication.


  [1] P. R. Amestoy, I. S. Duff and  J.-Y. L'Excellent (1998),
  Multifrontal parallel distributed symmetric and unsymmetric solvers,
  in Comput. Methods in Appl. Mech. Eng., 184,  501-520 (2000).
  An early version appeared as a Technical Report ENSEEIHT-IRIT (1998)
  and is available at http://www.enseeiht.fr/apo/MUMPS/.

  [2] P. R. Amestoy, I. S. Duff, J. Koster and  J.-Y. L'Excellent,
  A fully asynchronous multifrontal solver using distributed dynamic
  scheduling, SIAM Journal of Matrix Analysis and Applications,
  Vol 23, No 1, pp 15-41 (2001).
  An  early version appeared as a Technical Report ENSEEIHT-IRIT,
  RT/APO/99/2 (1999)  and is available at http://www.enseeiht.fr/apo/MUMPS/.

  None of the text from the Copyright notice up to and
  including this line shall be removed or altered in any way.

*/
/* $Id: dmumps_orderings.c,v 1.17 2003/10/01 14:11:01 jylexcel Exp $ */
/*
 * This file contains interfaces to external ordering packages.
 * At the moment, PORD (J. Schulze) and SCOTCH are interfaced.
 */

#include <dmumps_orderings.h>

#if defined(pord)

#include <space.h>

/************************************************************
dmumps_pord is used in usedbyanalysis.F
	permutation and inverse permutation not set in output,
	but may be printed in default file: "perm_pord" and "iperm_pord",
	if associated part uncommneted.
	But, if uncommetnted a bug occurs in psl_ma41_analysi.F
******************************************************************/	

/*********************************************************/

int dmumps_pord
( 
	int 	nvtx, 
	int	nedges,
	int	*xadj_pe,
	int	*adjncy,
	int	*nv
)
{
/**********************************
Argument Comments:

input:
-----
- nvtx 		: dimension of the Problem (N)
- nedges	: number of entries (NZ)
- adjncy	: non-zeros entries (IW input)

input/output:
-------------
- xadj_pe		: pointer through beginning of column non-zeros entries (PTRAR)
- on exit, "father array" (PE)

ouput:
------
- nv		: "son array" (NV)
*************************************/

	graph_t    	*G;
 	elimtree_t 	*T;
	timings_t  cpus[12];
	options_t  options[] = { SPACE_ORDTYPE, SPACE_NODE_SELECTION1,
                    SPACE_NODE_SELECTION2, SPACE_NODE_SELECTION3,
                    SPACE_DOMAIN_SIZE, 0 };
	int *ncolfactor, *ncolupdate, *parent, *vtx2front;
  	int *first, *link, nfronts, J, K, u, vertex, vertex_root, count;
      
      /**************************************************
       declaration to uncomment if printing ordering 
      ***************************************************

         FILE *fp1, *fp2;
         int  *perm,  *iperm;  
      */

/*** decalage des indices couteux dans un premier temps:
****  A modifier dans une version ulterieure de MA41GD  */
  
	for (u = nvtx; u >= 0; u--)
	{
		xadj_pe[u] = xadj_pe[u] - 1;
	}
	for (K = nedges-1; K >= 0; K--)
	{
		adjncy[K] = adjncy[K] - 1;
	}
	
	
	
 /* initialization of the graph */	
	
	mymalloc(G, 1, graph_t);
	G->xadj			= xadj_pe;
	G->adjncy		= adjncy;
	
	
	mymalloc(G->vwght, nvtx, int);
	G->nvtx = nvtx;
  	G->nedges = nedges;
  	G->type = UNWEIGHTED;
  	G->totvwght = nvtx;
  	for (u = 0; u < nvtx; u++)
    		G->vwght[u] = 1;
			

  /* main function of the Ordering */
	T = SPACE_ordering(G, options, cpus);

	
  	nfronts = T->nfronts;
  	ncolfactor = T->ncolfactor;
  	ncolupdate = T->ncolupdate;
  	parent = T->parent;
  /*    firstchild = T->firstchild; */
	vtx2front = T->vtx2front;
	
	
	
    /* -----------------------------------------------------------
     store the vertices/columns of a front in a bucket structure
     ----------------------------------------------------------- */
  	mymalloc(first, nfronts, int);
  	mymalloc(link, nvtx, int);


  	for (J = 0; J < nfronts; J++)
    	  first[J] = -1;
  	
	for (u = nvtx-1; u >= 0; u--)
   	{ 
	  J = vtx2front[u];
	  link[u] = first[J];
	  first[J] = u;
   	}

  /* -----------------------------------------------------------
     fill the two arrays corresponding to the MUMPS tree structure
     ----------------------------------------------------------- */
  count = 0;
   
  for (K = firstPostorder(T); K != -1; K = nextPostorder(T, K))
     {
  	vertex_root = first[K];
	if (vertex_root == -1)
          {
            /* JY: I think this cannot happen */
            printf(" Internal error in dmumps_pord (cf JY), %d\n",K);
            exit(-1);
          }

	/* for the principal column of the supervariable */
	if (parent[K] == -1)
          xadj_pe[vertex_root] = 0; /* root of the tree */
	else
	  xadj_pe[vertex_root] = - (first[parent[K]]+1);
	nv[vertex_root] = ncolfactor[K] + ncolupdate[K];
	count++;
	for (vertex = link[vertex_root]; vertex != -1; vertex = link[vertex])
        /* for the secondary columns of the supervariable */
	{ 
        	xadj_pe[vertex] = - (vertex_root+1);
		nv[vertex] = 0;
		count++;
        }
  
  }
  

/*******************************************
begining of printing ordering
********************************************

************************************************************ 
This part compile but there is an exucution mistake later in
psl_ma41_ordering when allocate id%POIDS.
The bug accurs when the following for loop is uncommented:
when perm and iperm are filled. 
************************************************************

   mymalloc(perm, nvtx, int); 
   mymalloc(iperm, nvtx, int);  

   for (K = firstPostorder(T); K != -1; K = nextPostorder(T, K))
   {
  	vertex_root = first[K];
	if (vertex_root != -1)
	{	
		perm[vertex_root] = count+1;
		iperm[count] = vertex_root+1;
		count++;
	}
	for (vertex = link[vertex_root]; vertex != -1; vertex = link[vertex])
	{ 
		perm[vertex] = count+1;
		iperm[count] = vertex+1;
		count++;
        }
    }
  
   if ((fp1 = fopen("perm_pord", "w")) == NULL)
   { fprintf(stderr, "\nError in function greg_pord\n"
          "  unable to open file perm_pord for writing\n");
     quit();
   }
   if ((fp2 = fopen("iperm_pord", "w")) == NULL)
   { fprintf(stderr, "\nError in function greg_pord\n"
          "  unable to open file iperm_pord for writing\n");
     quit();
   }
  
   for (K = 0; K < nvtx; K++)
   {
    fprintf(fp1, "%d\n", perm[K]);
   }
   for (K = 0; K < nvtx; K++)
   {
    fprintf(fp2, "%d\n", iperm[K]);
   }
  
   fclose(fp1);
   fclose(fp2);
  
   free(perm); 
   free(iperm);

*************************
end of printing ordering 
**************************/

  /* ----------------------
     free memory and return
     ---------------------- */
  free(first); free(link); 
  free(G->vwght);
  free(G);
  freeElimTree(T);
	
	return (0);
}

FORTRAN (
 DMUMPS_PORDF, dmumps_pordf, (
int 	*nvtx, 		\
int	*nedges,	\
int	*xadj,		\
int	*adjncy,	\
int	*nv,		\
int	*ncmpa
),
(nvtx, nedges, xadj, adjncy, nv, ncmpa))
{
*ncmpa = dmumps_pord (*nvtx, *nedges, xadj, adjncy, nv);
}
#endif



#if defined(scotch)
/* Fortran interface to SCOTCH */
FORTRAN (                                               \
        DMUMPS_395, dmumps_scotch, (                   \
        const int * const           n,                  \
        const int * const           iwlen,              \
        int * const                 petab,              \
        const int * const           pfree,              \
        int * const                 lentab,             \
        int * const                 iwtab,              \
        int * const                 nvtab,              \
        int * const                 elentab,            \
        int * const                 lasttab,            \
        int * const                 ncmpa),             \
	(n, iwlen, petab, pfree, lentab, iwtab, nvtab, elentab, lasttab, ncmpa))
{
	  *ncmpa = esmumps (*n, *iwlen, petab, *pfree, lentab, iwtab, nvtab, elentab, lasttab);
}
#endif

