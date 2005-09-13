/*

   THIS FILE IS PART OF MUMPS VERSION 4.5.3
   This Version was built on Tue Sep 13 14:45:37 2005

  This version of MUMPS is provided to you free of charge. It is public
  domain, based on public domain software developed during the Esprit IV
  European project PARASOL (1996-1999). 
  It has also been supported by CERFACS, ENSEEIHT-IRIT-CNRS, 
  LIP-INRIA-Rhone-Alpes, and PARALLAB.

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
/* $Id: dmumps_distribution.c,v 1.9 2003/04/25 12:25:46 jylexcel Exp $ */
/******************************************************************
 *   distribution.c    
 *   Tools to compute the distribution of the tree on the processors
 *   (see comments of the routine "dmumps_distri"
 *
 *   Created by Vincent Espirat, ENSEEIHT.
 *
 ******************************************************************
 *
 *  Modified by Patrick Amestoy, Nov. 1999.
 *
 *  Modified by Jean-Yves L'Excellent, April 97
 *
 *  - Avoid copy of na,ne,nd,frere,fils into local arrays
 *  - Use dynamic allocation instead of local arrays of size NMAX
 *  - Remove some unnecessary memory allocations, that were
 *    never freed
 *  - Suppress inclusion of the file cpara.h
 *  - Suppress the array tid since we are now using MPI
 *  - Changed symbol names to avoid conflicts
 *  - Various bug fixes
 *  - list is now freed
 *
 ******************************************************************/

#include <stdio.h>
#include <stdlib.h>
#define FALSE 0
#define TRUE 1
#define MINFLOPS 50000000
#define MONint int

/******************************************/
/* calcul du cout en nombre d'operations   */  
/******************************************/

/* calcul de cout de chaque sous-arbre */
    
void dmumps_frootcout(pos,ne,nd,frere,fils,tcout,ncout,k50)
   MONint pos;
   MONint k50;
   MONint *ne,*nd,*frere,*fils;
    double *tcout,*ncout;
{
  double cout;
 MONint posi;
 MONint nb,npiv,nfront;
  
  npiv = 1;
  posi = fils[pos];
  while(posi>0) 
     {
       posi  = fils[posi];
       npiv  = npiv+1;
     }
  nfront = nd[pos];
  if ( (nfront==npiv) && (npiv==1) ) 
   {
   cout = 0.0;
   }
  else
   {
   if ( k50 == 0 )
    {
     cout = 2.0 * ((double)nfront) * ((double)npiv) * ((double)(nfront-npiv-1)) +
      ((double)npiv) * (npiv+1.0) *(2.0*npiv+1)/3.0  + (2.0*nfront-npiv-1)* ((double)npiv) / 2.0;
    }
   else
    {
      cout = ((double)npiv) * (   ((double)nfront) * ((double)nfront) + 2.0*nfront - (nfront+1.0)*(npiv+1.0)
           +  (npiv+1.0) * (2.0*npiv+1)/6.0 );
    }
   }
  if (cout<0)
    {
     printf("Warning in distribution : cout = %lf, npiv=%d, nfront=%d",cout,npiv,nfront);
    }
  ncout[pos] = cout ;

  if (ne[pos] != 0 ) 
    {
     posi = fils[pos];
     while(posi>0) posi = fils[posi]; 
     posi = - posi;
     for ( nb = 0 ; nb < ne[pos] ; nb++) 
        { 
         dmumps_frootcout(posi,ne,nd,frere,fils,tcout,ncout,k50);
         cout = cout + tcout[posi];
         posi = frere[posi];
        } 
    }
  tcout[pos] = cout; 
 }

/* calcul du cout de l'arbre complet */

 void dmumps_fcout(n,ne,nd,frere,fils,tcout,ncout,k50)
   MONint n,k50;
   MONint *ne,*nd,*frere,*fils;
    double *tcout,*ncout;
{
 MONint pos;

  /* recherche des racines de l'arbre */
  for (pos =  1; pos < n+1 ; pos++)
      {
        if (frere[pos] ==(n+1)) 
          {
           tcout[pos]=0;
           ncout[pos]=0;
          }
        else 
          {
           /* noeud principal */
           if(frere[pos] == 0) dmumps_frootcout(pos,ne,nd,frere,fils,tcout,ncout,k50);
          }
      } 

}

/**********************************************************************/

/*********************************/
/* structure de la liste chainee */
/*********************************/

typedef struct tlist {
     double lcout;
    MONint num;
     struct tlist *suivant;} liste;
/************/
liste *dmumps_creerliste(cout,no)
  double cout;
 MONint no;
{
  liste *list0;
  list0 = (liste *)malloc(sizeof(liste));
  list0->lcout = cout;
  list0->num = no;
  list0->suivant = NULL;
  return list0;
}
/************/
void dmumps_ajouter(cout,no,list)  /* ajoute i a la liste non vide */
   double cout;
  MONint no;             
   liste *list;
{
 while(list->suivant) list = list->suivant;
 list->suivant = dmumps_creerliste(cout,no);
}
/************/
liste *dmumps_ajouter_avant(cout,no,list)  /* ajoute i a la liste non vide */
   double cout;
  MONint  no;             
   liste *list;
{
  liste *list1;
  list1 = dmumps_creerliste(cout,no);
  list1->suivant = list;
  return list1;
}

/**********************************************************************/

/**************************************************/
/* remplacement d'un noeud par une liste de noeud */
/**************************************************/
liste *dmumps_remplacer_aux(cout,list,list1,list2)
  double cout;
  liste *list,*list1,*list2;
{
 if (list->lcout==cout)
     {
      list1->suivant = list->suivant;
      return list2;
     }
 else
     {
      list->suivant = dmumps_remplacer_aux(cout,list->suivant,list1,list2);
      return list;
     }
}
/************/
 liste *dmumps_remplacer(cout,list,list1)
  double cout;
  liste *list;
  liste *list1;
{
  liste *list2;
/*list2 = (liste *)malloc(sizeof(liste));*/

  list2 = list1;
  while(list1->suivant) list1=list1->suivant;

  list = dmumps_remplacer_aux(cout,list,list1,list2);
  return list;
}
/**********************************************************************/

/******************************/
/* classement dans un tableau */
/******************************/
  void dmumps_classer(tab1,tab2,nb)
    double *tab1;
   MONint *tab2;
   MONint nb;
 {
   MONint i,fin;
    double aux1;
   MONint aux2;
 
    fin = FALSE;
    i=1;
    while(!fin)
     {
      if (i==nb) fin = TRUE;
      else
       { 
        if (tab1[i]<tab1[i+1])
         { aux1 = tab1[i+1];
           tab1[i+1] = tab1[i];
           tab1[i] = aux1;
           aux2 = tab2[i+1];
           tab2[i+1] = tab2[i];
           tab2[i] = aux2;
           if (i!=1) i--;
         }
        else i++;
       }
     }
}
/**********************************************************************/

/************/
 MONint dmumps_mintab(tabcout,nb)
  double *tabcout;
  MONint nb;
{
  MONint pos,i;
  double min;

  min = tabcout[1];
  pos = 1;
  for (i=2;i<nb+1;i++)
   {
    if (min>tabcout[i]) 
      {
        min = tabcout[i];
        pos = i;
      }
   }
  return pos;
}

/*******************************************************************/
/* rangement de la strate dans le tableau des processeurs procnode1 */
/*******************************************************************/

/* rangement simple pour la strate initiale pour un equilibrage */
/*  fait de maniere grossiere */

  void dmumps_rangersimp(list,n,nbproc,procnode1,tabcout,poids)
    liste *list;
   MONint nbproc,n;
    double *tabcout;
   MONint *procnode1;
   MONint *poids;

 {
   double *tab1;
  MONint *tab2;  /* tableau de travail cout et no noeud */
  MONint nb,pos;
/*     MONint j,posi,posj,fin */
  MONint i,mini; 

#if defined(check)
   printf(" Warning dmumps_rangersimp: simple strate0 mapping \n");
#endif
   tab1 = ( double * ) malloc( n * sizeof( double ) );
   tab2 = ( MONint * ) malloc( n * sizeof( MONint ) );
   nb=0;
   while(list)
    { 
     nb++;
     tab1[nb] = list->lcout;
     tab2[nb] = list->num;
     list = list->suivant;
    } 
   dmumps_classer(tab1,tab2,nb);
   if (nb<=nbproc) 
    {
     for (pos=1;pos<=nb;pos++)
     { tabcout[pos] = poids[pos]*tab1[pos];
       procnode1[tab2[pos]] = pos;
     }
    }
   else
    {
     for (i=1;i<=nb;i++)
     {
      mini = dmumps_mintab(tabcout,nbproc);
/* Modif done on November 1999:
   We want to avoid having roots that are leaves
   (ie diagonal elements of cost 0) to be all mapped
    on the same proc)  */
      if (tab1[i] == 0.0 ) 
      {
       procnode1[tab2[i]] = tab2[i]%nbproc + 1;
      }
      else
      {
       tabcout[mini] = tabcout[mini] + poids[mini]*tab1[i];
       procnode1[tab2[i]] = mini;
      }
     }
    }

/* autre maniere de ranger simplement en tenant compte du poids des noeuds */
/*  on fait un classement des charges */
/*  puis on dispose en faisant des aller-retour dans le tableau des processeur*/
/*   ex : pour 3 processeurs, 5 noeuds de la strate initiale  ranger */
/*         par ordre decroissant de plus gros cout au plus faible */
/*        le processeur 1 traite le noeud 1 */
/*        le processeur 2 traite les noeud 2 et 5 */
/*        le processeur 3 trait les noeuds 3 et 4 */
/*  on a une charge moyenne sur chaque processeur */
/*      for (pos=1;pos<=nbproc;pos++)*/
/*        { tabcout[pos] = poids[pos]*tab1[pos];*/
/*          procnode1[tab2[pos]] = pos;  */
/*        }*/
/*      if (nb<(2*nbproc+1))*/
/*       {*/
/*         for (j=nbproc;j>(2*nbproc-nb);j--)*/
/*           { tabcout[j] = tabcout[j] + poids[j]*tab1[pos] ;*/
/*             procnode1[tab2[pos++]] = j;*/
/*           }*/
/*       }*/
/*      else*/
/*       {*/
/*         posi = nb - nbproc ; */
/*         while (posi>0)*/
/*          {*/
/*           if (posi>=2*nbproc)*/
/*            {*/
/*             for (j=nbproc;j>0;j--)*/
/*              { tabcout[j] = tabcout[j] + poids[j]*tab1[pos] ;*/
/*                procnode1[tab2[pos++]] = j;*/
/*              } */
/*             for (j=1;j<=nbproc;j++)*/
/*              { tabcout[j] = tabcout[j] + poids[j]*tab1[pos] ;*/
/*                procnode1[tab2[pos++]] = j;*/
/*              }*/
/*             posi = posi - 2*nbproc;*/
/*            }*/
/*           else*/
/*            {*/
/*             if (posi<nbproc) */
/*              {*/
/*               for (j=nbproc;j>(nbproc-posi);j--)*/
/*                { tabcout[j] = tabcout[j] + poids[j]*tab1[pos] ;*/
/*                  procnode1[tab2[pos++]] = j;*/
/*                }*/
/*              }*/
/*             else*/
/*              {*/
/*               for (j=nbproc;j>0;j--)*/
/*                { tabcout[j] = tabcout[j] + poids[j]*tab1[pos] ;*/
/*                  procnode1[tab2[pos++]] = j;*/
/*                }*/
/*               posi = posi - nbproc;*/
/*               for (j=1;j<=posi;j++)*/
/*                { tabcout[j] = tabcout[j] + poids[j]*tab1[pos] ;*/
/*                  procnode1[tab2[pos++]] = j;*/
/*                }*/
/*              }*/
/*            posi = -1;*/
/*            }*/
/*          }*/
/*       }*/

 free(tab1);
 free(tab2);
 }
/************/

/* rangement de la strate initiale en tenant compte du coutmoy*/
/* rangement par groupement des noeuds les plus proches sur un meme processeur*/
  void dmumps_ranger(list,n,nbproc,procnode1,tabcout,poids)
    liste *list;
    MONint nbproc,n;
    double *tabcout;
    MONint *procnode1;
    MONint *poids;

 {
   double *tab1;
   MONint *tab2;
/* Tableaux de travail cout et noeud */
   MONint i,nb,pos,mini;
   double cout,coutmoy;
   liste *list1;

   tab1 = (double *) malloc( n*sizeof( double ) );
   tab2 = ( MONint *) malloc( n*sizeof( MONint ) );
/* list1 = (liste *)malloc(sizeof(liste));*/
   list1 = list;
   cout = 0.0;
   /* printf("dmumps_ranger: liste des couts : \n");*/
   while(list1)
    {
     cout = cout +list1->lcout;
     /*  printf("%lf (node %d) , ",list1->lcout,list1->num); */
     list1 = list1->suivant;
    }

/* coutmoy est le cout moyen de toute la charge a distribuee */
   coutmoy = cout/nbproc;
   /* printf("\n cout total=%lf, coutmoy=%lf\n",cout,coutmoy);*/
   pos = 0;
   list1 = list;
   while(list1&&(pos<nbproc))
    {
/* on range les noeuds sur le meme processeur jusqu'a avoir */
/* une charge de coutmoy, puis on passe a un autre processeur */
      pos++;
      tabcout[pos] = poids[pos]*list1->lcout;
      if ( (list1->lcout) == 0.0 ) 
        procnode1[list1->num] = ((list1->num)%nbproc) + 1;
      else
       procnode1[list1->num] = pos;
      list1 = list1->suivant;
      while ( list1 &&(tabcout[pos]+list1->lcout)<coutmoy) 
       {
/* Modif done on November 1999: 
   We want to avoid having roots that are leaves 
   (ie diagonal elements of cost 0) to be all mapped 
    on the same proc)  */
         if ( (list1->lcout) == 0.0 ) 
         {
          procnode1[list1->num] = ((list1->num)%nbproc) + 1;
         }
         else
         {
         tabcout[pos] = tabcout[pos] + poids[pos]*list1->lcout;
         procnode1[list1->num] = pos;
         }
         list1 = list1->suivant;
       }
    }

/* rangement des noeuds non pris car depassement de coutmoy */
   nb=0;
   while(list1)
    {
     nb++;
     tab1[nb] = list1->lcout;
     tab2[nb] = list1->num;
     list1 = list1->suivant;
    }
   if (nb)
    {
     dmumps_classer(tab1,tab2,nb);
     for (i=1;i<=nb;i++)
     {
/* mini est le processeur de charge minimum */
      mini = dmumps_mintab(tabcout,nbproc);
      tabcout[mini] = tabcout[mini] + poids[mini]*tab1[i];
      if (tab1[i] == 0 )
         procnode1[tab2[i]] = tab2[i]%nbproc + 1;
      else
         procnode1[tab2[i]] = mini;
     }
    }
 free(tab1);
 free(tab2);
 /* printf("nbproc=%d  tabcout=",nbproc); 
 for(i=1;i<=nbproc;i++)
   printf("%lf",tabcout[i]);
 printf("\n"); */
 }
 

/************/
/* rangement d'une strate */
  void dmumps_rangerstrat(list,n,nbproc,procnode1,tabcout,poids)
    liste *list;
    MONint nbproc,n;
    double *tabcout;
    MONint *procnode1;
    MONint *poids;

{
  MONint i,nb,mini;
  double *tab1;
  MONint *tab2;

  tab1 = (double *) malloc( n * sizeof( double ) );
  tab2 = ( MONint *) malloc( n * sizeof( MONint ) );
  nb = 0;
  while(list)
  { 
   nb++;
   tab1[nb] = list->lcout;
   tab2[nb] = list->num;
   list = list->suivant;
  }
  dmumps_classer(tab1,tab2,nb);
  for (i=1;i<=nb;i++)
    {
     mini =dmumps_mintab(tabcout,nbproc);
     tabcout[mini] = tabcout[mini] + poids[mini]*tab1[i];
     procnode1[tab2[i]] = mini;
    }

 free(tab1);
 free(tab2);
}  

/************/
 MONint dmumps_equilibrage(tabcout,nbproc)
    double *tabcout;
    MONint nbproc;

  {
   double maxi,mini;
   MONint i;
   
   /* printf("Entering Equilibrage, nbproc= %d\n", nbproc); */
   for (i=2,maxi=tabcout[1],mini=tabcout[1];i<(nbproc+1);i++)
      {
       if (tabcout[i]>maxi) maxi=tabcout[i];
       if (tabcout[i]<mini) mini=tabcout[i];
       }
    if ( maxi < MINFLOPS )  return TRUE;
    if ( nbproc < 47 )
    {
       if (maxi>(2*mini))
         {
           for(i=1;i<nbproc+1;i++) tabcout[i] = 0; 
           return FALSE;
         }
       else return TRUE; 
     }
     /* We do no want 
        to go too much down in the tree 
        second level of // will distribute work */
     else if ( nbproc < 128)
     {
       /* 4.1.5 : if (maxi>(5*mini)) */
       if (maxi>(4*mini))
         {
           for(i=1;i<nbproc+1;i++) tabcout[i] = 0;
           return FALSE;
         }
       else return TRUE; 
     }
     else if ( nbproc < 256)
     {
       /* 4.1.5 : if (maxi>(10*mini)) */
       if (maxi>(6*mini))
         {
           for(i=1;i<nbproc+1;i++) tabcout[i] = 0;
           return FALSE;
         }
       else return TRUE; 
     }
     else if ( nbproc < 512)
     {
       /* 4.1.5 if (maxi>(15*mini)) */
       if (maxi>(8*mini))
         {
           for(i=1;i<nbproc+1;i++) tabcout[i] = 0;
           return FALSE;
         }
       else return TRUE; 
     }
     else
     {
       /* 4.1.5 if (maxi>(30*mini)) */
       if (maxi>(10*mini))
         {
           for(i=1;i<nbproc+1;i++) tabcout[i] = 0;
           return FALSE;
         }
       else return TRUE; 
     }

  }

/************************************************/
/* rangement dans procnode1 un sous-arbre entier*/
/************************************************/

  void dmumps_sousarbreaux(procnode1,noeud,noproc,frere,fils)
  MONint *procnode1,*frere,*fils;
  MONint noeud,noproc;
 
 {
   procnode1[noeud] = noproc;
   noeud = fils[noeud];
   if (noeud)
   {
     while (noeud>0)
      { procnode1[noeud] = noproc;
        noeud = fils[noeud];
      }
     noeud = -noeud;
     while (noeud >0 )
      {
       dmumps_sousarbreaux(procnode1,noeud,noproc,frere,fils); 
       noeud = frere[noeud];
      }
   }
  }
   
/************/
  

  void dmumps_sousarbre(list,procnode1,frere,fils)
   liste *list;
   MONint *procnode1,*frere,*fils;

 {
   MONint noeud,noproc;
  
   while(list)
    {
     noeud = list->num;
     list = list->suivant;
     noproc = procnode1[noeud];

     noeud = fils[noeud];
     if (noeud) 
     {
      while (noeud>0)
       { procnode1[noeud] = noproc;
         noeud = fils[noeud];
       }
      noeud = -noeud;
      while (noeud >0 )
       {
        dmumps_sousarbreaux(procnode1,noeud,noproc,frere,fils);
        noeud = frere[noeud];
       }
     } 
    }
}


/**************************************/
/* construction de la strate initiale */
/**************************************/
   void dmumps_strateinit(list,procnode1,n,nbproc,nbroot,ne,frere,fils,tcout,tabcout,
                   poids)
    MONint n,nbproc,nbroot;
    MONint *ne,*frere,*fils,*procnode1;
    double *tabcout, *tcout;
    MONint *poids;
    liste **list;
{
  double cout;
  MONint noeud,pere;
  MONint fin,i;
  liste *list1,*list2;

  /* printf("Entering strate init\n"); */
/*list2 = (liste *)malloc(sizeof(liste));*/

  list2 = *list;
  fin = FALSE;

/* rangement initial */
  if (nbroot>=nbproc)
    {
     dmumps_ranger(*list,n,nbproc,procnode1,tabcout,poids);
     if (dmumps_equilibrage(tabcout,nbproc)) fin=TRUE;
    }

/*  iteration du rangement */
  while(!fin)
  {
/* determination du noeud de plus gros cout */
   list2 = *list;
   cout  = list2->lcout;
   noeud = list2->num;
   /* printf("Node %d, Cost %lf \n",list2->num,list2->lcout); */
   list2 = list2->suivant;
   while(list2)
    {
    /* printf("Node %d, Cost %lf \n",list2->num,list2->lcout); */
     if (list2->lcout>cout)
        { cout  = list2->lcout;
          noeud = list2->num;
        }
     list2 = list2->suivant;
    }
/* remplacement du noeud de plus gros cout par ses fils */
   pere = noeud;
   noeud = fils[noeud];
   while(noeud>0) noeud = fils[noeud];
   if (noeud)
      {
       noeud = - noeud;
       list1 = dmumps_creerliste(tcout[noeud],noeud);
       /* printf("Noeud %d eclate en %d \n",pere,ne[pere]);*/
       for (i=2;i<(ne[pere]+1);i++)
        {
         noeud = frere[noeud];
         list1 = dmumps_ajouter_avant(tcout[noeud],noeud,list1);
        }
       *list=dmumps_remplacer(cout,*list,list1);
       dmumps_ranger(*list,n,nbproc,procnode1,tabcout,poids);
/* test de l'equilibrage de la charge */
       if (dmumps_equilibrage(tabcout,nbproc))
           {
              fin=TRUE ;
             /* printf("Equilibrage ok \n");*/
           }
      }
   else 
      {
/* si le noeud a remplacer est une feuille on fait un dernier */
/*   rangement et si il n'y a pas d'equilibrage   */
/*    on fait un equilibrage simple */
       fin=TRUE;
       dmumps_ranger(*list,n,nbproc,procnode1,tabcout,poids);
/* November 1999 : we keep the previous mapping 
       if (!dmumps_equilibrage(tabcout,nbproc)) 
        {
         dmumps_rangersimp(*list,n,nbproc,procnode1,tabcout,poids);
        }
*/
      }
  }

    dmumps_sousarbre(*list,procnode1,frere,fils) ;

}

/**********************************************************************/
  MONint dmumps_appartient(noeud,list)
   MONint noeud;
   liste *list;
{
  while(list)
   {
    if (list->num==noeud) return TRUE;
    list = list->suivant;
   }
  return FALSE;
}
/*************/
     
 MONint dmumps_perevalideaux(noeud,nstk,frere)
  MONint noeud;
  MONint *nstk,*frere;
{

 if (noeud<0) return TRUE;
 else
  {
   if (nstk[noeud]==0) return FALSE;
   else return dmumps_perevalideaux(frere[noeud],nstk,frere);
  }
}
/************/
/* dmumps_perevalide determine si on peut prendre un pere */
/*  un pere est valide si tous ses fils sont places */

 MONint dmumps_perevalide(noeud,nstk,frere,fils)
  MONint noeud;
  MONint *nstk,*frere,*fils;
{
  while(noeud>0) noeud = fils[noeud];
  noeud = -noeud;
  return dmumps_perevalideaux(noeud,nstk,frere);
}

/************************************/
/* construction des autres strates  */
/************************************/

   void dmumps_strate(n,procnode1,list,nbproc,ncout,frere,fils,tabcout,nstk,poids)
    MONint n,nbproc;
    double *ncout,*tabcout;
    MONint *frere,*fils,*procnode1;
    MONint *nstk;
    MONint *poids;
    liste **list;
 {
   liste *list1,*list2;
   MONint noeud,debut;

   debut = TRUE;
/* list2 = (liste *)malloc(sizeof(liste));*/
   list2 = *list;
   while(list2)
    {
     noeud = list2->num;
     list2 = list2->suivant;

     noeud = frere[noeud];
     while(noeud>0) noeud = frere[noeud];
     noeud = -noeud;

/* ajout du noeud a la strate suivante si son pere est valide */
     if (noeud)
     {
      if (dmumps_perevalide(noeud,nstk,frere,fils))
      {
       if (debut) 
        {
         list1=dmumps_creerliste(ncout[noeud],noeud);
         debut=FALSE;
        }
       else
       {
         if(!dmumps_appartient(noeud,list1)) 
           dmumps_ajouter(ncout[noeud],noeud,list1);
       }
      }
     }
    }
/* rangement de la strate sur les processeurs */
   if (!debut) dmumps_rangerstrat(list1,n,nbproc,procnode1,tabcout,poids) ;
   else 
    {  procnode1[0] = 1;
       list1 = NULL;
    }
   *list = list1;


}
    
   
/**********************************************************************/
/***********************************/
/*  affichage de l'arbre distribue */
/***********************************/

/**********************/
/* ecriture */
 void dmumps_aff_liste(list)
  liste *list;
{
  MONint nb,nb_vid,i;
 
  nb = list->lcout;
  nb_vid = list->num;
  if (nb!=nb_vid)
  {
   printf("> ");
   list = list->suivant;
   for (i=1;i<nb;i++)
    {
     if (list->num) printf("%d[%lf] ",list->num,list->lcout);
     else printf("# ");
     list = list->suivant;
    }
   if (list->num) printf("%d[%lf]\n",list->num,list->lcout);
   else printf("\n");
  }
 } 
/***********************/
/* creation de la liste */
 liste *dmumps_rep_liste(list,procnode,frere,fils,ne)
 MONint *procnode,*frere,*fils,*ne;
 liste *list;
{
 liste *list1,*list2;
 MONint nb,nb_vid,i;
 MONint noed,noeud,pere;
 MONint debut;

/* list2 = (liste *)malloc(sizeof(liste));*/
 
 list2 = list;
 list1 = dmumps_creerliste(0.0,1);
 nb = 0;
 nb_vid = 0;
 debut = TRUE;

 list2 = list2->suivant;
 while (list2)
  {
   pere = list2->num;
   list2 = list2->suivant;
   if (!debut) 
    {
     if(pere) {dmumps_ajouter(0.0,0,list1); nb++;nb_vid++;}
    }
   debut = FALSE;
   if(pere)
   {
    noeud = fils[pere]; 

    if(noeud)
     {
      if (noeud>0) 
        {
         dmumps_ajouter(procnode[noeud]+100.0,noeud,list1);
         nb++;
        }
      else
       {
        noeud = - noeud;
        noed  = noeud;
        noed = frere[noed];
        while(noed>0) noed = frere[noed];
        nb = nb + ne[-noed];
        for (i=1;i<=ne[-noed];i++)
         {
          dmumps_ajouter(1.0*procnode[noeud],noeud,list1);
          noeud = frere[noeud];
         }
       }
     }
   }
  }
 list1->lcout = nb; 
 list1->num = nb_vid; 
 return list1;
}

#if defined(DEBUG)
/***********************/
 void dmumps_represent_aux(pos,procnode,frere,fils,ne)
 MONint *procnode,*frere,*fils,*ne;
{
 liste *list;
 list = dmumps_creerliste(1.0,0);
 dmumps_ajouter(1.0*procnode[pos],pos,list);
 while(list->lcout)
  {
    dmumps_aff_liste(list);
    list = dmumps_rep_liste(list,procnode,frere,fils,ne);
  }
}
/***********************/
 void dmumps_represent(n,procnode,frere,fils,ne)
 MONint *procnode,*frere,*fils,*ne;
 MONint n;
{
 MONint i;
 for (i=1;i<=n;i++)
  {
   if (frere[i] == 0) 
     {
      dmumps_represent_aux(i,procnode,frere,fils,ne);
      printf("\n");
     }
  }
}
#endif
  
/**********************************************************************/

/*******************************************/
/* MAIN routine   ---  called from Fortran */
/* We duplicate here symbol names to avoid */
/* preprocessing depending on the machine  */
/*******************************************/


/* SUN/SGI interface to C */
void dmumps_distri_(pn,pnbproc,ne,na,nd,frere,fils,procnode,
             poids,ssarbre,nbsa,pk50)

    MONint *ne,*na,*nd,*frere,*fils,*poids,*nbsa,*pk50;
    MONint *pn,*pnbproc;
    MONint *procnode;
    MONint *ssarbre;
{
/* input data not modified */
/* *************************/
/* pn : matrix order */
/* pnbproc: number of processors */
/* pne,pna,pnd,pfrere,pfils: tableaux definissant l'arbre d'elimination */ 
/* poids contient la valuation des processeurs */


/* output data */
/* *************/
/* procnode (i), i=1,N  est le tid du processeur traitant le noeud i */
/* ssarbre contient les noeuds de la strate initiale, */
/*         les noeuds sont des racines de sous arbre , 
           la strate initiale sera utilisee par le solve */

/* local variables         */
/* *************************/
/* tid contient le rang des processeurs*/
/* procnode1 contient la distribution des noeuds/variables 
             (un noeud par variable) sur les processeurs 
             procnode1(i)=numero de processeur 
    tous les noeuds (non amalgames sont consideres --> i E [1,N]    */
/* tabcout contient la charge par processeur */
/* combien contient le nombre de variables eliminees par processeur */
/*         il sert pour les statistiques  */
/* nstk definit si un noeud/variable est place ou non */
/*      nstk[i]= 0 le noeud/variable i n'est pas place */
/* tcout[1..N] contient le nb de flops du noeud i et du sous arbre associe  */
/* ncout[1..N] contient le nb de flops du noeud i  */


/* Parametre/variable controlant le deroulement de l'algorithme */
/* *******************/
/* hazard est une constante 
           si hazard = 0 alors distribution cyclique no de variable sur 
                               processeurs 
           sinon  distribution equitable */


    MONint *procnode1, *nstk, *combien;
    double *tcout,*ncout,*tabcout;
    MONint n,nbproc,k50;
    MONint nbroot,noeud;
    MONint i,j,debut,nbs;
    liste *list,*list1;
   
    MONint hazard;


/* initialisation des tableaux */
   n = *pn;
   nbproc = *pnbproc;
   k50 = *pk50;

/* allocate a few dynamic arrays ( + 1 just in case, avoids some
checkings
*/

    procnode1 = ( MONint *)malloc( ( n + 1 )       * sizeof( MONint ) );
    tabcout   = (double *)malloc( ( nbproc + 1 )  * sizeof( double ) );
    tcout     = (double *)malloc( ( n + 1 )       * sizeof( double ) );
    ncout     = (double *)malloc( ( n + 1 )       * sizeof( double ) );
    nstk      = ( MONint *)malloc( ( n + 1 )       * sizeof( MONint ) );
    combien   = ( MONint *)malloc( ( nbproc + 1 )  * sizeof( MONint ) );


/* list1 = (liste *)malloc(sizeof(liste)); */

   list = NULL;


   for(i=1;i<=n;i++) nstk[i] = 0;
   for(i=0;i<=n-1;i++) ssarbre[i] = 0;
   for(i=1;i<=nbproc;i++) tabcout[i] = 0;

     ne--;
     na--;
     nd--;
     frere--;
     fils--;
     poids--;

    hazard = 1;
if (hazard!=0)
 { 
/****************************************/
/* distribution par strates */

   debut = TRUE;
   if (nbproc>1)
     {
/* calcul des couts des noeuds de l'arbre */
       dmumps_fcout(n,ne,nd,frere,fils,tcout,ncout,k50);

       if (na[n] < 0)
             nbroot = -na[n] - 1;
       else 
             nbroot = na[n];

/* creation de la liste initiale des noeuds racines */
       for (i= 1;i<(n+1);i++)
           if (frere[i] == 0) 
            {
              if (debut) {
                          list  = dmumps_creerliste(tcout[i],i);
                          list1 =list;
                          debut = FALSE;
                         } 
              else {
                   list1->suivant = dmumps_creerliste(tcout[i],i);
                   list1          = list1->suivant;
                   }
            } 

/* calcul de la strate initiale */
       dmumps_strateinit(&list,procnode1,n,nbproc,nbroot,ne,frere,fils,tcout,
                  tabcout,poids);

/* tableau ssarbre  */
       list1 = list;
       nbs=0;



#if defined (check)
printf("Subtree list\n");
#endif

       while(list1)
        {
         nstk[list1->num] = 1;
         ssarbre[nbs++] = list1->num;

/* INSERT HERE SUBTREE COSTS AND MAPPING */


#if defined (check)
printf("One subtree of total cost %lf mapped on %d \n",list1->lcout,procnode1[list1->num]);
#endif


         list1=list1->suivant;
        }

/*

#if defined(check)
printf("\n\n");
#endif

*/


       *nbsa = nbs;

/* calcul des autres strates */

       procnode1[0] = 0;
       while(procnode1[0]==0)
       {
         dmumps_strate(n,procnode1,&list,nbproc,ncout,frere,fils,tabcout,nstk,poids);
         list1 = list;
         while(list1)
          {
           noeud = fils[list1->num];
           while(noeud>0)
            {
             if (nstk[noeud]==0)
              {procnode1[noeud] = procnode1[list1->num];
               nstk[noeud] = 1;
              }
             noeud = fils[noeud];
            }
           nstk[list1->num] = 1;
           list1=list1->suivant;
          }
       }
     }
    else
     {
      nbs = 0;
      for (i=1;i<=n;i++)
       {
        procnode1[i] = 1;
/* Use zero subtrees instead of haveing all nodes being subtress
   in the case where there is only one processor
        ssarbre[nbs++] = i;  */
       }
       
       *nbsa = nbs;
      }
}
/****************************************/
/* distribution cyclique ( aleatoire )*/
else
{
 j = 1;
 for (i=1;i<=n;i++)
  {
   procnode1[i] = j++; 
   if (j==(nbproc+1)) j = 1;
  }
}
    for (i=n-1;i>=0;i--)
     {
/*    procnode[i] = tid[procnode1[i+1]]; */
      procnode[i] = procnode1[i+1];
     } 
    
/****************************************************************/
/* affichage de l'arbre avec sa disposition sur les processeurs */
/****************************************************************/
#if defined(DEBUG)
    for (i=1;i<=nbproc;i++) combien[i] = 0;
    for (i=1;i<=n;i++)
       { combien[procnode1[i]] = combien[procnode1[i]] + 1;}
    for (i=1;i<=nbproc;i++) printf("proc %d: %d noeuds\n",i,combien[i]);
#endif
#if defined(DEBUG)
    dmumps_represent(n,procnode1,frere,fils,ne);                     
#endif
/****************************************************************/

/* liberation de la liste list */

while (list)
    { list1 = list->suivant;
      free(list);
      list = list1;
    }

/* liberation des tableaux */
/*

ne, na, nd, fils, frere, poids are not freed since they are
passed as arguments from Fortran (before, were local copies)

free(tid);
free(ne);
free(na);
free(nd);
free(fils);
free(frere);
free(poids); */

free(procnode1);
free(tabcout);
free(combien);
free(tcout);
free(nstk);
free(ncout);
}


/* HP or IBM interface to C */
void dmumps_distri(pn,pnbproc,ne,na,nd,frere,fils,procnode,poids,ssarbre,nbsa,pk50)
    MONint *ne,*na,*nd,*frere,*fils,*poids,*nbsa,*pk50;
    MONint *pn,*pnbproc;
    MONint *procnode;
    MONint *ssarbre;
{
    dmumps_distri_(pn,pnbproc,ne,na,nd,frere,fils,procnode,
             poids,ssarbre,nbsa,pk50);
}

/* CRAY interface to C */
void DMUMPS_DISTRI(pn,pnbproc,ne,na,nd,frere,fils,procnode,poids,ssarbre,nbsa,pk50)
    MONint *ne,*na,*nd,*frere,*fils,*poids,*nbsa,*pk50;
    MONint *pn,*pnbproc;
    MONint *procnode;
    MONint *ssarbre;
{
    dmumps_distri_(pn,pnbproc,ne,na,nd,frere,fils,procnode,
             poids,ssarbre,nbsa,pk50);
}

/* ALPHA interface to C */
void dmumps_distri__(pn,pnbproc,ne,na,nd,frere,fils,procnode,poids,ssarbre,nbsa,pk50)
    MONint *ne,*na,*nd,*frere,*fils,*poids,*nbsa,*pk50;
    MONint *pn,*pnbproc;
    MONint *procnode;
    MONint *ssarbre;
{
    dmumps_distri_(pn,pnbproc,ne,na,nd,frere,fils,procnode,
             poids,ssarbre,nbsa,pk50);
}

