/*****************************************/
/*  centerlist.h contains routines and   */
/*  definitions for centerlist           */
/*  started 3/3/01                       */
/*****************************************/

/*****************************************/
/*   Important centerlist quantities     */
/*

1) There is  one centerlist per tuple size.
   These are designated Center[g] g=1..NTS

2) Each centerlist contains a fixed number
   of centers.  This is a circular list
   and the positions are accessed by a
   mod function on the true center index.
   Each center occupies a
   center list only temporarily.
   The locations are designated
   Center[g][k] k=0..Centerlistsize[g]-1

   The Centerlistsize is designed so that
   when a the closest possible tuple match
   (i.e. two adjacent inverted repeat tuples)
   defines a new center, then the oldest center
   that formerly occupied the position of the
   new center in the centerlist has been
   processed or has been available for the
   delta test of a processed center.

   Centerlistsize[g]=Tuplemaxdistance[g]-Tuplesize[g]+1+Delta[g];

   where Delta[g] is the random walk parameter
   for Tuplemaxdistance[g]. 

3) The maximum distance between matching tuples
   for a center is Tuplemaxdistance[g].

4) Each center has an array which holds information
   about the runs of matches with that center.
   The size of this array is Centerlistmatcharraysize[g]=
   Tuplemaxdistance[g]/(2*(Tuplesize[g]+1))+3
   which leaves one or two extra cells than needed.

5) Note that the true center times 2 is used
   to avoid indices with 1/2 in them.
   First center is at 2*Tuplesize[g]+1;
   Last center is at 2*Length+1;
                                                  */
/**************************************************/

#include <assert.h>


int Centerlistsize[MAXTUPLESIZES+1];
int Centerlistmaxdistance[MAXTUPLESIZES+1];
int Centerlistmatcharraysize[MAXTUPLESIZES+1];
int Firstcenter,Lastcenter;


/* new definition for centerentry and centerlist 4-14-01 */
/* centerentry becomes matchentry */

struct matchentry{
	int l_index;           /* first index (low) of match */
	int h_index;           /* last index (high) of match */
	int num_matches;            /* number of matches in this entry */
	int tot_matches;            /* number of matches from start of list including this entry */
};

struct centerlist{
    int index;                   /* index of center times 2 */ 
    int k_run_sums_criteria;          /* criteria stored in list */
    int delta;                        /* this is +- delta */
    int numentries;                   /* number of entries in linked list of matches */
    int nummatches;                   /* number of matches in linked list of matches */
		int intervalmatches[MAXNUMINTERVALS]; /* matches in each interval size from the last match */
    int end_of_list;                  /* index of last entry in match list */
    int interval_lo[MAXNUMINTERVALS]; /* index of match entry at low end of interval */
    int interval_hi[MAXNUMINTERVALS]; /* index of match entry at high end of interval */
    int temp_interval_sum[MAXNUMINTERVALS]; /* holds a sum of matches from an interval matches
                                           calculation for use shortly after computing */
    struct matchentry *match;         /* pointer to array of match entries */
    int linked;                      /* a flag */
    int linkdown,               /* a true center value or -2, next entry in a linked list */
             linkup;                 /* a true center value or Maxcenter[g]+1, last entry in a linked list */
} *Center[MAXTUPLESIZES+1];

/* note that there must be a centerlist for each tuple size
   because the same center will be used by different tuple
   sizes
*/

//struct centerentry{
//    int location;  /* index of match */
//    int size;      /* size of match ending at index */
//    int total;     /* total matches from beginning of list including present entry */
//};

//struct centerlist{
//    int index; /* specifies actual center index in sequence */ 
//    int
//         k_run_sums_criteria, /* criteria stored in list */
//         waiting_time_criteria,	
//         lo_d_range,       /* this is +- delta */
//         hi_d_range;
//    int numentries,      /* number of entries in linked list of matches */
//        nummatches;      /* number of matches in linked list of matches */
//    int highindex;  /* last position in list of match entries */
//    int linked;          /* a flag */
//    int linkdown,   /* ??? next entry in a linked list */
//             linkup;     /* ??? last entry in a linked list */
//    struct centerentry *entry; /* linked list of matches */
//} *Center[MAXTUPLESIZES+1];

/* note that there must be a centerlist for each tuple size
   because the same center will be used by different tuple
   sizes
*/


/****************************************************/
/* new_centerlist creates memory for centerlist and */
/* linked list of matches                           */
/* the same Center will be on a list for each       */         
/* tuple size                                       */
/****************************************************/
//extern void memory_stats_print( char * message, int size );

struct centerlist *new_centerlist(int g)
/* g is in [1..NTS] */
{
  int k;
  struct centerlist *objptr=(struct centerlist *)
	scalloc(Centerlistsize[g],sizeof(struct centerlist));

  /* added by Gelfand on 10/30/03 */
  memory_stats_print("\n new_centerlist: requesting", Centerlistsize[g]*sizeof(struct centerlist) );
  //if ( NULL == objptr ) { debugerror("\nMemory error. Could not allocate centerlist of size %d bytes. Aborting!", Centerlistsize[g]*sizeof(struct centerlist)); exit(0); }

  for(k=0;k<Centerlistsize[g];k++)   /* zero entry of centerlist is not used */
  {
	 objptr[k].match=(struct matchentry *)
	  scalloc(Centerlistmatcharraysize[g]+1,sizeof(struct matchentry));

     /* added by Gelfand on 10/30/03 */
     memory_stats_print("\n new_centerlist: requesting (for matches)", (Centerlistmatcharraysize[g]+1)*sizeof(struct matchentry) );
     //if ( NULL == objptr[k].match ) { debugerror("\nMemory error. Could not allocate array of matches of size %d bytes. Aborting!", (Centerlistmatcharraysize[g]+1)*sizeof(struct matchentry)); exit(0); }
  }
  return(objptr);
}


/****************************************************/
/* clear_center zeros out information in the given  */
/* center                                           */
/****************************************************/
void clear_center(struct centerlist *objptr, int k, int g)
/* k is in [0..Centerlistsize[g]-1] */
{
  int h;

  	objptr[k].index=-2; /* -2 means no index */
	  objptr[k].numentries=0;
	  objptr[k].nummatches=0;
    objptr[k].end_of_list=-1;
    for(h=0;h<MAXNUMINTERVALS;h++)
			objptr[k].intervalmatches[h]=0;
}


/****************************************************/
/* clear_centerlist zeros out information in the    */
/* centerlist                                       */
/****************************************************/

void clear_centerlist(struct centerlist *objptr,int g)
/* g is in [1..NTS] */
{
  int k;
  for(k=0;k<Centerlistsize[g];k++)
   clear_center(objptr,k,g);
}
/****************************************************/
/* ClearCenterlists clearts all centerlists         */
/****************************************************/

void ClearCenterlists(void)
{
  int g;
  for(g=1;g<=NTS;g++)
  {      
    /* zero out the values in Center */
    clear_centerlist(Center[g],g);
  }

}
/****************************************************/
/* InitCenterlists initializes the centerlists      */
/****************************************************/

void InitCenterlists(void)
{
  int g;
  for(g=1;g<=NTS;g++)
  {  
	   /* define the size of the circular centerlist */ 
     Centerlistsize[g]=Tuplemaxdistance[g]-Tuplesize[g]+1+Delta[g];

     /* define the number of match entries that can be stored in a center */
     Centerlistmatcharraysize[g]=
       Tuplemaxdistance[g]/(2*(Tuplesize[g]+1))+3;

    /* allocate memory for Center */   
     Center[g]=new_centerlist(g);
    
    /* zero out the values in Center */
    clear_centerlist(Center[g],g);
  }
  /*
     debugmessage("\n\nCenterlists\n  tupsize index   tuple size   centerlistsize   matcharraysize");
  for(g=1;g<=NTS;g++)
  {
     debugmessage("\n      %d              %d              %4d              %3d",
       g,Tuplesize[g],Centerlistsize[g],Centerlistmatcharraysize[g]);
  }
  */
}


/*****************************************************/
/**** print CenterMatchlist **************************/
/* used for testing **********************************/
/*****************************************************/

void PrintCenterMatchlist(int g, int c)
{
	
  int k,end_of_list;
  
  end_of_list=Center[g][c].end_of_list;
  fprintf(stderr,"\n          Center[%3d][%3d]: (h_index,num_matches,tot_matches)",g,c);
  for(k=0;k<=end_of_list;k++)
  {

    fprintf(stderr,"(%3d,%3d,%3d) ",
      Center[g][c].match[k].h_index,Center[g][c].match[k].num_matches,
      Center[g][c].match[k].tot_matches);
  }
}

