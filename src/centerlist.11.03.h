/*****************************************/
/*  centerlist.h contains routines and   */
/*  definitions for centerlist           */
/*  started 3/3/01                       */
/* modified 11.03 Gary Benson            */
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
						
	4) Each center has a linked list which holds information
	about the runs of matches with that center.
							
	5) Note that the true center times 2 is used
	to avoid indices with 1/2 in them.
	First center is at 2*Tuplesize[g]+1;
	Last center is at 2*Length+1;
*/
/**************************************************/

#include <assert.h>

/*
int Centerlistsize[MAXTUPLESIZES+1];
int Centerlistmaxdistance[MAXTUPLESIZES+1];
int Centerlistmatcharraysize[MAXTUPLESIZES+1];
int Firstcenter,Lastcenter;
*/

/* new definition for centerentry and centerlist 4-14-01 */
/* centerentry becomes matchentry */

/* modified 11.03 Gary Benson                       */
struct matchentry3{
	int l_index;           /* first index (low) of match */
	int h_index;           /* last index (high) of match */
	double num_matches;            /* number of matches in this entry */
	double tot_matches;            /* number of matches from start of list including this entry */
	struct matchentry3 *next;
};


/* modified 11.03 Gary Benson                       */
struct centerlist3{
	int index;                   /* index of center times 2 */ 
	int k_run_sums_criteria;          /* criteria stored in list */
	int delta;                        /* this is +- delta */
	int numentries;                   /* number of entries in linked list of matches */
	struct matchentry3 *match;        /* pointer to first entry in match list */
	struct matchentry3 *end_of_list;  /* pointer to last entry in match list */
  struct matchentry3 *interval_lo[MAXNUMINTERVALS]; /* pointers to first entry
																		at low end of each interval */
	//int linked;                       /* a flag */
    int linkdown,                /* a true center value or -2, next entry in a linked list */
			     linkup;                  /* a true center value or Maxcenter[g]+1, last entry in a linked list */
    //char waitingToAlign;
} *Center3[MAXTUPLESIZES+1];

/* note that there must be a centerlist for each tuple size
because the same center will be used by different tuple
sizes
*/



/****************************************************/
/* new_centerlist creates memory for centerlist and */
/* linked list of matches                           */
/* the same Center will be on a list for each       */         
/* tuple size                                       */
/* modified 11.03 Gary Benson                       */
/****************************************************/

struct centerlist3 *new_centerlist3(int g)
/* g is in [1..NTS] */
{
  struct centerlist3 *objptr=(struct centerlist3 *)
		scalloc(Centerlistsize[g],sizeof(struct centerlist3));

  memory_stats_print("\n new_centerlist3: requesting ", Centerlistsize[g] * sizeof(struct centerlist3) );

  return(objptr);
}

/****************************************************/
/* clear_center zeros out information in the given  */
/* center                                           */
/* modified 11.03 Gary Benson                       */
/****************************************************/
void clear_center3(struct centerlist3 *objptr, int k)
/* k is in [0..Centerlistsize[g]-1] */
{
  int h;
  struct matchentry3 *temp;
  
  objptr[k].index=-2; /* -2 means no index */
  objptr[k].numentries=0;
  objptr[k].end_of_list=NULL;
  //objptr[k].waitingToAlign = 0;
  temp=objptr[k].match;
  while(temp!=NULL)
  {
    objptr[k].match=temp->next;
    sfree(temp);
    temp=objptr[k].match;
  }
  for(h=0;h<MAXNUMINTERVALS;h++)
    objptr[k].interval_lo[h]=NULL;
}


/****************************************************/
/* clear_centerlist zeros out information in the    */
/* centerlist                                       */
/* modified 11.03 Gary Benson                       */
/****************************************************/

void clear_centerlist3(struct centerlist3 *objptr,int g)
/* g is in [1..NTS] */
{
  int k;
  for(k=0;k<Centerlistsize[g];k++)
		clear_center3(objptr,k);
}
/****************************************************/
/* ClearCenterlists clearts all centerlists         */
/* modified 11.03 Gary Benson                       */
/****************************************************/

void ClearCenterlists3(void)
{
  int g;
  for(g=1;g<=NTS;g++)
  {      
    /* zero out the values in Center */
    clear_centerlist3(Center3[g],g);
  }
	
}
/****************************************************/
/* InitCenterlists initializes the centerlists      */
/* modified 11.03 Gary Benson                       */
/****************************************************/

void InitCenterlists3(void)
{
  int g;
  for(g=1;g<=NTS;g++)
  {
    /* define the size of the circular centerlist */ 
    Centerlistsize[g]=Tuplemaxdistance[g]-Tuplesize[g]+1+Delta[g];
		
    /* allocate memory for Center */   
    Center3[g]=new_centerlist3(g);

    if (Center3[g]==NULL ) {
		debugerror("ERROR(InitCenterlists3): cannot allocate a new center. Aborting!");
		exit(-1);
	}

    
    /* zero out the values in Center */
    clear_centerlist3(Center3[g],g);
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
/* modified 11.03 Gary Benson                       */
/*****************************************************/

void PrintCenterMatchlist3(int g, int c)
{
	
	struct matchentry3 *temp;
  
  temp=Center3[g][c].match;
  debugerror("\n          Center3[%3d][%3d]: (h_index,num_matches,tot_matches)",g,c);
  while(temp!=NULL)
  {
		
    debugerror("(%3d,%3d,%3d) ",
      temp->h_index, temp->num_matches, temp->tot_matches);
  }
}







/***************************************************************/
void shift_interval_pointers_and_discard_old_matches3(int location, int g, int c)
/* modified 11.03 Gary Benson        */
/* location is index of tuple        */
/* g is in [1..NTS]                  */
/* c is index in circular centerlist */
{
  struct centerlist3 *centergc;
  struct matchentry3 *temp, *lopointer;
	int h_temp,lowend,h;

  centergc=&Center3[g][c];
	
  /* assume that interval pointers can only move to the right along the list */
	/* but, a pointer might move past the last match (when testing criteria) */
		
  for(h=1;h<=Intervals[g].num_intervals;h++)
	{
		if((lopointer=centergc->interval_lo[h])!=NULL) /* check first for null interval_lo pointer */
		{
			lowend = location-Intervals[g].size[h]+1;
			//h_temp=centergc->interval_lo[h]->h_index;
			h_temp=lopointer->h_index;
			while(h_temp<lowend)
			{			
				//centergc->interval_lo[h]=centergc->interval_lo[h]->next;
				centergc->interval_lo[h]=lopointer=lopointer->next;
				//if(centergc->interval_lo[h]==NULL) break; 
				if(lopointer==NULL) break; 
				//{fprintf(stderr,"\n-16 center[%d][%d], h=%d, location=%d, , h_temp=%d, lowend=%d",
				//g,c,h,location,h_temp,lowend); exit(0);}
				h_temp=lopointer->h_index;
			}
		}
	}
  /* remove any match elements from linked list that fall before the low end of the 
	largest interval */
	/* note that this may remove all the matches */

	temp=centergc->match;
	while (temp!=centergc->interval_lo[Intervals[g].num_intervals]) /* which could be null */
	{
		centergc->match=temp->next;
		centergc->numentries--;
		if(centergc->numentries==0) 
		{
			centergc->end_of_list=NULL;
			for(h=1;h<Intervals[g].num_intervals;h++)
				if(centergc->interval_lo[h]!=NULL) 
					fprintf(stderr,"\nerror in shift intervals, empty list but interval_lo pointers not null");
		}
		if(centergc->numentries<0)fprintf(stderr,"\nBig error numentries goes negative in shift intervals");
		//fprintf(stderr,"\n removed entry, h=%d, h_index=%d, i=%d g=%d",h,temp->h_index,location,g);
		sfree(temp);
		temp=centergc->match;
		//if(temp==NULL) exit(0);
	}
}



void add_tuple_match_to_Center3(int location, int size, int d, int g, int c, int rcloc) 
/* modified 11.03 Gary Benson                   */
/* location is index of tuple                   */
/* size is size of tuple                        */
/* d is the true center                         */
/* g is in [1..NTS]                             */
/* c is index in circular centerlist            */


// added by Gelfand on 11/18/2004 Need for GTmatch
/* rcloc is the leftmost char of the RCtuple    */

{
	struct centerlist3 *centergc;
	struct matchentry3 *lastmatch, *temp;
	int h;
    int ph, phCount, loc;
    double gtPairs = 0.0;

	centergc=&Center3[g][c];

	/* test if this list is currently used for entries at center d    */
	/* if not, initialize the information, updaterunningminmaxcenters */
	/* has already cleared the entry                                  */

	if(centergc->index!=d) centergc->index=d;

	/* obtain pointer to last match */

	lastmatch=centergc->end_of_list;

	/* test if this match is a continuation of the preceding match */
	if((centergc->numentries!=0)
		&&(lastmatch->h_index==location-1)) 
		/* there are entries and this is an adjacent tuple, */
		/* just add on to last entry */
	{
		//debugmessage("\n        add to existing entry");
		lastmatch->h_index++;
		//lastmatch->tot_matches++;
		//lastmatch->num_matches++;

        /* if gtmatch, do not give a full score to GT match */
        if (rcloc) {
          
            if ((Sequence[location]=='G' && Sequence[rcloc]=='T') || (Sequence[location]=='T' && Sequence[rcloc]=='G')) 
                                gtPairs++;

            if (gtPairs!=0.0) {
                lastmatch->tot_matches+=GTDetectMatch;
                lastmatch->num_matches+=GTDetectMatch;
            } else {
                lastmatch->tot_matches++;
                lastmatch->num_matches++;
            }

        } else {
                lastmatch->tot_matches++;
                lastmatch->num_matches++;
        }

	}
	else  /* need a new entry here */
	{

		//debugmessage("\n        new entry");
		temp=(struct matchentry3 *)scalloc(1,sizeof(struct matchentry3));
/*	    if (temp==NULL) {
		    debugerror("ERROR(add_tuple_match_to_Center3): cannot allocate %d bytes of memory. Aborting!",sizeof(struct matchentry3));
		    exit(0);
	    }*/
        memory_stats_print("\n add_tuple_match_to_Center3: requesting ", sizeof(struct matchentry3) );


        /* if gtmatch, find the number of GT pairs between the tuples and give partial scores to these pairs */
        if (rcloc) {
            for (loc=location,ph=rcloc,phCount=0; phCount<size; loc--,ph++,phCount++ ) {
                if ((Sequence[loc]=='G' && Sequence[ph]=='T') || (Sequence[loc]=='T' && Sequence[ph]=='G')) 
                    gtPairs++; 

              //  fprintf(stderr,"\n loc: %c  ph: %c", Sequence[loc], Sequence[ph]);
            }
            //fprintf(stderr,"\n"); fflush(stderr);
        }


        /* debug 
        if (gtPairs>0) {
            fprintf(stderr,"\nGotcha!!!\n");
            exit(1);
        }*/


		temp->h_index=location;
		temp->l_index=location-size+1;
		temp->num_matches = (double)size - gtPairs + gtPairs*GTDetectMatch;
		//temp->next=NULL; /* shouldn't need this with scalloc */

		if(centergc->numentries!=0) /* list is not empty, append to list */ 
		{
			if(centergc->end_of_list==NULL) 
				fprintf(stderr,"\nerror in add tuple 3, numentries!=0 but end of list is null");
			temp->tot_matches=centergc->end_of_list->tot_matches+temp->num_matches;
			centergc->end_of_list->next=temp;
			centergc->end_of_list=temp;
			for(h=1;h<=Intervals[g].num_intervals;h++)
			{
				/* initialize all interval pointers that currently point to null because the 
				last match on the list lies outside these intervals */
				if(centergc->interval_lo[h]==NULL)
					centergc->interval_lo[h]=temp;
				//fprintf(stderr,"\n interval pointer %d h_index=%d, i=%d, g=%d, c=%d",h,centergc->interval_lo[h]->h_index,location,g,c);
			}

		}
		else /* list is empty, temp becomes first element */
		{
			temp->tot_matches=temp->num_matches;
			centergc->match=centergc->end_of_list=temp;
			/* initialize all interval pointers for the first element in list */
			for(h=1;h<=Intervals[g].num_intervals;h++)
			{
				centergc->interval_lo[h]=temp;
				//fprintf(stderr,"\n interval pointer %d h_index=%d, i=%d, g=%d, c=%d",h,centergc->interval_lo[h]->h_index,location,g,c);
			}
		}
		centergc->numentries++;
	}

	/* shift interval pointers */
	//printf("\n\nfrom add match");
	shift_interval_pointers_and_discard_old_matches3(location,g,c);
	//printf("\nto add match");


}

/***************************************************************/
void Update_Runningmaxmincenters3(int i, int g)
/* modified 11.03 Gary Benson            */
/* moves the maxcenter and mincenter as the index of the sequence moves */
/* i is index in sequence */
/* g is in 1..NTS */
{
	int ormc,cormc;


	/* define the new Runningmaxcenter and Runningmincenter */
	
	Runningmaxcenter3[g]=2*i-2*Tuplesize[g]+1;  
	Runningmincenter3[g]=2*i-Tuplemaxdistance[g]-Tuplesize[g]+1;  
	//printf("\n    g=%3d  Runningmaxcenter=%3d  Runningmincenter=%3d  Oldrunningmaxcenter=%3d  Oldrunningmincenter=%3d",
	//	g,Runningmaxcenter[g],Runningmincenter[g],Oldrunningmaxcenter[g],Oldrunningmincenter[g]);
	
	/* clear centers Oldruningmincenter to Runningmincenter-1 */
	if(Runningmincenter3[g]>=Mincenter3[g]) /* Mincenter is the smallest possible center */
	{
		/* ormc is a true center */
		ormc=Oldrunningmincenter3[g];
		while(ormc!=Runningmincenter3[g])
		{
			/* cormc is an index in the centerlist */
			cormc=ormc%Centerlistsize[g];
			
			/* remove information stored in center */
			
			if(Center3[g][cormc].index!=-2) /* if it is -2, then this center is already empty */
				clear_center3(Center3[g],cormc);
			ormc++;
		}
		Oldrunningmincenter3[g]=ormc;
	}
	
	/* no clearing necessary for runningmaxcenter */

	Oldrunningmaxcenter3[g]=Runningmaxcenter3[g];
}

void Clear_MaxMincenters3(void) {

   int g;

  for(g=1;g<=NTS;g++)
  {
    Maxcenter3[g]=2*Length-2*Tuplesize[g]+1;
    Mincenter3[g]=2*Tuplesize[g]+1;
    Runningmaxcenter3[g]=0;
    Runningmincenter3[g]=0;
    Oldrunningmaxcenter3[g]=0;
    Oldrunningmincenter3[g]=0;
	/*debugmessage("\n\nInit_MaxMincenters\n  g  Maxcenter  Mincenter  Runningmax  Runningmin  Oldrunmax  Oldrunmin");
    debugmessage("\n%3d %6d  %9d  %9d   %9d   %9d",g,Maxcenter[g],Mincenter[g],
      Runningmaxcenter[g],Runningmincenter[g],Oldrunningmincenter[g]);
	*/
  }

}

void Init_MaxMincenters3(void)
{

  Maxcenter3=(int *)scalloc(MAXTUPLESIZES+1,sizeof(int));
  Mincenter3=(int *)scalloc(MAXTUPLESIZES+1,sizeof(int));
  Runningmaxcenter3=(int *)scalloc(MAXTUPLESIZES+1,sizeof(int));
  Runningmincenter3=(int *)scalloc(MAXTUPLESIZES+1,sizeof(int));
  Oldrunningmaxcenter3=(int *)scalloc(MAXTUPLESIZES+1,sizeof(int));
  Oldrunningmincenter3=(int *)scalloc(MAXTUPLESIZES+1,sizeof(int));
  
  memory_stats_print("\n Init_MaxMincenters3(Maxcenter): requesting", (MAXTUPLESIZES+1) * sizeof(int) );
  memory_stats_print("\n Init_MaxMincenters3(Mincenter): requesting", (MAXTUPLESIZES+1) * sizeof(int) );
  memory_stats_print("\n Init_MaxMincenters3(Runningmaxcenter3): requesting", (MAXTUPLESIZES+1) * sizeof(int) );
  memory_stats_print("\n Init_MaxMincenters3(Runningmincenter3): requesting", (MAXTUPLESIZES+1) * sizeof(int) );
  memory_stats_print("\n Init_MaxMincenters3(Oldrunningmaxcenter3): requesting", (MAXTUPLESIZES+1) * sizeof(int) );
  memory_stats_print("\n Init_MaxMincenters3(Oldrunningmincenter3): requesting", (MAXTUPLESIZES+1) * sizeof(int) );

  /* added by Gelfand on 10/30/03 
  if ( NULL == Maxcenter3 || NULL == Mincenter3 || NULL == Runningmaxcenter3 || NULL == Runningmincenter3 || NULL == Oldrunningmaxcenter3 || NULL == Oldrunningmincenter3 ) 
    { debugerror("\nMemory error in Init_MaxMincenters3. Aborting!"); }
*/

  Clear_MaxMincenters3();

}

/**********************************************************/
double intervalmatches3(int i, int g, int c, int h)
/* created by Gary Benson 11/03 */
/* i is the location is the last examined location in the sequence */
/* g is the tuple number     */
/* c is the centerlist index */
/* h is the interval number  */

{
	double hang_over,end_tot,start_tot,start_nummat;
	struct centerlist3 *centergc;
	struct matchentry3 *lopointer;

	/* following code finds number of matches in the interval ending at i */
  centergc=&Center3[g][c];
	if(centergc->end_of_list==NULL) return(0);
	end_tot=centergc->end_of_list->tot_matches;
	lopointer=centergc->interval_lo[h];
	if(lopointer==NULL) return(0);
	start_tot=lopointer->tot_matches;
	start_nummat=lopointer->num_matches;
	/* comment out next line so routine counts only matches that occur within
	the interval instead of all matches in an entry that overlaps the low end 
	of the interval */
	//return(end_tot-start_tot+start_nummat);
	/* comment out next three lines to make this routine the same as 
	the version that uses circular lists for the matche */
	hang_over=i-Intervals[g].size[h]+1-lopointer->l_index;
	if(hang_over<0) hang_over=0;
	return(end_tot-start_tot+start_nummat-hang_over);
	//return(end_tot-start_tot+start_nummat); /* test without the hangover, 3.02. gelfand */
}

/**********************************************************/
int TestIntervalCriteria34(int g, int t, int h, int i)
/* modified 11.03 Gary Benson            */
/* g is tuplesize index */
/* t is true center  */
/* h is the index of the interval */
/* i is the location of last matched character in sequence */
{
	
	int s,c,delta,fl,flc,fh,fhc,fc;
	struct intervallist *intervalsg;
    double t_matches,sum, m;
	
	s=Centerlistsize[g];
	c=t%s;
	intervalsg=&Intervals[g];
	delta=intervalsg->delta[h];
	
	//printf("\nTestIntervalCriteria34(%d,%d,%d,%d)",g,t,h,i);
	
	/* get matches for t */
	t_matches=intervalmatches3(i,g,c,h);
	//printf("\n  t_matches(%d)=%d",c,t_matches);

	/* get matches in delta range below t */
	fl=max(t-2*delta,Mincenter3[g]);
	flc=fl%s;
	
	sum=0;
	fc=flc;
	while(fc!=c)
	{
		if(Center3[g][fc].numentries>0)
		{


			//printf("\n\nfrom test criteria");
			shift_interval_pointers_and_discard_old_matches3(i,g,fc);


			//printf("\nto test criteria");
			m=intervalmatches3(i,g,fc,h);
			//printf("\n  m=%d, (%d,%d,%d,%d)",m,i,g,fc,h);

			if(m>t_matches) return(0);  /* t is not best center */
			sum+=m;
		}
		//fc=(++fc%s);
		fc++;
		fc=fc%s;
	}
	
	sum+=t_matches;
	
	/* get matches in delta range above t */
	fh=min(t+2*delta,Maxcenter3[g]);
	fhc=fh%s;
	fc=(c+1)%s;
	while(fc!=fhc)
	{
		if(Center3[g][fc].numentries>0)
		{
			//printf("\n\nfrom test criteria");
			shift_interval_pointers_and_discard_old_matches3(i,g,fc);
			//printf("\nto test criteria");
			m=intervalmatches3(i,g,fc,h);
			//printf("\n  m=%d, (%d,%d,%d,%d)",m,i,g,fc,h);

			if(m>t_matches) return(0);  /* t is not best center */
			sum+=m;
		}
		//fc=(++fc%s);
		fc++;
		fc=fc%s;
	}
	//printf("\ntotal=%d, test=%d",sum,intervalsg->RNKP[h]);

	/* if sum is greater than rnkp criteria, return true */
	if(sum>=intervalsg->RNKP[h]) 
	{
		/* following 7 lines added by Gary Benson, 12/04/02 for testing */
		TestMaxIntervalSum.tupleindex=g;
		TestMaxIntervalSum.intervalindex=h;
		TestMaxIntervalSum.matches_in_interval=sum;
		TestMaxIntervalSum.lowcenter=fl;
		TestMaxIntervalSum.highcenter=fh;
		TestMaxIntervalSum.testcenter=t;
       

		return(1);
	}
	else return(0);
}
 



